#!/usr/bin/env python3
"""
Parse a RepeatMasker .out file and keep only rows intersecting HERV locus windows.

    python parse_rmsk_out.py --out-file hs1rep/hs1.repeatMasker.out.gz \
        --assembly t2t --db herv_db/herv_catalog.db --parquet rmsk_t2t_windows.parquet

COORDINATE CONVENTION -- the whole point of this module.

RepeatMasker .out is 1-BASED INCLUSIVE on both ends. The catalog is 0-BASED
HALF-OPEN (locus_coord.coord_convention says so explicitly, and v1.1 of the db
exists specifically because an earlier export was 1 bp high). So:

    start_0 = out_begin - 1
    end_0   = out_end            # inclusive end == half-open end

Getting this wrong shifts every repeat 1 bp right, which is invisible in a
graphic and wrong in every overlap calculation. Round-tripping is asserted in
self_test().

The 15 fields are positional, NOT tab-delimited -- the file is column-aligned
with variable whitespace, and field 9 ('+' or 'C') means strand where 'C' is the
reverse complement. Repeat names can contain spaces in rare cases (e.g.
'tRNA-Leu-TTA(m)' is fine, but some libraries emit parenthesised composites), so
we split on whitespace and validate the field count per row rather than trusting
a fixed split.
"""
from __future__ import annotations
import argparse, gzip, os, sqlite3, sys
from collections import defaultdict

import numpy as np
import pandas as pd

PAD = 1000          # window padding; matches PAD in detail.js and build_dashboard.py
N_FIELDS = 15


def _open(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def load_windows(db: str, assembly: str, pad: int = PAD) -> pd.DataFrame:
    """Locus windows for one assembly, padded. Returns 0-based half-open bounds."""
    con = sqlite3.connect(db)
    co = pd.read_sql("SELECT locus_uid,chrom,start,end,strand FROM locus_coord "
                     "WHERE assembly=?", con, params=[assembly])
    con.close()
    if co.empty:
        sys.exit(f"no locus_coord rows for assembly={assembly!r}")
    co["w0"] = (co.start - pad).clip(lower=0)
    co["w1"] = co.end + pad
    return co


def _index_windows(win: pd.DataFrame):
    """Per-chromosome sorted (w0, w1, uid) plus a w0 array for bisect."""
    by = defaultdict(list)
    for r in win.itertuples():
        by[r.chrom].append((r.w0, r.w1, r.locus_uid))
    out = {}
    for c, v in by.items():
        v.sort()
        out[c] = (np.array([a[0] for a in v]), np.array([a[1] for a in v]),
                  [a[2] for a in v], np.maximum.accumulate(np.array([a[1] for a in v])))
    return out


def parse(out_file: str, win: pd.DataFrame, keep_all_classes: bool = True):
    """Stream the .out file, emitting one row per (repeat, overlapping locus) pair.

    A repeat overlapping two nearby loci is emitted twice, once per locus -- the
    table is a locus-to-repeat relation, not a repeat list, so the dashboard can
    fetch by locus without an interval search in the browser.
    """
    idx = _index_windows(win)
    rows = []
    n_lines = n_kept = n_pairs = 0
    skipped_chrom = defaultdict(int)
    malformed = 0

    with _open(out_file) as fh:
        for ln in fh:
            f = ln.split()
            if len(f) != N_FIELDS:
                # header lines (3) and the occasional '*' continuation marker
                if f and f[0].isdigit():
                    malformed += 1
                continue
            n_lines += 1
            chrom = f[4]
            ent = idx.get(chrom)
            if ent is None:
                skipped_chrom[chrom] += 1
                continue
            w0a, w1a, uids, maxend = ent
            # 1-based inclusive -> 0-based half-open
            s = int(f[5]) - 1
            e = int(f[6])
            # candidate windows: those starting at or before e; scan back while
            # the running max end can still reach s
            hi = int(np.searchsorted(w0a, e, side="right"))
            if hi == 0 or maxend[hi - 1] <= s:
                continue
            hit = False
            for i in range(hi - 1, -1, -1):
                if maxend[i] <= s:
                    break
                if w1a[i] > s and w0a[i] < e:
                    rows.append((uids[i], chrom, s, e,
                                 "-" if f[8] == "C" else "+",
                                 f[9], f[10], int(f[0]), float(f[1]),
                                 max(s, w0a[i]), min(e, w1a[i])))
                    n_pairs += 1
                    hit = True
            n_kept += hit

    df = pd.DataFrame(rows, columns=["locus_uid", "chrom", "start", "end", "strand",
                                     "rep_name", "rep_class_family", "sw_score",
                                     "pct_div", "ov_start", "ov_end"])
    if not df.empty:
        cf = df.rep_class_family.str.split("/", n=1, expand=True)
        df["rep_class"] = cf[0]
        df["rep_family"] = cf[1].fillna(cf[0])
        df["overlap_bp"] = df.ov_end - df.ov_start
        df = df.drop(columns=["ov_start", "ov_end", "rep_class_family"])
    stats = {"data_lines": n_lines, "repeats_in_window": n_kept,
             "locus_repeat_pairs": n_pairs, "malformed": malformed,
             "unplaced_chroms": len(skipped_chrom),
             "unplaced_lines": sum(skipped_chrom.values())}
    return df, stats


def self_test():
    """Coordinate round-trip and overlap logic on a synthetic case."""
    win = pd.DataFrame({"locus_uid": ["A", "B"], "chrom": ["chrT"] * 2,
                        "start": [5000, 20000], "end": [6000, 21000],
                        "strand": ["+"] * 2})
    win["w0"] = win.start - PAD
    win["w1"] = win.end + PAD
    lines = [
        "  100  1.0 0.0 0.0  chrT   4001   4100 (0) + rA  SINE/Alu   1 100 (0) 1",  # in A window
        "  100  1.0 0.0 0.0  chrT   3999   4000 (0) + rB  SINE/Alu   1 100 (0) 2",  # ends before w0=4000 -> out
        "  100  1.0 0.0 0.0  chrT   7000   7001 (0) C rC  LINE/L1    1 100 (0) 3",  # at w1=7000 boundary
        "  100  1.0 0.0 0.0  chrT   6500   19500 (0) + rD LTR/ERV1   1 100 (0) 4",  # spans both
    ]
    p = "/tmp/_rmsk_selftest.out"
    open(p, "w").write("h1\nh2\n\n" + "\n".join(lines) + "\n")
    df, st = parse(p, win)
    got = {(r.locus_uid, r.rep_name, r.start, r.end, r.strand) for r in df.itertuples()}
    # rA: out 4001-4100 -> 0-based 4000-4100, overlaps A window [4000,7000)
    assert ("A", "rA", 4000, 4100, "+") in got, got
    # rB: out 3999-4000 -> 3999-4000, half-open end == w0 -> no overlap
    assert not any(r[1] == "rB" for r in got), got
    # rC: out 7000-7001 -> 6999-7001, overlaps A (6999 < 7000)
    assert ("A", "rC", 6999, 7001, "-") in got, got
    # rD spans both windows -> emitted twice
    assert sum(r[1] == "rD" for r in got) == 2, got
    assert st["data_lines"] == 4, st
    lens = df.set_index(["locus_uid", "rep_name"]).overlap_bp.to_dict()
    assert lens[("A", "rA")] == 100 and lens[("A", "rC")] == 1, lens
    print("self_test OK:", st)
    return True


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-file")
    ap.add_argument("--assembly", help="locus_coord.assembly value, e.g. t2t or hg38")
    ap.add_argument("--db", default="herv_db/herv_catalog.db")
    ap.add_argument("--parquet")
    ap.add_argument("--pad", type=int, default=PAD)
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args(argv)
    if a.self_test:
        return 0 if self_test() else 1
    if not (a.out_file and a.assembly and a.parquet):
        ap.error("--out-file, --assembly and --parquet are required unless --self-test")
    self_test()          # always run before touching real data
    win = load_windows(a.db, a.assembly, a.pad)
    print(f"windows: {len(win):,} loci on {win.chrom.nunique()} chroms "
          f"| {int((win.w1-win.w0).sum())/1e6:.0f} Mb", flush=True)
    df, st = parse(a.out_file, win)
    df["assembly"] = a.assembly
    df["source"] = os.path.basename(a.out_file)
    print("stats:", st)
    print(f"rows: {len(df):,} | loci with >=1 repeat: {df.locus_uid.nunique():,}")
    print(df.rep_class.value_counts().head(12).to_string())
    df.to_parquet(a.parquet, index=False)
    print("wrote", a.parquet, f"{os.path.getsize(a.parquet)/1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
