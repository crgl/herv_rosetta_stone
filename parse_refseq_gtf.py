#!/usr/bin/env python3
"""
Parse a UCSC ncbiRefSeq GTF into per-transcript exon models, in the SAME record
shape and SAME coordinate convention as the UCSC-API GENCODE models that
build_dashboard.py produces.

COORDINATE CONVENTION -- the whole point of this file
-----------------------------------------------------
GTF is 1-BASED INCLUSIVE.  The UCSC getData/track API returns txStart/exonStarts
0-BASED HALF-OPEN, and the catalog and gene_models.json.gz are 0-based half-open.
So every start here is decremented by 1 and every end is left alone:

    gtf (1-based incl.)  chr1  101  200      -> 100 bp
    bundle (0-based h-o) chr1  100  200      -> 100 bp   len == end - start

Getting this wrong shifts every exon one base left or right, which is invisible
in aggregate and wrong in every graphic.  --self-test asserts it, plus the
degenerate single-base feature and a multi-exon transcript round-trip.

Why RefSeq and not GENCODE for T2T
----------------------------------
hit_gene's t2t rows were computed from hs1.ncbiRefSeq.gtf.gz.  A T2T gene lane
built from anything else would disagree with the overlap table it sits next to --
the same failure mode as the GENCODE V44/V50 mismatch (see build_dashboard.py).
There is no GENCODE annotation for hs1 at UCSC in any case.
"""
from __future__ import annotations
import argparse, gzip, os, sys
from collections import defaultdict

import pandas as pd


def parse_gtf(path: str, want_features=("exon",)) -> pd.DataFrame:
    """
    -> DataFrame[name, name2, chrom, txStart, txEnd, strand, exonStarts, exonEnds]
    name = transcript_id, name2 = gene_name, matching the UCSC API field names
    used by the GENCODE path so downstream code needs no branch.
    """
    tx: dict[str, dict] = {}
    op = gzip.open if path.endswith(".gz") else open
    nbad = 0
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in want_features:
                continue
            try:
                s1, e1 = int(f[3]), int(f[4])
            except ValueError:
                nbad += 1
                continue
            # 1-based inclusive -> 0-based half-open
            s0, e0 = s1 - 1, e1
            attrs = f[8]
            tid = _attr(attrs, "transcript_id")
            if not tid:
                nbad += 1
                continue
            gname = _attr(attrs, "gene_name") or _attr(attrs, "gene_id") or tid
            d = tx.get(tid)
            if d is None:
                d = tx[tid] = {"name": tid, "name2": gname, "chrom": f[0],
                               "strand": f[6], "starts": [], "ends": []}
            d["starts"].append(s0)
            d["ends"].append(e0)
    if nbad:
        print(f"  skipped {nbad:,} malformed/attribute-less lines", file=sys.stderr)

    rows = []
    for d in tx.values():
        order = sorted(range(len(d["starts"])), key=lambda i: d["starts"][i])
        st = [d["starts"][i] for i in order]
        en = [d["ends"][i] for i in order]
        rows.append((d["name"], d["name2"], d["chrom"], st[0], max(en),
                     d["strand"], st, en))
    return pd.DataFrame(rows, columns=["name", "name2", "chrom", "txStart",
                                       "txEnd", "strand", "exonStarts", "exonEnds"])


def _attr(attrs: str, key: str) -> str | None:
    i = attrs.find(key + ' "')
    if i < 0:
        return None
    j = i + len(key) + 2
    k = attrs.find('"', j)
    return attrs[j:k] if k > 0 else None


def self_test() -> bool:
    import tempfile
    gtf = (
        # 1-based inclusive input                       expected 0-based half-open
        'chr1\ts\texon\t101\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GENE1";\n'
        'chr1\ts\texon\t301\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GENE1";\n'
        'chr1\ts\ttranscript\t101\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GENE1";\n'
        'chr2\ts\texon\t50\t50\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_name "GENE2";\n'
        'chr3\ts\texon\t900\t1000\t.\t+\t.\tgene_id "G3"; transcript_id "T3";\n'
    )
    with tempfile.NamedTemporaryFile("w", suffix=".gtf", delete=False) as fh:
        fh.write(gtf)
        p = fh.name
    try:
        G = parse_gtf(p).set_index("name")
        ok = True

        def chk(cond, msg):
            nonlocal ok
            print(("  ok   " if cond else "  FAIL ") + msg)
            ok = ok and cond

        t1 = G.loc["T1"]
        # exon 101-200 (1-based, 100 bp) -> 100-200 (0-based h-o, still 100 bp)
        chk(t1.exonStarts == [100, 300], f"start decremented: {t1.exonStarts} == [100, 300]")
        chk(t1.exonEnds == [200, 400], f"end unchanged: {t1.exonEnds} == [200, 400]")
        chk(all(e - s == 100 for s, e in zip(t1.exonStarts, t1.exonEnds)),
            "exon length preserved (100 bp both conventions)")
        chk((t1.txStart, t1.txEnd) == (100, 400), f"tx extent spans exons: {(t1.txStart, t1.txEnd)}")
        chk(t1.name2 == "GENE1", "gene_name captured")
        chk(len(G) == 3, f"transcript features ignored, only exons grouped: {len(G)} == 3")

        t2 = G.loc["T2"]
        # a single-base GTF feature (50..50) is 1 bp; half-open must be 49..50
        chk((t2.exonStarts, t2.exonEnds) == ([49], [50]),
            f"single-base feature: {(t2.exonStarts, t2.exonEnds)} == ([49], [50])")
        chk(t2.exonEnds[0] - t2.exonStarts[0] == 1, "single-base feature is 1 bp, not 0")

        chk(G.loc["T3"].name2 == "G3", "falls back to gene_id when gene_name absent")
        return ok
    finally:
        os.unlink(p)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf")
    ap.add_argument("--parquet")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()

    if a.self_test:
        return 0 if self_test() else 1
    if not (a.gtf and a.parquet):
        ap.error("--gtf and --parquet are required (or use --self-test)")
    if not self_test():                      # never touch real data on a broken parser
        print("self-test failed, refusing to parse", file=sys.stderr)
        return 1

    G = parse_gtf(a.gtf)
    G.to_parquet(a.parquet, index=False)
    print(f"transcripts {len(G):,} | chroms {G.chrom.nunique()} | "
          f"exons {int(G.exonStarts.map(len).sum()):,} -> {a.parquet} "
          f"({os.path.getsize(a.parquet)/1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
