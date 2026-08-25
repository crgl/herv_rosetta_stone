#!/usr/bin/env python3
"""
Build the static HERV catalog dashboard bundle from herv_catalog.db.

    python build_dashboard.py --db herv_db/herv_catalog.db --out dash
    python build_dashboard.py --db ... --skip-genes        # reuse cached gene models
    python build_dashboard.py --db ... --gencode wgEncodeGencodeCompV44

Produces  <out>/data/{search_index.json.gz, shard_meta.json, loci/<b>.json.gz,
gene_models.json.gz}.  index.html and detail.js are static and are NOT written by
this script -- copy them in (or keep them under version control alongside it).

Stages, each independently skippable:
  1  search index    identifier + classifier keyspaces  ->  search_index.json.gz
  2  locus shards    per-locus detail, djb2(uid) % NB    ->  loci/<b>.json.gz
  3  gene models     GENCODE fetch from UCSC, cached     ->  gene_models.json.gz
  4  validate        resolution + strict-JSON + coverage

--- Three invariants this script exists to protect -------------------------------

DJB2, NOT hash().  Shard buckets are djb2(locus_uid) % NB, implemented identically
here and in index.html.  Python's built-in hash() is salted per process: a bundle
built with it sends every browser lookup to the wrong shard.  Caught in testing;
do not "simplify" back to hash().

allow_nan=False.  json.dump writes bare NaN for float NaN.  Python's json.load
accepts it, JSON.parse in every browser rejects it.  A contaminated shard fetches
with HTTP 200 and then throws inside an async function -- the page silently does
nothing.  This shipped once (332/400 shards, from loci with no cytoband).  Every
dump here goes through _clean() with allow_nan=False so a recurrence is a
build-time exception instead of silent corruption.

GENCODE VERSION.  hit_gene was computed against GENCODE V50, not V44.  Proof:
locus HERVL000012 (chr10:37,930,137-37,934,250) records an intronic ZNF25 overlap
of 4,112 bp == its full span; ZNF25 spans 37,949,572-37,976,647 in V20..V49 (15 kb
away) but 37,916,978-37,976,658 in V50.  Building against V44 leaves 8,080 loci
with recorded overlaps and nothing to draw, which looks like a catalog defect and
is not.  V44 agreement 78.8% / coverage 63%; V50 agreement 98.3% / coverage 92.8%.
Change --gencode only if you also intend the graphic to disagree with hit_gene.
"""
from __future__ import annotations
import argparse, gzip, json, math, os, re, sqlite3, sys, time
from collections import defaultdict

import numpy as np
import pandas as pd

# Subramanian 2011 (HML-2 only) contributes four types. `subramanian_id` uses the
# HML2_ clade prefix so its band names never collide as bare strings with the
# catalog's own combined_id -- but 6 of them still name a DIFFERENT locus than the
# identically-spelled combined_id does (see table subramanian_alias_conflict).
# The three K-series are kept distinct from ERVmap's K-numbers: ERVmap `K-10` and
# `ERVK-10` are different loci, so merging the series would fabricate identities.
IDENT = ["combined_id", "versioned_id", "telescope_id", "hervd_id",
         "ervmap_id", "ervmap_alt_name",
         "subramanian_id", "hml2_k_number", "hml2_k_designation",
         "hgnc_ervk_symbol",
         "missillac_id", "rbrt_id"]
CLASS = ["dfam_int_model", "dfam_accession", "repbase_name"]
DETAIL_TABLES = ["locus_coord", "locus_segment", "locus_structure", "hit_geve",
                 "hit_hervarium_domain", "hit_hervarium_int", "hit_gene",
                 "aln_crossgenome"]
# ---------------------------------------------------------------------------
# Catalog completeness manifest.
#
# herv_catalog.db is stored with `working_data` retention: only the latest copy
# survives, so a rebuild that omits a hand-loaded table cannot be recovered
# from artifacts. That happened once -- a leaner 24-table rebuild became the
# latest, dropping the Subramanian 2011 load (196 locus_alias rows +
# hml2_provirus_detail). The dashboard kept building and shipped four
# identifier groups with zero members for two releases.
#
# Each entry is (table, expected_min_rows, predicate_sql). The build asserts
# these before doing any work and fails loudly on a shortfall. Update the
# expected counts deliberately when a layer legitimately grows.
CATALOG_MANIFEST = [
    ("locus",                39_733, None),
    ("locus_alias",         460_090, None),
    ("locus_alias",              196, "source='subramanian2011'"),
    ("hml2_provirus_detail",      91, None),
    ("subramanian_alias_conflict", 6, None),
    ("locus_repeat",        680_174, None),
    ("locus_coord",           36_746, "assembly='hg38'"),
    ("resource_registry",         1, "resource_key='subramanian2011'"),
    ("missillac_record",     149_057, None),
    ("missillac_locus_map",   33_465, None),
    ("missillac_locus_map",   26_143, "is_primary=1"),
    ("missillac_lineage_group",   83, None),
    ("locus_alias",           52_286, "source='ERV Navigator'"),
    ("resource_registry",          1, "resource_key='ervnav'"),
    # Splice-tier ladder. The row count is the whole catalog because the table
    # carries the negative and non-assessable states explicitly rather than by
    # absence -- a locus missing from it would be indistinguishable from one
    # with no junction evidence, which is the distinction the layer exists for.
    ("locus_splice_tier",      39_733, None),
    ("locus_splice_tier",      35_396, "tier NOT IN ('none','not_assessable')"),
    ("group_lineage_dominant",     93, None),
    ("resource_registry",           1, "resource_key='snaptron_splice_tier'"),
]
# Every manifest entry is fatal by default -- that is the point. A build from a
# deliberately leaner catalog must say so explicitly with --allow-incomplete,
# which downgrades shortfalls to warnings and stamps the bundle metadata so an
# incomplete build is identifiable after the fact.

PAD = 1000          # graphic window padding, must match PAD in detail.js
N_BUCKETS = 400
TILE = 2_000_000    # UCSC query tile size
DEFAULT_TRACK = "wgEncodeGencodeCompV50"


# ---------------------------------------------------------------- primitives

def _has_table(con, name: str) -> bool:
    """True if `name` exists as a table or view in the connected database."""
    return con.execute(
        "SELECT 1 FROM sqlite_master WHERE type IN ('table','view') AND name=?",
        (name,)).fetchone() is not None


# ---------------------------------------------------------------------------
# Evidence-layer input files.
#
# CATALOG_MANIFEST protects tables INSIDE the catalog. It cannot see the
# side-car parquet/db files that supply the arc, mappability, CAGE and TU
# layers -- those are passed as paths, and every one of the loaders is
# deliberately tolerant so a public build without them still works.
#
# That tolerance hid a real regression: --jx-parquet defaulted to a filename
# that never existed, so a default invocation logged one quiet
# "no Snaptron parquet" line and shipped a bundle with zero arcs. The catalog
# manifest passed, the validator passed, and the loss was invisible until a
# human noticed the arcs were gone from the rendered page.
#
# Each entry is (arg_name, attr, layer description). A missing file is now
# reported as a block and is FATAL unless --allow-missing-layers is passed,
# which downgrades to warnings and stamps the bundle metadata.
OPTIONAL_INPUTS = [
    ("--jx-parquet", "jx_parquet", "packed junction arcs (splicing evidence)"),
    ("--iv-parquet", "iv_parquet", "mappability interval blocks + junction anchor bit"),
    ("--f5-parquet", "f5_parquet", "FANTOM5 CAGE lane"),
    ("--tu-db",      "tu_db",      "transcription-unit panel"),
]


def check_inputs(args, allow_missing: bool = False) -> dict:
    """Report presence of every evidence-layer input file.

    Fatal on any miss unless allow_missing is set. This is the file-level twin
    of check_catalog: the tolerant loaders stay tolerant (a deliberate lean
    build is still possible), but silence is no longer the default.
    """
    report, missing = [], []
    for arg, attr, desc in OPTIONAL_INPUTS:
        path = getattr(args, attr, "") or ""
        ok = bool(path) and os.path.exists(path)
        report.append({"arg": arg, "path": path, "present": ok, "layer": desc})
        log(f"{'  ' if ok else '!!'} {arg:<14s} {('ok  ' if ok else 'MISSING')} "
            f"{path or '(unset)':<44s} {desc}")
        if not ok:
            missing.append(f"{arg} -> {path or '(unset)'}: {desc}")

    if missing:
        if not allow_missing:
            log("")
            log("EVIDENCE LAYER INPUT MISSING -- refusing to build:")
            for m in missing:
                log("   " + m)
            log("")
            log("These layers are additive, so the build would otherwise succeed")
            log("and ship a bundle silently missing them -- which is exactly how")
            log("the v0.7 arc regression happened. Supply the files, or pass")
            log("--allow-missing-layers if a lean bundle is deliberate.")
            raise SystemExit("evidence layer input missing; pass --allow-missing-layers")
        log(f"!! building WITHOUT {len(missing)} evidence layer(s) (--allow-missing-layers)")
    return {"inputs": report, "complete": not missing,
            "allow_missing_layers": bool(allow_missing)}


def check_catalog(con, allow_incomplete: bool = False) -> dict:
    """Assert CATALOG_MANIFEST against the connected catalog.

    Returns a manifest report dict (embedded in bundle metadata). Raises
    SystemExit on any shortfall unless allow_incomplete is set.
    """
    report, problems = [], []
    for table, want, pred in CATALOG_MANIFEST:
        label = table + (f" [{pred}]" if pred else "")
        if not _has_table(con, table):
            report.append({"check": label, "expected": want, "got": None,
                           "status": "MISSING TABLE"})
            problems.append(f"{label}: table absent (expected >= {want:,} rows)")
            continue
        sql = f"SELECT COUNT(*) FROM {table}" + (f" WHERE {pred}" if pred else "")
        got = con.execute(sql).fetchone()[0]
        ok = got >= want
        report.append({"check": label, "expected": want, "got": got,
                       "status": "ok" if ok else "SHORT"})
        if not ok:
            problems.append(f"{label}: {got:,} rows, expected >= {want:,}")

    for r in report:
        mark = "  " if r["status"] == "ok" else "!!"
        log(f"{mark} {r['check']:<44s} {str(r['got'] or '-'):>9s} / {r['expected']:>9,}")

    if problems:
        if not allow_incomplete:
            log("")
            log("CATALOG INCOMPLETE -- refusing to build:")
            for p in problems:
                log("   " + p)
            log("")
            log("The catalog is stored with `working_data` retention, so a rebuild")
            log("that drops a table cannot be restored from artifact history. Fix")
            log("the catalog, or pass --allow-incomplete if this is deliberate.")
            raise SystemExit(2)
        log(f"!! {len(problems)} manifest shortfall(s), continuing under "
            f"--allow-incomplete")
    return {"checks": report, "problems": problems,
            "complete": not problems, "allow_incomplete": bool(allow_incomplete)}


def djb2(s: str) -> int:
    """Shard hash. Mirrored byte-for-byte in index.html; see module docstring."""
    h = 5381
    for ch in s:
        h = ((h * 33) + ord(ch)) & 0xFFFFFFFF
    return h


def _clean(o):
    """Coerce non-finite floats to None so the result is legal JSON."""
    if isinstance(o, float):
        return None if (math.isnan(o) or math.isinf(o)) else o
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    return o


def dump_gz(obj, path: str) -> int:
    """Write gzipped JSON. allow_nan=False turns NaN leakage into an exception.

    mtime=0 is load-bearing: gzip stamps the wall clock into its header by
    default, so two builds of identical content produced different checksums
    and could not be compared byte-for-byte. Fixing the stamp makes the
    release/debug index comparison a plain cmp.
    """
    payload = json.dumps(_clean(obj), separators=(",", ":"),
                         allow_nan=False, default=str).encode()
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as fh:
            fh.write(payload)
    return os.path.getsize(path)


def dump_json(obj, path: str) -> int:
    with open(path, "w") as fh:
        json.dump(_clean(obj), fh, separators=(",", ":"), allow_nan=False, default=str)
    return os.path.getsize(path)


MAPPABILITY_PAD = 1000   # bp of flank whose mappability blocks are shipped
PAD_F5 = 1000            # bp of flank whose FANTOM5 CAGE peak marks are shipped


def log(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------- stage 1

def build_search_index(con, out: str) -> pd.DataFrame:
    al = pd.read_sql("SELECT locus_uid,alias,alias_type,assignment,is_current "
                     "FROM locus_alias", con)
    loc = pd.read_sql('SELECT locus_uid,combined_id,versioned_id,"group",band,origin '
                      'FROM locus', con)

    # identifier keyspaces: alias -> [[uid,...], any_current]
    # a list, not a scalar: retired positional strings get re-minted onto other
    # loci, so ~11.5k combined_id strings legitimately resolve to >1 locus. The UI
    # badges these and offers a disambiguation table rather than guessing.
    ident = {}
    for at in IDENT:
        sub = al[al.alias_type == at]
        g = (sub.groupby("alias")
                .agg(u=("locus_uid", lambda s: sorted(set(s))), cur=("is_current", "max"))
                .reset_index())
        ident[at] = {r.alias: [r.u, int(r.cur)] for r in g.itertuples()}
        amb = sum(len(v[0]) > 1 for v in ident[at].values())
        log(f"  {at:16s} {len(g):>7,} strings | ambiguous {amb:>6,}")

    # classifier keyspaces: low cardinality, resolve to a locus LIST
    clsidx = {}
    for at in CLASS:
        sub = al[al.alias_type == at]
        clsidx[at] = {k: sorted(set(v))
                      for k, v in sub.groupby("alias").locus_uid.apply(list).items()}
        med = int(np.median([len(v) for v in clsidx[at].values()])) if clsidx[at] else 0
        log(f"  {at:16s} {len(clsidx[at]):>7,} values  | median loci/value {med}")

    fuzzy = build_fuzzy_index(al, loc)

    # `ident` is NOT shipped: every string in it is present in `fuzzy` tagged
    # with its type, so the UI serves single-keyspace search by filtering fuzzy
    # postings on type. Shipping both would duplicate ~0.9 MB. It is still built
    # above because validate() checks resolution against it.
    #
    # locus_uid likewise needs no keyspace of its own -- fuzzy["uids"] is the
    # sorted uid list, and the UI matches against it directly when the user
    # explicitly selects the internal-key option (last in the dropdown).
    listmeta = build_list_meta(con, loc, fuzzy["uids"])

    misidx = build_missillac_index(con)

    n = dump_gz({"classifier": clsidx, "fuzzy": fuzzy, "listmeta": listmeta,
                 "missillac": misidx,
                 "uid2cid": dict(zip(loc.locus_uid, loc.combined_id)),
                 "uid2group": dict(zip(loc.locus_uid, loc["group"])),
                 "meta": {"n_loci": len(loc), "ident_types": IDENT, "class_types": CLASS}},
                f"{out}/data/search_index.json.gz")
    log(f"  search_index.json.gz  {n/1e6:.2f} MB")
    return loc, al


# Normalisation for the cross-keyspace fuzzy search. Two levels, deliberately:
#
#   nrm()  strips every non-alphanumeric and uppercases.  HERV-K108 -> HERVK108
#   nrmc() additionally drops a leading HERV/ERV.          HERV-K108 -> K108
#
# nrmc is the ONLY way a user typing "HERVK108" reaches the catalog string
# "K108R", but it is NOT a safe silent merge: K-10 and ERVK-10 are different
# loci (HML2_1q22 vs HML2_5q33.3) and ERVK-10 is the same locus as K-11 -- the
# two K-series are independent numbering systems, not offset by a constant.
# 991 nrmc keys merge loci that nrm keeps apart. So nrmc hits are emitted as a
# separate, lower-ranked tier that the UI badges "prefix-collapsed"; they are
# never folded into the nrm keyspace.
_NON_ALNUM = re.compile(r"[^A-Z0-9]")
_ERV_PREFIX = re.compile(r"^(HERV|ERV)")


def nrm(s: str) -> str:
    return _NON_ALNUM.sub("", str(s).upper())


def nrmc(s: str) -> str:
    return _ERV_PREFIX.sub("", nrm(s))


def build_list_meta(con, loc: pd.DataFrame, uids: list) -> dict:
    """Row data for the paged family lists, one entry per locus, uid-index aligned.

    A classifier hit can name 6,815 loci (ERVLE). Rendering that by fetching each
    locus from its shard touches essentially all 400 shards (~22 MB) because uids
    are hash-spread; this table is 0.64 MB for the whole catalog and lets the list
    page, filter and sort entirely client-side with no fetches at all.

    Columns are dictionary-encoded ints where the domain is small. Order matches
    fuzzy["uids"] exactly, so a posting's uid index addresses this table directly.
    """
    lm = loc.set_index("locus_uid").reindex(uids)
    co = (pd.read_sql("SELECT locus_uid,chrom,start,end FROM locus_coord "
                      "WHERE assembly='hg38'", con)
            .drop_duplicates("locus_uid").set_index("locus_uid").reindex(uids))

    groups = sorted(lm["group"].dropna().unique())
    origins = sorted(lm["origin"].fillna("").unique())
    chroms = sorted(co["chrom"].fillna("").unique())
    gi = {v: i for i, v in enumerate(groups)}
    oi = {v: i for i, v in enumerate(origins)}
    ci = {v: i for i, v in enumerate(chroms)}

    rows = []
    for cid, grp, band, org, ch, st, en in zip(
            lm.combined_id, lm["group"], lm.band, lm["origin"].fillna(""),
            co.chrom.fillna(""), co.start, co.end):
        rows.append([cid or "", gi.get(grp, -1), band if pd.notna(band) else "",
                     ci.get(ch, -1),
                     int(st) if pd.notna(st) else -1,
                     int(en) if pd.notna(en) else -1,
                     oi.get(org, -1)])
    log(f"  listmeta         {len(rows):>7,} rows    | client-side paging, no shard fetch")
    return {"groups": groups, "origins": origins, "chroms": chroms, "rows": rows}


def build_fuzzy_index(al: pd.DataFrame, loc: pd.DataFrame) -> dict:
    """Cross-keyspace postings: normalised key -> [[uid_idx, type_idx, original, is_current]].

    Covers every IDENT keyspace EXCEPT locus_uid, which is not an alias type and
    is exposed only via its own explicitly-selected dropdown option.

    ervmap_alt_name packs several names into one string ('K108R, ERVK-6'), so it
    is split on commas before indexing -- without the split, searching K108R
    misses it entirely.
    """
    uids = sorted(loc.locus_uid)
    uid_ix = {u: i for i, u in enumerate(uids)}
    type_ix = {t: i for i, t in enumerate(IDENT)}

    post = defaultdict(list)
    seen = set()
    sub = al[al.alias_type.isin(IDENT)]
    for r in sub.itertuples():
        # index the packed string AND its parts, so this index is a strict
        # superset of the per-keyspace `ident` index (which stores only the
        # packed form) and can replace it outright.
        toks = [str(r.alias)]
        if r.alias_type == "ervmap_alt_name":
            toks += [t.strip() for t in str(r.alias).split(",")]
        for tok in toks:
            if not tok:
                continue
            sig = (tok, r.alias_type, r.locus_uid)
            if sig in seen:
                continue
            seen.add(sig)
            post[nrm(tok)].append([uid_ix[r.locus_uid], type_ix[r.alias_type],
                                   tok, int(r.is_current)])

    keys = sorted(post)
    key_ix = {k: i for i, k in enumerate(keys)}
    # Collapsed keyspace -> the nrm keys it covers, so the collapsed tier reuses
    # the tier-1/2 postings instead of duplicating them.
    #
    # EVERY key is entered under its collapsed form, not just the ones that
    # actually carry a HERV/ERV prefix. The user's case is exactly the reverse of
    # the obvious one: they type "HERVK108" (prefixed) and the catalog string is
    # "K108R" (unprefixed). If only prefixed keys were collapsed, "K108R" would
    # never appear under collapsed key "K108" and the query would miss.
    coll = defaultdict(list)
    for k in keys:
        coll[_ERV_PREFIX.sub("", k)].append(key_ix[k])
    ckeys = sorted(coll)

    n_amb = sum(len({p[0] for p in v}) > 1 for v in post.values())
    log(f"  fuzzy            {len(keys):>7,} keys    | postings {sum(len(v) for v in post.values()):,}"
        f" | multi-locus keys {n_amb:,}")
    log(f"  fuzzy-collapsed  {len(ckeys):>7,} keys    | (lower-ranked tier, badged in UI)")
    return {"uids": uids, "types": IDENT, "keys": keys,
            "post": [post[k] for k in keys],
            "ckeys": ckeys, "cmap": [coll[c] for c in ckeys]}


def build_missillac_index(con) -> dict:
    """Search lane for ERV Navigator records that resolve to NO catalog locus.

    The ~120k solo LTRs with no overlapping locus are deliberately not loci (they
    would quadruple the catalog and re-letter the positional keyspace), so they
    carry no locus_alias row and cannot appear in the uid-indexed fuzzy postings.
    Without this lane a user pasting a perfectly valid Missillac ID gets nothing
    back, which reads as a broken index rather than an out-of-scope record.

    Postings here resolve to a record row, not a locus: the UI shows category,
    coordinates and an outbound ERV Navigator link, badged as 'no catalog locus'.
    Mapped records are excluded -- they already resolve through the locus path.
    """
    if not _has_table(con, "missillac_record"):
        log("  missillac index       ABSENT (declared optional in manifest)")
        return {}
    rec = pd.read_sql("""
        SELECT r.missillac_id, r.rbrt_id, r.category, r.lineage_id,
               r.chrom_hg38, r.start_hg38, r.end_hg38, r.url
        FROM missillac_record r
        WHERE r.missillac_id NOT IN (SELECT missillac_id FROM missillac_locus_map)
        ORDER BY r.missillac_id""", con)
    cats = sorted(rec.category.dropna().unique())
    cat_ix = {c: i for i, c in enumerate(cats)}
    rows = [[r.missillac_id, r.rbrt_id, cat_ix.get(r.category, -1), r.lineage_id,
             r.chrom_hg38, (None if pd.isna(r.start_hg38) else int(r.start_hg38)),
             (None if pd.isna(r.end_hg38) else int(r.end_hg38)), r.url]
            for r in rec.itertuples()]
    post = defaultdict(list)
    for i, r in enumerate(rec.itertuples()):
        for tok in (r.missillac_id, r.rbrt_id):
            if tok:
                post[nrm(str(tok))].append(i)
    keys = sorted(post)
    log(f"  missillac unmapped {len(rows):>7,} records | {len(keys):,} fuzzy keys")
    return {"cats": cats, "rows": rows, "keys": keys,
            "post": [post[k] for k in keys]}


# ---------------------------------------------------------------- stage 2

def _btx():
    """Import the sibling transcription module regardless of cwd."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import build_transcription as btx
    return btx


def _tx_payload(tu_db: str, jx_parquet: str, want_arcs: bool = True):
    """(locus tx records, locus arcs) -- both {} when the inputs are absent.

    Kept tolerant on purpose: the catalog database alone is enough to build a
    working bundle, and the transcription layer is additive.

    want_arcs=False skips the debug arc layer outright. Only --debug-local puts
    `arcs` in a shard, and nothing else reads it (the validator does not), so in a
    public build computing it is pure cost -- it dominated the shard stage."""
    if not tu_db or not os.path.exists(tu_db):
        log("  transcription        no TU database -- skipping tx/arc panels")
        return {}, {}
    btx = _btx()
    cx = sqlite3.connect(tu_db)
    tx = btx.locus_tx(cx)
    log(f"  locus tx {len(tx):>6,} loci   | FANTOM5 TSS + splice status")
    arcs = {}
    if not want_arcs:
        log("  locus arcs           skipped (public build -- packed `pjx` layer only)")
    elif jx_parquet and os.path.exists(jx_parquet):
        arcs = btx.locus_arcs(cx, jx_parquet)
        capped = sum(1 for v in arcs.values() if v["n_total"] > v["shown"])
        log(f"  locus arcs {len(arcs):>6,} loci | Snaptron hg38 | {capped:,} capped "
            f"at top {btx.CAP} by support")
    else:
        log("  locus arcs           no Snaptron parquet -- tabulation only, no arcs")
    return tx, arcs


def _splice_tier_payload(con) -> dict:
    """{locus_uid: {"t": tier, "d": detail, "sc": max_sample_count}}.

    Asserted present rather than tolerated: unlike the drawn evidence lanes,
    this is a per-locus summary line that a reader will take as a statement
    about the locus. A silently absent layer here would render as "no splicing
    evidence" on 35,396 loci that have it -- indistinguishable from a real
    negative, which is the failure mode the tri-state tier exists to prevent.
    """
    if not _has_table(con, "locus_splice_tier"):
        raise SystemExit("locus_splice_tier missing -- build the splice layer "
                         "or remove it from CATALOG_MANIFEST deliberately")
    df = pd.read_sql("SELECT locus_uid, tier, detail, max_sc FROM locus_splice_tier", con)
    log(f"  splice tier          {len(df):,} loci | "
        + " ".join(f"{k}={v:,}" for k, v in df.tier.value_counts().items()))
    return {r.locus_uid: {"t": r.tier, "d": r.detail, "sc": r.max_sc}
            for r in df.itertuples()}


def _fantom5_payload(con, pk_parquet: str = ""):
    """Per-locus FANTOM5 CAGE evidence: window counts from locus_fantom5, plus
    the individual dominant-TSS positions so the lane draws real peak marks
    rather than a summary bar.

    This resource shares primary data with FANTOM CAT (CAT clusters are built
    from FANTOM5 CAGE), so the renderer labels it as non-independent: agreement
    between the two CAGE lanes is one observation, not two.

    Single tier by design -- the distributed peak set is already thresholded
    (only 171 of 209,911 peaks fall below 1 TPM in every library), so sample
    breadth is carried as a continuous value instead of a tier name. A locus
    with assessable=0 (alt/random contig carrying no peaks) yields an explicit
    null so "not measured here" stays distinct from "measured as zero".
    """
    try:
        lf = pd.read_sql("SELECT * FROM locus_fantom5", con)
    except Exception as exc:
        log(f"  fantom5: absent ({exc.__class__.__name__}) -- lane omitted")
        return {}, {}
    stats = defaultdict(dict)
    for r in lf.itertuples():
        if int(r.assessable) == 0:
            stats[r.locus_uid][str(r.asm)] = {"no_data": 1}
            continue
        stats[r.locus_uid][str(r.asm)] = {
            "no_data": 0,
            "body": int(r.body_sense_n), "fkb": int(r.fkb_sense_n),
            "up": int(r.up_sense_n), "anti": int(r.body_anti_n),
            "max_tpm": (None if pd.isna(r.body_sense_max_tpm)
                        else round(float(r.body_sense_max_tpm), 2)),
            "breadth": (None if pd.isna(r.body_sense_breadth_n_ge1)
                        else int(r.body_sense_breadth_n_ge1)),
            "d5p": (None if pd.isna(r.tss5p_dist) else int(r.tss5p_dist)),
            "top": (None if (r.body_sense_top_sample is None
                             or pd.isna(r.body_sense_top_sample))
                    else str(r.body_sense_top_sample)[:60]),
            "new": int(r.body_sense_n_new_peaks),
        }
    marks = defaultdict(dict)
    if pk_parquet and os.path.exists(pk_parquet):
        pk = pd.read_parquet(pk_parquet)
        idx = {}
        for (asm, c), g in pk.groupby(["assembly", "chrom"], observed=True):
            g = g.sort_values("tss")
            idx[(str(asm), str(c))] = g
        co = pd.read_sql(
            "SELECT locus_uid,assembly,chrom,start,end,strand FROM locus_coord", con)
        for r in co.itertuples():
            if pd.isna(r.chrom) or pd.isna(r.start):
                continue
            g = idx.get((str(r.assembly), str(r.chrom)))
            if g is None:
                continue
            s = max(0, int(r.start) - PAD_F5)
            e = int(r.end) + PAD_F5
            tp = g.tss.to_numpy()
            i = int(np.searchsorted(tp, s, "left"))
            j = int(np.searchsorted(tp, e, "left"))
            if j <= i:
                continue
            sub = g.iloc[i:j]
            # [pos, sense(1/0), log10-ish tpm, breadth, hg38-native flag]
            marks[r.locus_uid][str(r.assembly)] = [
                [int(t), int(sd == str(r.strand)), round(float(mt), 1), int(nb), int(sc == "hg38")]
                for t, sd, mt, nb, sc in zip(sub.tss, sub.strand, sub.max_tpm,
                                             sub.n_ge1, sub.src)]
    return dict(stats), dict(marks)


def _packed_jx_payload(con, jx_parquet: str, iv_parquet: str = "") -> dict:
    """{locus_uid: {"w": [uint64, ...], "n": pre-cap junction count}} or {}.

    Kept tolerant like _tx_payload: the packed layer is additive, and a build
    without the Snaptron parquet still produces a working bundle.
    """
    if not jx_parquet or not os.path.exists(jx_parquet):
        log("  packed jx            no Snaptron parquet -- no packed arcs")
        return {}
    import pack_junctions as pjx
    loci = pd.read_sql("SELECT locus_uid, chrom, start, end, strand FROM locus_coord "
                       "WHERE assembly='hg38'", con).dropna(subset=["chrom", "start", "end"])
    loci["start"] = loci.start.astype(int)
    loci["end"] = loci.end.astype(int)
    # pm151 supplies the anchor-mappable bit. Panmask is hg38-only, which is fine
    # here: the srav3h junction set is hg38-only too. Without it the bit would be
    # 0 everywhere, indistinguishable from "measured unmappable", so the build
    # asserts below that it was actually measured.
    anchor_iv = None
    if iv_parquet and os.path.exists(iv_parquet):
        _iv = pd.read_parquet(iv_parquet)
        anchor_iv = _iv[(_iv.resource == "pm151") & (_iv.assembly == "hg38")][
            ["chrom", "start", "end"]].copy()
        anchor_iv["chrom"] = anchor_iv.chrom.astype(str)
    words, st, keep = pjx.build(jx_parquet, loci, canonical_only=True,
                                anchor_intervals=anchor_iv)
    pre = keep.groupby("locus_uid").size().to_dict()
    # Words go out as DECIMAL STRINGS, not JSON numbers: a 64-bit value exceeds
    # Number.MAX_SAFE_INTEGER (2^53), so JSON.parse would round the top bits away
    # and corrupt the donor offset. The JS side reads them with BigInt.
    out = {u: {"w": [str(x) for x in w], "n": int(pre.get(u, len(w)))}
           for u, w in words.items()}
    log(f"  packed jx {st['loci']:>6,} loci | {st['junctions']:,} junctions | "
        f"{st['antisense']:,} antisense | {st['saturated']} saturated | "
        f"max {st['max_per_locus'] * 8} B/locus")
    log(f"  anchor bit           {st['anchor_ok']:,} set ({100 * st['anchor_ok'] / max(st['junctions'], 1):.1f}%)"
        f" | {st['anchor_nodata']:,} without pm151 coverage")
    if anchor_iv is not None and st["anchor_nodata"]:
        # A clear bit must mean "neither anchor easy", never "no data". Measured
        # over all 1,425,453 distinct junctions this count is zero, so a nonzero
        # value means the interval reference no longer covers the junction set.
        raise AssertionError(
            f"FAIL: {st['anchor_nodata']:,} junctions have an anchor window on a "
            f"contig pm151 does not cover -- a clear anchor bit would be ambiguous")
    if st["saturated"]:
        # 23/24-bit signed fields carry the observed offset range with headroom;
        # a nonzero count means an assumption in pack_junctions has been violated.
        log(f"  WARNING: {st['saturated']} junction offsets saturated a packed field")
    return out


def build_tu_shards(cat_con, out: str, tu_db: str, jx_parquet: str, nb: int):
    """Stage 5, --debug-local only: per-unit shards for the TU detail view.

    Deliberately gated. This roughly doubles bundle size, and a static host with a
    size ceiling should carry the locus-centric public view only. The TU view is
    the internal tool for adjudicating splits, so it lives in the debug bundle
    where it can grow (UMAP neighbourhoods, per-unit read evidence) without
    threatening what ships publicly.

    Shards use the same djb2(key) % nb convention as the locus shards -- keyed on
    tu_id -- so index.html resolves both with one code path."""
    btx = _btx()
    cx = sqlite3.connect(tu_db)
    tx = btx.tu_tx(cx)
    snap = btx.tu_snaptron(cx)
    arcs = btx.tu_arcs(cx, cat_con, jx_parquet) if os.path.exists(jx_parquet) else {}
    # members come from their own hg38/t2t-fallback pass, not from the arc payload:
    # Snaptron is hg38-only but 15% of units exist only in t2t
    mem = btx.tu_members(cx, cat_con)
    ids = sorted(set(tx) | set(snap) | set(arcs) | set(mem))
    log(f"  units {len(ids):>6,} | tx {len(tx):,} | snaptron {len(snap):,} | arcs {len(arcs):,}")
    n_t2t = sum(1 for v in mem.values() if v["asm"] == "t2t")
    n_nom = sum(1 for v in mem.values() if v.get("no_members"))
    log(f"  member tracks {len(mem):>6,} | hg38 {len(mem)-n_t2t:,} + t2t fallback {n_t2t:,}")
    n_rc = sum(1 for v in arcs.values() if v.get("recount"))
    if n_rc:
        log(f"  {n_rc:,} units: recomputed arc classes differ from the upstream "
            f"summary counts -- both shown on the page")
    if n_nom:
        log(f"  {n_nom:,} units have NO tu_locus_member row (upstream gap) -- "
            f"drawn as unmapped extent")
    unanchored = sum(1 for v in arcs.values() if v.get("unanchored"))
    if unanchored:
        log(f"  {unanchored:,} chimeric units have window junctions but none "
            f"element-to-element -- shown as unanchored")

    shards = {}
    for t in ids:
        a = arcs.get(t, {})
        m = mem.get(t, {})
        shards.setdefault(djb2(t) % nb, {})[t] = {
            "tu_id": t, "tx": tx.get(t, {}), "snaptron": snap.get(t, {}),
            "arcs": a, "members": m.get("members", []),
            "map_asm": m.get("asm"), "map_chrom": m.get("chrom"),
            "extent": m.get("extent"), "no_members": bool(m.get("no_members")),
        }
    d = f"{out}/data/tu"
    os.makedirs(d, exist_ok=True)
    for stale in os.listdir(d):
        os.remove(f"{d}/{stale}")
    sizes = [dump_gz(v, f"{d}/{b}.json.gz") for b, v in shards.items()]
    dump_json({"n_buckets": nb, "hash": "djb2_mod", "n_units": len(ids)},
              f"{out}/data/tu_meta.json")
    # the TU search keyspace: tu_id -> member locus uids, for the exact-match box
    dump_gz({t: [x[0] for x in mem.get(t, {}).get("members", [])] for t in ids},
            f"{out}/data/tu_index.json.gz")
    log(f"  tu shards {len(shards)} | {sum(sizes)/1e6:.1f} MB | "
        f"max {max(sizes)/1024:.0f} KB")


def _mappability_payload(con, loc: pd.DataFrame, iv_parquet: str = ""):
    """Per-locus mappability: summary stats from the catalogue, plus the raw
    intervals clipped to the drawn window so the lane can be drawn as blocks.

    Two resources on hg38 (pm151 Panmask "easy", Umap k100 unique) and one on
    T2T (Umap k100 -- Panmask has no T2T release).  Intervals are clipped to
    locus +/- MAPPABILITY_PAD rather than shipped whole: an easy region can run
    to 75 kb, and only the drawn window is ever needed.

    A resource with `no_data` for a locus (alt/scaffold contigs, which none of
    the three interval sets cover) yields an explicit null-fraction row so the
    renderer can distinguish "not measured here" from "measured as unmappable".
    """
    try:
        lm = pd.read_sql("SELECT * FROM locus_mappability", con)
    except Exception as exc:
        log(f"  mappability: absent ({exc.__class__.__name__}) -- lane omitted")
        return {}, {}
    stats = defaultdict(dict)
    for r in lm.itertuples():
        stats[r.locus_uid][f"{r.resource}_{r.assembly}"] = {
            "frac": (None if pd.isna(r.mappable_frac) else round(float(r.mappable_frac), 4)),
            "cov": (None if pd.isna(r.covered_bp) else int(r.covered_bp)),
            "blocks": (None if pd.isna(r.n_blocks) else int(r.n_blocks)),
            "longest_unmap": (None if pd.isna(r.longest_unmappable) else int(r.longest_unmappable)),
            "t5": (None if pd.isna(r.term5_mappable) else int(r.term5_mappable)),
            "t3": (None if pd.isna(r.term3_mappable) else int(r.term3_mappable)),
            "no_data": int(r.no_data),
        }
    blocks = defaultdict(dict)
    if iv_parquet and os.path.exists(iv_parquet):
        iv = pd.read_parquet(iv_parquet)
        idx = {}
        for (rs, asm, c), g in iv.groupby(["resource", "assembly", "chrom"], observed=True):
            g = g.sort_values("start")
            idx[(str(rs), str(asm), str(c))] = (g.start.to_numpy(), g.end.to_numpy())
        # `loc` is the search-index frame and carries no coordinates; the padded
        # windows come from locus_coord, which is keyed by (locus_uid, assembly)
        # and is what every other coordinate lane in this builder reads.
        co = pd.read_sql("SELECT locus_uid,assembly,chrom,start,end FROM locus_coord", con)
        for r in co.itertuples():
            if pd.isna(r.chrom) or pd.isna(r.start):
                continue
            s = max(0, int(r.start) - MAPPABILITY_PAD)
            e = int(r.end) + MAPPABILITY_PAD
            asm = str(r.assembly)
            for rs in ("pm151", "umap100"):
                t = idx.get((rs, asm, str(r.chrom)))
                if t is None:
                    continue
                st, en = t
                i = int(np.searchsorted(en, s, "right"))
                j = int(np.searchsorted(st, e, "left"))
                if j <= i:
                    blocks[r.locus_uid][f"{rs}_{asm}"] = []
                    continue
                blocks[r.locus_uid][f"{rs}_{asm}"] = [
                    [int(max(a, s)), int(min(b, e))] for a, b in zip(st[i:j], en[i:j])]
    return dict(stats), dict(blocks)


def build_shards(con, out: str, loc: pd.DataFrame, al: pd.DataFrame, nb: int,
                 tu_db: str = "", jx_parquet: str = "", debug: bool = False,
                 iv_parquet: str = "", f5_parquet: str = ""):
    grp = pd.read_sql('SELECT * FROM "group"', con).set_index("group")
    # Group-level dominant Navigator lineage. Carried WITH its share and the
    # group's lineage count, never bare: lineage and group are many-to-many
    # (median 8 lineages per group, max 51), and only 21 of 93 groups have a
    # lineage covering >=80% of their loci. A bare value would read as an
    # equivalence the data does not support.
    if _has_table(con, "group_lineage_dominant"):
        _gld = pd.read_sql("SELECT grp, dom_lineage, dom_frac, n_lineages "
                           "FROM group_lineage_dominant", con).set_index("grp")
        grp = grp.join(_gld, how="left")
    else:
        log("  group lineage        ABSENT (declared optional in manifest)")
    sf = pd.read_sql("SELECT * FROM superfamily", con).set_index("superfamily")
    tabs = {t: pd.read_sql(f"SELECT * FROM {t}", con) for t in DETAIL_TABLES}
    dbest = pd.read_sql("SELECT locus_uid,dfam_accession,consensus_name,pct_identity,"
                        "cons_cov,sw_score,aln_quality,is_best_by_sw_score "
                        "FROM hit_dfam_aln WHERE is_best_by_sw_score=1", con)

    byuid = {k: {u: d for u, d in v.groupby("locus_uid")} for k, v in tabs.items()}
    dby = {u: d for u, d in dbest.groupby("locus_uid")}
    # Subramanian 2011 panel. This table was silently lost once when a leaner
    # catalog rebuild became the `working_data` latest (see CATALOG_MANIFEST):
    # a bare try/except reported "absent" and the build carried on shipping an
    # empty identifier group. Presence is now a manifest assertion; only a
    # genuine query error is tolerated here, and it is re-raised.
    if _has_table(con, "hml2_provirus_detail"):
        hd = pd.read_sql("SELECT * FROM hml2_provirus_detail "
                         "WHERE locus_uid IS NOT NULL", con)
        hdby = {u: d for u, d in hd.groupby("locus_uid")}
        log(f"  hml2_provirus_detail {len(hd):>4,} rows | {len(hdby)} loci")
    else:
        hdby = {}
        log("  hml2_provirus_detail  ABSENT (declared optional in manifest)")
    # ERV Navigator (Missillac) panel. Records are a standalone resource layer,
    # not loci: the ~120k solo LTRs with no catalog counterpart live in
    # missillac_record with no mapping row, so they are searchable and linkable
    # without touching the locus keyspace. Only mapped records reach a panel.
    # Every overlapping pair is kept; is_primary marks the best Jaccard, and
    # primary_margin exposes a thin call rather than hiding it.
    if _has_table(con, "missillac_locus_map") and _has_table(con, "missillac_record"):
        mn = pd.read_sql("""
            SELECT m.locus_uid, m.missillac_id, m.rbrt_id, m.category,
                   m.ovl_bp, m.jaccard, m.is_primary, m.primary_margin,
                   m.n_loci_for_record, m.strand_agree,
                   r.lineage_id, r.clade, r.url,
                   r.chrom_hg19, r.start_hg19, r.end_hg19, r.lift_status,
                   g.majority_group, g.purity, g.name_match_group,
                   g.name_match_level, g.name_agrees_with_coords, g.confidence
            FROM missillac_locus_map m
            JOIN missillac_record r USING (missillac_id)
            LEFT JOIN missillac_lineage_group g ON g.lineage_id = r.lineage_id
            ORDER BY m.locus_uid, m.is_primary DESC, m.jaccard DESC""", con)
        mnby = {u: d for u, d in mn.groupby("locus_uid")}
        log(f"  missillac panel      {len(mn):>4,} rows | {len(mnby)} loci")
    else:
        mnby = {}
        log("  missillac_locus_map  ABSENT (declared optional in manifest)")
    alby = {u: d for u, d in al.groupby("locus_uid")}
    rby = repeats_for(con, loc.locus_uid)
    # transcription evidence (optional: needs the TU database and the Snaptron
    # parquet; a public build without them simply omits the panels)
    _txl, _arcl = _tx_payload(tu_db, jx_parquet, want_arcs=debug)
    # Bit-packed junction reference for the NON-debug bundle.
    #
    # `arcs` above is the debug payload: full coordinates, strand and canonical
    # flag per junction, capped at CAP by support. At ~36 MB across the shard set
    # it is too large to ship, and it is built ONLY under --debug-local. `pjx`
    # carries the same evidence as one 64-bit word per junction (see
    # pack_junctions for the layout) at ~128 B/locus worst case, and is what the
    # public bundle draws from. Exactly one of the two is populated per build.
    _pjx = _packed_jx_payload(con, jx_parquet, iv_parquet)
    # Mappability lane (pm151 Panmask + Umap k100). Summary stats always ship;
    # the clipped interval blocks need the standalone parquet reference.
    _mapst, _mapbl = _mappability_payload(con, loc, iv_parquet)
    # FANTOM5 CAGE lane REMOVED at v0.8 (user decision): the hit rate was too
    # low to justify the vertical space it occupied in the graphic. The catalog
    # tables (locus_fantom5) and the parquet are untouched, and _fantom5_payload
    # is retained below, so restoring the lane is a renderer change plus one
    # call here -- no recomputation. The --f5-parquet flag is likewise retained.
    # Per-locus splice-evidence tier (Snaptron srav3h). One short string per
    # locus, computed at build time from the FULL filtered pair set rather than
    # the display-capped arc set -- see locus_splice_tier's registry row.
    _splice = _splice_tier_payload(con)
    # v0.1 transcriptional-unit layer for the LOCUS page's unit panel.
    #
    # These tables live in the TU database, not the catalogue, and the crosswalk is
    # named locus_v1_to_tu_v01 (not locus_tu). The original read looked for both in
    # `con` under the wrong name, so the except branch fired on every build and the
    # locus-page unit panel had never once rendered -- a silent skip, logged as if
    # the layer were legitimately absent. Read from tu_db, and log the miss loudly
    # enough that a genuinely absent layer is distinguishable from a wrong lookup.
    _tu, _ltby = None, {}
    if tu_db and os.path.exists(tu_db):
        tcon = sqlite3.connect(tu_db)
        try:
            _tu = pd.read_sql("SELECT * FROM tu", tcon).set_index("tu_id")
            _lt = pd.read_sql("SELECT * FROM locus_v1_to_tu_v01", tcon)
            _ltby = {u: d for u, d in _lt.groupby("locus_uid")}
            log(f"  tu layer {len(_tu):>6,} units | {len(_ltby):,} loci crosswalked")
        except Exception as e:
            log(f"  tu layer  UNAVAILABLE in {tu_db}: {type(e).__name__}: {e}")
        finally:
            tcon.close()
    else:
        log("  tu layer             no --tu-db given -- omitting unit panel")

    def tu_for(u):
        if _tu is None:
            return {}
        d = _ltby.get(u)
        if d is None or not len(d):
            return {}
        r = d.iloc[0]
        out = {"in_v01": bool(r.get("in_v01")),
               "n_tu_overlap_hg38": (None if pd.isna(r.get("n_tu_overlap")) else int(r["n_tu_overlap"])),
               "dom_frac_hg38": (None if pd.isna(r.get("dom_frac")) else round(float(r["dom_frac"]), 3)),
               "n_tu_overlap_t2t": (None if pd.isna(r.get("n_tu_t2t")) else int(r["n_tu_t2t"])),
               "dom_frac_t2t": (None if pd.isna(r.get("dom_frac_t2t")) else round(float(r["dom_frac_t2t"]), 3)),
               "units": []}
        seen = set()
        for col, asm in (("tu_id_hg38", "hg38"), ("tu_id_t2t", "t2t")):
            tid = r.get(col)
            if tid is None or pd.isna(tid) or tid in seen:
                continue
            seen.add(tid)
            if tid in _tu.index:
                rec = _tu.loc[tid].to_dict()
                rec["tu_id"] = tid
                rec["matched_via"] = asm
                out["units"].append(json.loads(pd.Series(rec).to_json()))
        return out

    GKEYS = ("superfamily", "herv_class", "n_loci", "intModel", "repbase_class",
             "hervd_family", "dfam_accession", "dominant_ltr",
             "frac_with_flanking_ltr", "n_with_hervarium_domain", "extension_verdict",
             "dom_lineage", "dom_frac", "n_lineages")

    def recs(d, cols=None, drop=("locus_uid", "resource_key")):
        if d is None:
            return []
        d = d.drop(columns=[c for c in drop if c in d.columns])
        if cols:
            d = d[[c for c in cols if c in d.columns]]
        return json.loads(d.to_json(orient="records"))

    shards = defaultdict(dict)
    for r in loc.itertuples():
        u = r.locus_uid
        g = grp.loc[r.group].to_dict() if r.group in grp.index else {}
        gs = {k: g.get(k) for k in GKEYS}
        s = gs.get("superfamily")
        shards[djb2(u) % nb][u] = {
            "uid": u, "combined_id": r.combined_id, "versioned_id": r.versioned_id,
            "group": r.group, "band": r.band, "origin": r.origin,
            "aliases": recs(alby.get(u), ["alias", "alias_type", "assignment", "is_current"]),
            "group_info": gs,
            "superfamily_info": (sf.loc[s].to_dict() if s in sf.index else {}),
            "coord": recs(byuid["locus_coord"].get(u)),
            "structure": (recs(byuid["locus_structure"].get(u)) or [{}])[0],
            "segments": recs(byuid["locus_segment"].get(u)),
            "geve": recs(byuid["hit_geve"].get(u),
                         drop=("locus_uid", "resource_key", "telescope_id")),
            "domains": recs(byuid["hit_hervarium_domain"].get(u),
                            drop=("locus_uid", "resource_key", "telescope_id")),
            "hervarium_int": recs(byuid["hit_hervarium_int"].get(u)),
            "genes": recs(byuid["hit_gene"].get(u)),
            "tu": tu_for(u),
            # FANTOM5 TSS/splice status per assembly, and Snaptron arcs (hg38 only,
            # capped -- see build_transcription for the tri-state and cap semantics)
            "tx": _txl.get(u, {}),
            # Two junction payloads, and only ONE ships in a given build.
            #
            #   pjx   -- bit-packed, 8 B/junction, always present. What the public
            #            bundle draws from.
            #   arcs  -- full uncapped-coordinate records with exact sample counts,
            #            --debug-local ONLY. It is ~36 MB across the shard set, which
            #            is the entire cost the packing exists to avoid; shipping both
            #            would make the packed layer pure overhead.
            #
            # The renderer prefers `arcs` when present, so a debug build still shows
            # exact values while the public build shows the log-quantised ones.
            # `_arcl` is empty unless debug, so this is {} in a public build.
            "arcs": _arcl.get(u, {}),
            "pjx": _pjx.get(u, {}),
            "crossgenome": (recs(byuid["aln_crossgenome"].get(u)) or [{}])[0],
            # Subramanian 2011 HML-2 detail: age / ORFs / polymorphism. Present for
            # 87 loci only -- the detail view omits the panel when absent.
            "hml2_detail": recs(hdby.get(u)),
            # ERV Navigator: intactness category, Vargiu-style lineage path and a
            # direct locus link per Missillac record overlapping this locus.
            # Multiple rows are expected and meaningful -- our internal-only spans
            # can sit inside one longer Navigator extent.
            "missillac": recs(mnby.get(u)),
            "dfam_best": recs(dby.get(u)),
            "repeats": rby.get(u, []),
            "mappability": {"stats": _mapst.get(u, {}), "blocks": _mapbl.get(u, {})},
            # "fantom5" removed at v0.8 -- lane dropped for low hit rate.
            "splice": _splice.get(u, {}),
        }

    os.makedirs(f"{out}/data/loci", exist_ok=True)
    for stale in os.listdir(f"{out}/data/loci"):
        os.remove(f"{out}/data/loci/{stale}")
    sizes = [dump_gz(v, f"{out}/data/loci/{b}.json.gz") for b, v in shards.items()]
    dump_json({"n_buckets": nb, "hash": "djb2_mod"}, f"{out}/data/shard_meta.json")
    log(f"  shards {len(shards)} | {sum(sizes)/1e6:.1f} MB | mean {np.mean(sizes)/1024:.0f} KB"
        f" | max {max(sizes)/1024:.0f} KB")
    log(f"  loci written {sum(len(v) for v in shards.values()):,} of {len(loc):,}")


# ---------------------------------------------------------------- stage 3

def gene_windows(con, assembly: str = "hg38") -> pd.DataFrame:
    """
    Loci needing gene models on `assembly`: a coordinate on that assembly AND
    >=1 hit_gene row for ANY genome, padded by PAD.

    Deliberately NOT filtered to hit_gene.genome == assembly.  Doing that drops
    191 hg38 loci, and 32 of those really do have hg38 transcripts in window --
    unnamed GENCODE genes (bare ENSG accessions) that hit_gene has no hg38 row
    for.  Filtering by genome would silently delete their gene lane.  The looser
    rule costs only empty entries and keeps hg38 at the verified 22,524 loci.
    """
    gu = pd.read_sql("SELECT DISTINCT locus_uid FROM hit_gene", con).locus_uid
    co = pd.read_sql("SELECT locus_uid,chrom,start,end FROM locus_coord "
                     "WHERE assembly=?", con, params=[assembly])
    need = co[co.locus_uid.isin(set(gu))].copy()
    need["w0"] = (need.start - PAD).clip(lower=0)
    need["w1"] = need.end + PAD
    return need


def fetch_gencode(need: pd.DataFrame, track: str, cache: str) -> pd.DataFrame:
    """Tile-fetch a GENCODE track from the UCSC API. Cached per tile; resumable."""
    import requests
    os.makedirs(cache, exist_ok=True)
    tiles = sorted({(c, t) for c, a, b in zip(need.chrom, need.w0, need.w1)
                    for t in range(a // TILE, (b // TILE) + 1)})
    keep = ("name", "name2", "chrom", "txStart", "txEnd", "strand",
            "exonStarts", "exonEnds")
    rows, fails = [], []
    for i, (c, t) in enumerate(tiles):
        fp = f"{cache}/{c}_{t}.json"
        if os.path.exists(fp):
            rows += json.load(open(fp))
            continue
        got = None
        for k in range(3):
            try:
                r = requests.get("https://api.genome.ucsc.edu/getData/track", timeout=120,
                                 params={"genome": "hg38", "track": track, "chrom": c,
                                         "start": t * TILE, "end": (t + 1) * TILE})
                r.raise_for_status()
                j = r.json()
                v = j.get(track, j)
                v = v if isinstance(v, list) else v.get(c, [])
                got = [{kk: g.get(kk) for kk in keep} for g in v]
                break
            except Exception:
                if k == 2:
                    fails.append((c, t))
                else:
                    time.sleep(2 * (k + 1))
        if got is not None:
            json.dump(got, open(fp, "w"))
            rows += got
        if i % 300 == 0:
            log(f"    {i}/{len(tiles)} tiles | {len(rows):,} rows | {len(fails)} fails")
    log(f"    done {len(tiles)} tiles | {len(rows):,} rows | {len(fails)} fails")
    if fails:
        log(f"    WARNING {len(fails)} tiles failed after 3 attempts; rerun to retry "
            f"(cache makes it cheap)")
    return pd.DataFrame(rows).drop_duplicates(subset=["name", "chrom", "txStart"])


def _parse_exons(v):
    if v is None:
        return []
    if isinstance(v, (list, np.ndarray)):
        return [int(z) for z in v if str(z) != ""]
    return [int(z) for z in str(v).strip(",").split(",") if z != ""]


def _window_join(need: pd.DataFrame, G: pd.DataFrame) -> dict:
    """Transcripts overlapping each locus window. G must already have list exons."""
    bych = defaultdict(list)
    for r in G.itertuples():
        bych[r.chrom].append((r.txStart, r.txEnd, r.name, r.name2, r.strand,
                              r.exonStarts, r.exonEnds))
    for c in bych:
        bych[c].sort()

    gm = defaultdict(list)
    for r in need.itertuples():
        arr = bych.get(r.chrom)
        if not arr:
            continue
        for tx in arr:
            if tx[0] > r.w1:
                break          # sorted by txStart, nothing further can overlap
            if tx[1] <= r.w0:
                continue
            gm[r.locus_uid].append({"name": tx[2], "name2": tx[3], "strand": tx[4],
                                    "txStart": tx[0], "txEnd": tx[1],
                                    "exonStarts": list(tx[5]), "exonEnds": list(tx[6])})
    return gm


def _agreement(con, gm: dict, genome: str, source: str, label: str, floor: float = 90.0):
    """Regression guard: do the built models recover hit_gene's gene names?"""
    hgg = pd.read_sql("SELECT locus_uid,ref_gene_name FROM hit_gene "
                      "WHERE genome=? AND source=?", con, params=[genome, source])
    want = {u: set(d.ref_gene_name.dropna()) for u, d in hgg.groupby("locus_uid")}
    a = p = z = 0
    for u, w in want.items():
        i = w & {t["name2"] for t in gm.get(u, ())}
        a += (i == w); p += bool(i) and i != w; z += not i
    tot = a + p + z
    if not tot:
        return True
    log(f"  hit_gene {genome} agreement (n={tot:,}): all {100*a/tot:.1f}% | "
        f"partial {100*p/tot:.1f}% | none {100*z/tot:.1f}%")
    if 100 * a / tot < floor:
        log(f"  WARNING {genome} agreement below {floor:.0f}% -- is {label} the "
            f"annotation hit_gene was computed against?")
        return False
    return True


def build_gene_models(con, out: str, track: str, cache: str, parquet: str,
                      hs1_gtf: str, hs1_parquet: str):
    """
    Build the assembly-keyed gene bundle: {"hg38": {uid: [tx]}, "t2t": {uid: [tx]}}.

    The two assemblies come from DIFFERENT annotation sets, and must:
      hg38 -> GENCODE (UCSC API, `track`), because hit_gene's hg38 rows are gencode
      t2t  -> hs1.ncbiRefSeq.gtf.gz,       because hit_gene's t2t rows are refseq
    There is no GENCODE for hs1 at UCSC, and using a different source for either
    reproduces the V44/V50 disagreement this script's docstring warns about.
    """
    bundle, ok = {}, True

    # ---- hg38 / GENCODE
    need = gene_windows(con, "hg38")
    log(f"  hg38: loci needing models {len(need):,}")
    if os.path.exists(parquet):
        G = pd.read_parquet(parquet)
        log(f"  reusing {parquet} ({len(G):,} transcripts)")
    else:
        G = fetch_gencode(need, track, cache)
        G.to_parquet(parquet, index=False)
        log(f"  wrote {parquet} ({len(G):,} transcripts)")
    G["exonStarts"] = G.exonStarts.map(_parse_exons)
    G["exonEnds"] = G.exonEnds.map(_parse_exons)
    G = G[G.exonStarts.map(len) > 0]
    gm = _window_join(need, G)
    log(f"  hg38: loci with >=1 transcript {len(gm):,} ({100*len(gm)/len(need):.1f}%)")
    bundle["hg38"] = gm
    ok &= _agreement(con, gm, "hg38", "gencode", track)

    # ---- t2t / RefSeq
    needt = gene_windows(con, "t2t")
    log(f"  t2t: loci needing models {len(needt):,}")
    if os.path.exists(hs1_parquet):
        T = pd.read_parquet(hs1_parquet)
        log(f"  reusing {hs1_parquet} ({len(T):,} transcripts)")
    elif os.path.exists(hs1_gtf):
        import parse_refseq_gtf
        T = parse_refseq_gtf.parse_gtf(hs1_gtf)
        T.to_parquet(hs1_parquet, index=False)
        log(f"  wrote {hs1_parquet} ({len(T):,} transcripts)")
    else:
        log(f"  SKIP t2t gene lane: neither {hs1_parquet} nor {hs1_gtf} present")
        log(f"       fetch: curl -O https://hgdownload.soe.ucsc.edu/goldenPath/"
            f"hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz")
        T = None
    if T is not None:
        gmt = _window_join(needt, T)
        log(f"  t2t: loci with >=1 transcript {len(gmt):,} ({100*len(gmt)/len(needt):.1f}%)")
        bundle["t2t"] = gmt
        ok &= _agreement(con, gmt, "t2t", "refseq", os.path.basename(hs1_gtf))

    n = dump_gz(bundle, f"{out}/data/gene_models.json.gz")
    log(f"  gene_models.json.gz  {n/1e6:.2f} MB  "
        f"(assembly-keyed: {', '.join(bundle)})")
    if not ok:
        log("  WARNING gene stage completed with agreement warnings above")


# ---------------------------------------------------------------- stage 4

def repeats_for(con, uid_order):
    """locus_repeat rows keyed by locus_uid, both assemblies.

    Rendered selectively: the graphic draws INTERSPERSED by default and hides
    LOW_INFO behind a toggle. Everything is shipped -- the filter lives in the
    page, so changing your mind costs a page edit, not a rebuild.
    """
    # locus_repeat is OPTIONAL: it is derived from UCSC rmsk rather than from the
    # crosswalk build, so a catalogue can legitimately predate it. Aborting the
    # whole bundle over a missing accessory lane is the wrong failure mode -- the
    # v0.3 build only succeeded because a workspace copy happened to carry the
    # table, which made the saved catalogue non-reproducing. Warn and skip.
    try:
        r = pd.read_sql(
            "SELECT locus_uid,assembly,chrom,start,end,strand,rep_name,rep_class,"
            "rep_family,pct_div,n_loci FROM locus_repeat", con)
    except Exception:
        log("  repeats              locus_repeat absent -- RepeatMasker lane omitted")
        return {}
    by = {}
    for t in r.itertuples(index=False):
        by.setdefault(t.locus_uid, []).append({
            "assembly": t.assembly, "start": int(t.start), "end": int(t.end),
            "strand": t.strand, "rep_name": t.rep_name, "rep_class": t.rep_class,
            "rep_family": t.rep_family, "pct_div": float(t.pct_div),
            "n_loci": int(t.n_loci)})
    for v in by.values():
        v.sort(key=lambda d: (d["assembly"], d["start"]))
    return by


def validate(out: str, debug: bool = False) -> bool:
    ok = True
    idx = json.load(gzip.open(f"{out}/data/search_index.json.gz", "rt"))
    nb = json.load(open(f"{out}/data/shard_meta.json"))["n_buckets"]
    cache = {}

    def shard(b):
        if b not in cache:
            cache[b] = json.load(gzip.open(f"{out}/data/loci/{b}.json.gz", "rt"))
        return cache[b]

    # every dump is legal JSON under a strict parser (browsers are strict)
    strict = lambda c: (_ for _ in ()).throw(ValueError(f"non-finite literal {c!r}"))
    bad = []
    files = ([f"{out}/data/search_index.json.gz", f"{out}/data/gene_models.json.gz"]
             + [f"{out}/data/loci/{b}.json.gz" for b in range(nb)
                if os.path.exists(f"{out}/data/loci/{b}.json.gz")])
    for f in files:
        try:
            json.loads(gzip.open(f, "rb").read().decode(), parse_constant=strict)
        except Exception as e:
            bad.append((os.path.basename(f), str(e)[:80]))
    if bad:
        ok = False
        log(f"  FAIL {len(bad)} files are not strict JSON: {bad[:3]}")
    else:
        log(f"  strict JSON: {len(files)} files clean")

    # EVERY identifier must resolve to a locus actually present in its shard.
    # Exhaustive, not sampled: a strided sample missed a deliberately deleted
    # locus during testing, which is precisely the failure this check is for.
    present = {}
    for b in range(nb):
        if os.path.exists(f"{out}/data/loci/{b}.json.gz"):
            for u in shard(b):
                present[u] = b
    tested = miss = 0
    examples = []
    fz = idx["fuzzy"]
    uids = fz["uids"]
    for t in idx["meta"]["class_types"]:
        for k, v in idx["classifier"][t].items():
            for u in v:
                tested += 1
                b = djb2(u) % nb
                if present.get(u) != b:
                    miss += 1
                    if len(examples) < 3:
                        examples.append((t, k, u, f"in shard {present.get(u)}, want {b}"))
    # fuzzy postings carry uid INDICES into fz["uids"] -- an off-by-one here would
    # silently point every hit at the wrong locus, so resolve through the index
    # exactly as the page does rather than trusting the uid list.
    for key, plist in zip(fz["keys"], fz["post"]):
        for ui, ti, orig, cur in plist:
            tested += 1
            u = uids[ui]
            b = djb2(u) % nb
            if present.get(u) != b:
                miss += 1
                if len(examples) < 3:
                    examples.append(("fuzzy", key, u, f"in shard {present.get(u)}, want {b}"))
    # the collapsed tier must reference real key indices
    nk = len(fz["keys"])
    badc = [c for c, ids in zip(fz["ckeys"], fz["cmap"])
            if any(i < 0 or i >= nk for i in ids)]
    if badc:
        ok = False
        log(f"  FAIL {len(badc)} collapsed keys reference out-of-range postings: {badc[:3]}")
    log(f"  resolutions tested {tested:,} (exhaustive) | misses {miss}"
        + (f" {examples}" if miss else ""))
    if miss:
        ok = False
        log("  FAIL a resolvable identifier points at a locus absent from its shard. "
            "Either the locus is missing, or the shard hash disagrees with djb2 "
            "-- check that index.html was not rebuilt with Python's hash().")

    total = sum(len(shard(b)) for b in range(nb)
                if os.path.exists(f"{out}/data/loci/{b}.json.gz"))
    if total != idx["meta"]["n_loci"]:
        ok = False
        log(f"  FAIL shards hold {total:,} loci, index claims {idx['meta']['n_loci']:,}")
    else:
        log(f"  locus count consistent: {total:,}")

    # graphic coverage: hg38 preferred, t2t fallback. detail.js draws whenever EITHER
    # exists, so counting hg38 only (as this did before the t2t lane) understates it.
    n_hg = n_t2 = n_none = n_rep = n_repasm = 0
    n_pjx = n_arc = n_both = 0
    for b in range(nb):
        if not os.path.exists(f"{out}/data/loci/{b}.json.gz"):
            continue
        for d in shard(b).values():
            asms = {c["assembly"] for c in d["coord"]}
            if "hg38" in asms:
                n_hg += 1
            elif "t2t" in asms:
                n_t2 += 1
            else:
                n_none += 1
            reps = d.get("repeats", [])
            if reps:
                n_rep += 1
                # every repeat must name an assembly the locus actually has a
                # coordinate on, else it would be drawn against the wrong window
                if not {r["assembly"] for r in reps} <= asms:
                    n_repasm += 1
            has_p = bool((d.get("pjx") or {}).get("w"))
            has_a = bool((d.get("arcs") or {}).get("jx"))
            n_pjx += has_p
            n_arc += has_a
            n_both += has_p and has_a
    log(f"  graphic: hg38 {n_hg:,} | t2t fallback {n_t2:,} | none {n_none:,}")
    log(f"  repeats present for {n_rep:,} loci")

    # Junction payloads. A --debug-local build legitimately carries both (exact
    # records for inspection, packed for testing the shipped path). A PUBLIC build
    # carrying `arcs` means the debug gate regressed and the packed layer has become
    # pure overhead -- the specific bug this check exists to catch, since the page
    # renders correctly either way and nothing else would notice the size.
    log(f"  junction payload: packed {n_pjx:,} loci | debug arcs {n_arc:,} loci"
        + ("  (debug build -- both expected)" if debug else ""))
    if n_arc and not debug:
        ok = False
        log(f"  FAIL {n_arc:,} loci ship debug `arcs` in a non-debug build -- the "
            f"--debug-local gate has regressed; packed layer is redundant overhead")
    if n_repasm:
        ok = False
        log(f"  FAIL {n_repasm:,} loci carry repeats for an assembly they have no "
            f"coordinate on -- window would be wrong")

    gmf = f"{out}/data/gene_models.json.gz"
    if os.path.exists(gmf):
        gm = json.load(gzip.open(gmf, "rt"))
        # the bundle is assembly-keyed: {"hg38": {...}, "t2t": {...}}.  A flat
        # {uid: [...]} bundle is the pre-two-assembly shape and detail.js will
        # find no genes at all with it, so fail loudly rather than ship it.
        if not (isinstance(gm, dict) and gm and
                all(k in ("hg38", "t2t") for k in gm)):
            ok = False
            log("  FAIL gene_models.json.gz is not assembly-keyed "
                "(expected top-level 'hg38'/'t2t' keys)")
        else:
            log("  gene models: " + " | ".join(
                f"{k} {len(v):,} loci" for k, v in gm.items()))
            for k, v in gm.items():
                bad = [u for u, txs in list(v.items())[:2000]
                       if any(len(t["exonStarts"]) != len(t["exonEnds"]) for t in txs)]
                if bad:
                    ok = False
                    log(f"  FAIL {k}: {len(bad)} loci have exonStarts/exonEnds "
                        f"length mismatch, e.g. {bad[:3]}")
    return ok


# ---------------------------------------------------------------- driver

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="herv_db/herv_catalog.db")
    ap.add_argument("--out", default="dash")
    ap.add_argument("--buckets", type=int, default=N_BUCKETS,
                    help=f"shard count (default {N_BUCKETS}); must match nothing in "
                         f"index.html -- it reads shard_meta.json")
    ap.add_argument("--gencode", default=DEFAULT_TRACK,
                    help=f"UCSC track (default {DEFAULT_TRACK}; see docstring before changing)")
    ap.add_argument("--gene-cache", default="genecache50", help="per-tile fetch cache dir")
    ap.add_argument("--hs1-gtf", default="hs1genes/hs1.ncbiRefSeq.gtf.gz",
                    help="UCSC hs1 RefSeq GTF for the t2t gene lane")
    ap.add_argument("--hs1-parquet", default="hs1_refseq_models.parquet",
                    help="parsed hs1 transcript cache (reused if present)")
    ap.add_argument("--gene-parquet", default="gencode_v50_models.parquet",
                    help="reused if present, else written after fetch")
    ap.add_argument("--tu-db", default="herv_tu_v0.1.db",
                    help="TU/transcription database; panels are skipped if absent")
    # The default was "snaptron_srav3h_jx_sc10.parquet" -- a filename that has
    # never existed. Every build that did not pass --jx-parquet explicitly
    # therefore resolved it to nothing, logged one tolerant line, and shipped a
    # bundle with no arcs. v0.7 did exactly that. The default now names the file
    # that is actually on disk, and OPTIONAL_INPUTS below makes a miss loud.
    ap.add_argument("--jx-parquet", default="snaptron_srav3h_herv_junctions.parquet",
                    help="Snaptron srav3h canonical HERV junctions (hg38); "
                         "supplies the packed-arc layer")
    ap.add_argument("--f5-parquet", default="fantom5_peaks_herv.parquet",
                    help="FANTOM5 CAGE dominant-TSS positions for HERV windows")
    ap.add_argument("--iv-parquet", default="mappability_intervals.parquet",
                    help="standalone merged mappability interval reference; supplies "
                         "the drawn blocks for the mappability lane (stats come "
                         "from the catalogue and ship without it)")
    ap.add_argument("--debug-local", action="store_true",
                    help="also write the per-TU shard set and TU detail route "
                         "(stage 5). Roughly doubles bundle size -- intended for "
                         "local/internal use, not the size-capped public host.")
    ap.add_argument("--skip-index", action="store_true")
    ap.add_argument("--skip-shards", action="store_true")
    ap.add_argument("--skip-genes", action="store_true")
    ap.add_argument("--validate-only", action="store_true")
    ap.add_argument("--allow-incomplete", action="store_true",
                    help="downgrade CATALOG_MANIFEST shortfalls from fatal to "
                         "warnings (for deliberate builds from a leaner "
                         "catalog). The bundle is stamped incomplete.")
    ap.add_argument("--allow-missing-layers", action="store_true",
                    help="downgrade missing evidence-layer input files "
                         "(OPTIONAL_INPUTS) from fatal to warnings, for a "
                         "deliberately lean bundle. Stamped in metadata.")
    a = ap.parse_args(argv)

    if not os.path.exists(a.db):
        sys.exit(f"database not found: {a.db}")
    os.makedirs(f"{a.out}/data/loci", exist_ok=True)
    con = sqlite3.connect(a.db)

    if a.validate_only:
        log("[validate]")
        return 0 if validate(a.out, debug=a.debug_local) else 1

    log("[0/4] catalog manifest")
    manifest_report = check_catalog(con, allow_incomplete=a.allow_incomplete)
    # File-level twin of the catalog manifest. Skipped when the stages that
    # consume these files are not running, so --skip-shards stays usable.
    input_report = ({"inputs": [], "complete": None, "skipped": "no shard stage"}
                    if a.skip_shards else
                    check_inputs(a, allow_missing=a.allow_missing_layers))

    t0 = time.time()
    loc = al = None
    if not a.skip_index:
        log("[1/4] search index")
        loc, al = build_search_index(con, a.out)
    if not a.skip_shards:
        log("[2/4] locus shards")
        if loc is None:
            loc = pd.read_sql('SELECT locus_uid,combined_id,versioned_id,"group",band,'
                              'origin FROM locus', con)
            al = pd.read_sql("SELECT locus_uid,alias,alias_type,assignment,is_current "
                             "FROM locus_alias", con)
        build_shards(con, a.out, loc, al, a.buckets, a.tu_db, a.jx_parquet,
                     debug=a.debug_local, iv_parquet=a.iv_parquet,
                     f5_parquet=a.f5_parquet)
    if a.debug_local:
        log("[5/5] TU shards (--debug-local)")
        build_tu_shards(con, a.out, a.tu_db, a.jx_parquet, a.buckets)
    if not a.skip_genes:
        log(f"[3/4] gene models ({a.gencode})")
        build_gene_models(con, a.out, a.gencode, a.gene_cache, a.gene_parquet,
                          a.hs1_gtf, a.hs1_parquet)
    # Provenance: record what the catalog contained at build time, so an
    # incomplete bundle is identifiable without re-deriving it from the data.
    dump_json({**manifest_report, "input_manifest": input_report},
              f"{a.out}/data/catalog_manifest.json")
    log("[4/4] validate")
    good = validate(a.out, debug=a.debug_local)
    log(f"\n{'OK' if good else 'FAILED'} in {time.time()-t0:.0f}s -> {a.out}/")
    if good:
        log(f"serve with:  cd {a.out} && python -m http.server 8000")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main())
