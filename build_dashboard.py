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

IDENT = ["combined_id", "versioned_id", "telescope_id", "hervd_id",
         "ervmap_id", "ervmap_alt_name"]
CLASS = ["dfam_int_model", "dfam_accession", "repbase_name"]
DETAIL_TABLES = ["locus_coord", "locus_segment", "locus_structure", "hit_geve",
                 "hit_hervarium_domain", "hit_hervarium_int", "hit_gene",
                 "aln_crossgenome"]
PAD = 1000          # graphic window padding, must match PAD in detail.js
N_BUCKETS = 400
TILE = 2_000_000    # UCSC query tile size
DEFAULT_TRACK = "wgEncodeGencodeCompV50"


# ---------------------------------------------------------------- primitives

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
    """Write gzipped JSON. allow_nan=False turns NaN leakage into an exception."""
    with gzip.open(path, "wt") as fh:
        json.dump(_clean(obj), fh, separators=(",", ":"), allow_nan=False, default=str)
    return os.path.getsize(path)


def dump_json(obj, path: str) -> int:
    with open(path, "w") as fh:
        json.dump(_clean(obj), fh, separators=(",", ":"), allow_nan=False, default=str)
    return os.path.getsize(path)


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

    n = dump_gz({"classifier": clsidx, "fuzzy": fuzzy, "listmeta": listmeta,
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


# ---------------------------------------------------------------- stage 2

def build_shards(con, out: str, loc: pd.DataFrame, al: pd.DataFrame, nb: int):
    grp = pd.read_sql('SELECT * FROM "group"', con).set_index("group")
    sf = pd.read_sql("SELECT * FROM superfamily", con).set_index("superfamily")
    tabs = {t: pd.read_sql(f"SELECT * FROM {t}", con) for t in DETAIL_TABLES}
    dbest = pd.read_sql("SELECT locus_uid,dfam_accession,consensus_name,pct_identity,"
                        "cons_cov,sw_score,aln_quality,is_best_by_sw_score "
                        "FROM hit_dfam_aln WHERE is_best_by_sw_score=1", con)

    byuid = {k: {u: d for u, d in v.groupby("locus_uid")} for k, v in tabs.items()}
    dby = {u: d for u, d in dbest.groupby("locus_uid")}
    alby = {u: d for u, d in al.groupby("locus_uid")}
    rby = repeats_for(con, loc.locus_uid)

    GKEYS = ("superfamily", "herv_class", "n_loci", "intModel", "repbase_class",
             "hervd_family", "dfam_accession", "dominant_ltr",
             "frac_with_flanking_ltr", "n_with_hervarium_domain", "extension_verdict")

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
            "crossgenome": (recs(byuid["aln_crossgenome"].get(u)) or [{}])[0],
            "dfam_best": recs(dby.get(u)),
            "repeats": rby.get(u, []),
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
    r = pd.read_sql(
        "SELECT locus_uid,assembly,chrom,start,end,strand,rep_name,rep_class,"
        "rep_family,pct_div,n_loci FROM locus_repeat", con)
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


def validate(out: str) -> bool:
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
    log(f"  graphic: hg38 {n_hg:,} | t2t fallback {n_t2:,} | none {n_none:,}")
    log(f"  repeats present for {n_rep:,} loci")
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
    ap.add_argument("--skip-index", action="store_true")
    ap.add_argument("--skip-shards", action="store_true")
    ap.add_argument("--skip-genes", action="store_true")
    ap.add_argument("--validate-only", action="store_true")
    a = ap.parse_args(argv)

    if not os.path.exists(a.db):
        sys.exit(f"database not found: {a.db}")
    os.makedirs(f"{a.out}/data/loci", exist_ok=True)
    con = sqlite3.connect(a.db)

    if a.validate_only:
        log("[validate]")
        return 0 if validate(a.out) else 1

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
        build_shards(con, a.out, loc, al, a.buckets)
    if not a.skip_genes:
        log(f"[3/4] gene models ({a.gencode})")
        build_gene_models(con, a.out, a.gencode, a.gene_cache, a.gene_parquet,
                          a.hs1_gtf, a.hs1_parquet)
    log("[4/4] validate")
    good = validate(a.out)
    log(f"\n{'OK' if good else 'FAILED'} in {time.time()-t0:.0f}s -> {a.out}/")
    if good:
        log(f"serve with:  cd {a.out} && python -m http.server 8000")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main())
