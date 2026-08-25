#!/usr/bin/env python3
"""
Transcription-evidence payload for the HERV dashboard: FANTOM5 TSS / splice-site
status per locus and per transcriptional unit, plus Snaptron srav3h junction arcs.

Imported by build_dashboard.py (stage 2 folds the locus payload into the locus
shards; stage 5 writes the TU shards under --debug-local).

--- Semantics this module exists to keep straight ---------------------------------

TRI-STATE, NOT BOOLEAN.  `assessable` is carried separately from evidence counts.
A locus with assessable=False has no FANTOM coverage at all; a locus with
assessable=True and tss_level=None was looked at and found silent. Collapsing
those two into "no evidence" was the mistake that made an earlier draft of this
layer unreadable -- the UI must be able to say "not assessable" and "assessed,
none found" as different things.

LEVEL IS THE STRICTEST TIER STILL SUPPORTED.  tss_level / sj_level are computed by
walking permissive -> robust -> stringent in that order and keeping the last tier
that had sense evidence, so "stringent" implies the looser tiers also fired. Do
not reorder TIERS: it silently inverts the meaning.

DIRECTION IS SENSE-RELATIVE, NOT GENOMIC.  *_anti flags mean evidence on the
strand opposite the element's own strand. The upstream layer resolved element
strand already; nothing here re-derives it from coordinates.

ARC COUNTS ARE PRE-CAP.  n_total (and n_cross / n_within for chimeric units) are
the true counts in the window; `jx` is the capped subset actually shipped. The UI
must render "showing N of M" from both, never len(jx) alone -- a static bundle
cannot carry the 1.43M junctions at samples_count>=10 and pretending otherwise
would misstate the evidence.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

TIERS = ("permissive", "robust", "stringent")
PAD = 1000          # window around an element, bp; matches the FANTOM layer
CAP = 25            # arcs per locus / per non-chimeric unit, ranked by support
CAP_CROSS = 25      # chimeric units: cross-element arcs
CAP_WITHIN = 25     # chimeric units: within-element arcs
SC_FLOOR = 10       # samples_count floor, from the threshold calibration


# ------------------------------------------------------------------ state calls

def _tss_state(r):
    """(strictest tier with sense TSS in body or upstream window, any antisense)."""
    lvl = None
    for t in TIERS:
        if (r.get(f"tss_body_{t}_sense_n") or 0) > 0 or (r.get(f"tss_up_{t}_sense_n") or 0) > 0:
            lvl = t
    anti = any((r.get(f"tss_body_{t}_anti_n") or 0) > 0 for t in TIERS)
    return lvl, anti


def _sj_state(r):
    """(strictest tier with a sense donor or acceptor, any antisense)."""
    lvl = None
    for t in TIERS:
        if any((r.get(f"sj_{w}_{t}_sense_{e}_n") or 0) > 0
               for w in ("body", "up") for e in ("donor", "acceptor")):
            lvl = t
    anti = any((r.get(f"sj_body_{t}_anti_donor_n") or 0) > 0
               or (r.get(f"sj_body_{t}_anti_acceptor_n") or 0) > 0 for t in TIERS)
    return lvl, anti


def _num(v, cast=int):
    return None if v is None or pd.isna(v) else cast(v)


def _row_payload(r):
    """Compact per-assembly record. Counts/distances are reported at the robust
    tier -- the middle rung -- while *_level records how far up the ladder the
    evidence actually reaches."""
    tl, ta = _tss_state(r)
    sl, sa = _sj_state(r)
    return {
        "assessable": bool(r.get("fantom_assessable")),
        "tss_level": tl, "tss_anti": bool(ta),
        "sj_level": sl, "sj_anti": bool(sa),
        "tss_body_n": int(r.get("tss_body_robust_sense_n") or 0),
        "tss_up_n": int(r.get("tss_up_robust_sense_n") or 0),
        "tss_anti_n": int(r.get("tss_body_robust_anti_n") or 0),
        "tss5p_dist": _num(r.get("tss5p_robust_dist")),
        "tss5p_score": _num(r.get("tss5p_robust_score"), float),
        "sj_donor_n": int(r.get("sj_body_robust_sense_donor_n") or 0),
        "sj_acceptor_n": int(r.get("sj_body_robust_sense_acceptor_n") or 0),
        "sj_span_n": int(r.get("sjspan_robust_sense_n") or 0),
    }


def locus_tx(cx) -> dict:
    """{locus_uid: {asm: record}} from locus_transcription."""
    lt = pd.read_sql("SELECT * FROM locus_transcription", cx)
    out = {}
    for r in lt.to_dict("records"):
        out.setdefault(r["locus_uid"], {})[r["asm"]] = _row_payload(r)
    return out


def tu_tx(cx) -> dict:
    """{tu_id: {asm: record}} from tu_transcription."""
    tt = pd.read_sql("SELECT * FROM tu_transcription", cx)
    out = {}
    for r in tt.to_dict("records"):
        out.setdefault(r["tu_id"], {})[r.get("asm") or "hg38"] = _row_payload(r)
    return out


SNAP_COLS = ("jx_internal_n", "jx_boundary_n", "jx_intronic_ctx_n", "jx_n",
             "jx_internal_max_sc", "jx_boundary_max_sc", "jx_intronic_ctx_max_sc",
             "intronic_ctx_median_span", "cross_element_jx_n", "within_element_jx_n",
             "cross_max_sc", "within_max_sc", "n_mem_loci", "n_mem_groups",
             "is_chimeric", "snaptron_assessable")


def tu_snaptron(cx) -> dict:
    """{tu_id: summary} -- the three membership classes, cross/within counts and
    the split verdict. Class names are the audited ones: jx_intronic_ctx_* counts
    junctions whose intron *contains* the element (host-gene introns), which is
    not readthrough of it."""
    sj = pd.read_sql("SELECT * FROM tu_snaptron_junctions", cx)
    out = {}
    for r in sj.to_dict("records"):
        d = {c: _num(r.get(c)) for c in SNAP_COLS}
        d["split_verdict"] = r.get("split_verdict")
        d["group_call"] = r.get("group_call")
        out[r["tu_id"]] = d
    return out


# ------------------------------------------------------------------ arcs

def _load_jx(path: str) -> pd.DataFrame:
    jx = pd.read_parquet(path)
    if "sc" not in jx.columns and "samples_count" in jx.columns:
        jx = jx.rename(columns={"samples_count": "sc"})
    jx = jx[jx.sc >= SC_FLOOR]
    # Plain-overlap tiling returns a spanning junction once per tile it crosses.
    return jx.sort_values("sc", ascending=False).drop_duplicates(
        ["chrom", "start", "end", "strand"])


def _members(cx, cat_con) -> pd.DataFrame:
    """Member elements of each unit with coordinates and group call."""
    lm = pd.read_sql("SELECT locus_uid, assembly AS asm, tu_id FROM tu_locus_member", cx)
    lc = pd.read_sql("SELECT locus_uid, asm, combined_id, chrom, start, end, strand "
                     "FROM locus_transcription", cx)
    m = lm.merge(lc, on=["locus_uid", "asm"], how="left").dropna(
        subset=["chrom", "start", "end"])
    grp = pd.read_sql('SELECT locus_uid, "group" AS grp FROM locus', cat_con)
    m["grp"] = m.locus_uid.map(dict(zip(grp.locus_uid, grp.grp)))
    m["start"] = m.start.astype(int)
    m["end"] = m.end.astype(int)
    return m


def tu_members(cx, cat_con) -> dict:
    """Member-element track per unit, hg38 preferred with a t2t fallback.

    Deliberately independent of tu_arcs: Snaptron srav3h is hg38-only, but 3,669
    units (15%) exist only in t2t coordinates. Deriving the track from the arc
    payload left those units with no map at all -- the same hg38/t2t fallback the
    locus page already does, so the unit page does it too. The returned asm names
    which coordinates were used so the axis label cannot claim the wrong assembly."""
    mem = _members(cx, cat_con)
    out = {}
    for tid, g in mem.groupby("tu_id", sort=False):
        for asm in ("hg38", "t2t"):
            h = g[g.asm == asm]
            if len(h):
                out[tid] = {
                    "asm": asm, "chrom": h.chrom.iloc[0],
                    "members": [[r.locus_uid, r.combined_id or "", int(r.start),
                                 int(r.end), r.strand, r.grp or ""]
                                for r in h.sort_values("start").itertuples()],
                }
                break

    # 3,018 units (12%) carry coordinates in tu_transcription but have NO row in
    # tu_locus_member -- no locus mapping was ever made for them, in either
    # assembly (2,866 are on primary chromosomes, so this is not a scaffold
    # artifact). That is an upstream gap in the v0.1 unit build, not something to
    # paper over: they get the unit extent as a single unmapped bar, flagged
    # no_members, so the page states the mapping is absent rather than implying
    # the unit has no structure.
    tc = pd.read_sql("SELECT tu_id, asm, chrom, start, end, strand FROM tu_transcription",
                     cx).dropna(subset=["chrom", "start", "end"])
    for tid, g in tc[~tc.tu_id.isin(out)].groupby("tu_id", sort=False):
        for asm in ("hg38", "t2t"):
            h = g[g.asm == asm]
            if len(h):
                r = h.iloc[0]
                out[tid] = {"asm": asm, "chrom": r.chrom, "no_members": True,
                            "extent": [int(r.start), int(r.end), r.strand],
                            "members": []}
                break
    return out


def locus_arcs(cx, jx_path: str, asm: str = "hg38") -> dict:
    """{locus_uid: {n_total, shown, jx:[[start,end,sc,strand,canonical],...]}}.

    Snaptron srav3h is hg38-only, so arcs exist for one assembly; the t2t view
    states that rather than showing an empty lane."""
    jx = _load_jx(jx_path)
    lc = pd.read_sql("SELECT locus_uid, asm, chrom, start, end FROM locus_transcription",
                     cx).dropna(subset=["chrom", "start", "end"])
    lc = lc[lc.asm == asm]
    lc["start"] = lc.start.astype(int)
    lc["end"] = lc.end.astype(int)
    out = {}
    for ch, g in lc.groupby("chrom", sort=False):
        jc = jx[jx.chrom == ch]
        if not len(jc):
            continue
        js, je = jc.start.values, jc.end.values
        for r in g.itertuples():
            m = (je >= r.start - PAD) & (js <= r.end + PAD)
            n = int(m.sum())
            if not n:
                continue
            sub = jc[m].nlargest(CAP, "sc")
            out[r.locus_uid] = {
                "n_total": n, "shown": len(sub),
                "jx": [[int(x.start), int(x.end), int(x.sc), x.strand, int(x.canonical)]
                       for x in sub.itertuples()]}
    return out


def tu_arcs(cx, cat_con, jx_path: str, asm: str = "hg38") -> dict:
    """Per-unit arcs plus the member track the TU view draws.

    Chimeric units get a split budget -- up to CAP_CROSS cross-element and
    CAP_WITHIN within-element arcs -- so a general support ranking cannot bury
    the junctions that discriminate between splitting and keeping a unit merged.
    Every arc carries its class so the view can style the two differently."""
    jx = _load_jx(jx_path)
    tc = pd.read_sql("SELECT tu_id, asm, chrom, start, end, strand FROM tu_transcription",
                     cx).dropna(subset=["chrom", "start", "end"])
    tc = tc[tc.asm == asm]
    tc["start"] = tc.start.astype(int)
    tc["end"] = tc.end.astype(int)
    mem = _members(cx, cat_con)
    mem = mem[mem.asm == asm]
    by_tu = {t: g for t, g in mem.groupby("tu_id", sort=False)}
    _cq = pd.read_sql("SELECT tu_id, is_chimeric, cross_element_jx_n, within_element_jx_n "
                      "FROM tu_snaptron_junctions", cx)
    chim = set(_cq[_cq.is_chimeric == 1].tu_id)
    up_counts = {r.tu_id: (r.cross_element_jx_n or 0, r.within_element_jx_n or 0)
                 for r in _cq.dropna(subset=["cross_element_jx_n"]).itertuples()}
    out = {}
    for ch, gT in tc.groupby("chrom", sort=False):
        jc = jx[jx.chrom == ch]
        if not len(jc):
            continue
        js, je = jc.start.values, jc.end.values
        for r in gT.itertuples():
            lo, hi = r.start - PAD, r.end + PAD
            m = (je >= lo) & (js <= hi)
            n = int(m.sum())
            if not n:
                continue
            sub = jc[m]
            g = by_tu.get(r.tu_id)
            rec = {"n_total": n, "strand": r.strand, "win": [int(lo), int(hi)],
                   "chrom": ch,
                   "members": ([[x.locus_uid, x.combined_id, int(x.start), int(x.end),
                                 x.strand, x.grp] for x in g.itertuples()]
                               if g is not None else [])}
            if r.tu_id in chim and g is not None and len(g) > 1:
                # Classification is GROUP-based and order-independent, and both
                # properties are load-bearing.
                #
                # Member elements can be NESTED -- e.g. ERV1_3q25.1e has an
                # HERVFH19 internal segment lying inside an LTR23 span. An earlier
                # version resolved each junction end to the FIRST member row that
                # contained it, which made the call depend on row order: the same
                # junction read as cross or within depending on which nested member
                # was listed first. That produced arcs drawn "cross" on a unit whose
                # upstream summary said zero cross -- the page contradicting itself.
                #
                # Correct rule: collapse members to their GROUP, ask which groups
                # contain each end, and call the junction within-element when the
                # two ends share any group and cross-element only when they share
                # none. Reproduces the upstream cross counts exactly on the
                # hand-checked units and raises agreement over the element-index
                # rule (cross 0.888 vs 0.857 across 2,376 chimeric units).
                gl = sorted(set(g.grp.fillna("")))
                gsegs = [[(int(x.start), int(x.end))
                          for x in g[g.grp.fillna("") == k].itertuples()] for k in gl]

                def ghit(p, gsegs=gsegs):
                    return {i for i, sp in enumerate(gsegs)
                            if any(s <= p <= e for s, e in sp)}

                ha = [ghit(p) for p in sub.start.values]
                hb = [ghit(p) for p in sub.end.values]
                both = np.array([bool(a) and bool(b) for a, b in zip(ha, hb)])
                shared = np.array([bool(a & b) for a, b in zip(ha, hb)])
                cross = sub[both & ~shared]
                within = sub[both & shared]
                rec["n_cross"], rec["n_within"] = len(cross), len(within)
                # Pre-cap anchored total. On a chimeric unit the shown set is drawn
                # from the anchored junctions only, so shown < n_total has two very
                # different causes: the cap bit (shown == budget), or most window
                # junctions were unanchorable (n_anchored << n_total). Recording
                # both lets the panel say which -- "top 25 by support" is false
                # advertising for a unit with 195 in window and 2 anchored.
                rec["n_anchored"] = len(cross) + len(within)
                # The panel's summary counts come from the upstream Snaptron layer;
                # these arcs are recomputed here from the junction parquet. They
                # agree on most units but not all (the upstream derivation is not
                # recoverable in full -- window padding and any support/motif filter
                # it applied are unknown). Where they differ, the page must say so
                # rather than showing two numbers and letting the reader assume they
                # are the same quantity.
                up = up_counts.get(r.tu_id)
                if up is not None and (up[0] != len(cross) or up[1] != len(within)):
                    rec["recount"] = [int(up[0]), int(up[1])]
                rec["jx"] = (
                    [[int(x.start), int(x.end), int(x.sc), x.strand, int(x.canonical), "cross"]
                     for x in cross.nlargest(CAP_CROSS, "sc").itertuples()] +
                    [[int(x.start), int(x.end), int(x.sc), x.strand, int(x.canonical), "within"]
                     for x in within.nlargest(CAP_WITHIN, "sc").itertuples()])
                # 983 chimeric units have junctions in the padded window but none
                # with BOTH ends inside a member element -- they start or land in
                # the flanks. Both budgets come up empty and the panel would be
                # blank, which reads as "no junction data" when the truth is
                # "junctions here, none of them element-to-element". Show the
                # window's top arcs, classed "unanchored", and let n_cross=0 say
                # the split-relevant class is genuinely empty.
                if not rec["jx"]:
                    rec["jx"] = [[int(x.start), int(x.end), int(x.sc), x.strand,
                                  int(x.canonical), "unanchored"]
                                 for x in sub.nlargest(CAP, "sc").itertuples()]
                    rec["unanchored"] = True
            else:
                rec["jx"] = [[int(x.start), int(x.end), int(x.sc), x.strand,
                              int(x.canonical), "any"]
                             for x in sub.nlargest(CAP, "sc").itertuples()]
            rec["shown"] = len(rec["jx"])
            out[r.tu_id] = rec
    return out
