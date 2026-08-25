"""Bit-packed splice-junction reference for the non-debug dashboard.

One junction := one 64-bit big-endian word, so a locus's whole splicing
evidence is a fixed-stride byte range and the browser needs no per-junction
object.  Layout, MSB first:

    bits 63..42  (22)  donor offset, SIGNED, saturating
    bit  41      ( 1)  strand-mismatch flag: 1 = junction strand != locus strand
    bits 40..19  (22)  acceptor offset, SIGNED, saturating
    bit  18      ( 1)  anchor-mappable flag: 1 = at least one outer 50 bp anchor
                       window lies wholly inside a pm151 Panmask "easy" interval
    bits 17..16  ( 2)  SPARE (zero on write; a reader must ignore them)
    bits 15..8   ( 8)  sample count,  round(log10(sc)  * 32)   [3+5 fixed point]
    bits  7..0   ( 8)  coverage/sample, round(log10(cov) * 32)  [3+5 fixed point]

Both offsets are relative to the locus 5' end and sign-flipped on minus-strand
loci, so a negative donor offset always means "upstream of the element" no
matter which strand it sits on.  That is what makes the packed value
interpretable without also shipping the locus coordinates.

Why both offset fields are 22 bits
----------------------------------
The original layout was 3 bytes per offset, then 23+24 once the strand flag
needed a bit.  Neither was measured against the data that actually gets packed.
Over the real kept set (172,678 words), |donor offset| tops out at 1,392,836 and
|acceptor offset| at 1,187,155 -- both inside a 22-bit signed capacity of
2,097,151, with 33% headroom, and zero saturation.  21 bits does NOT suffice
(+/-1,048,575 overflows 8 words), so 22 is the floor, not a guess.

Narrowing 23+24 -> 22+22 reclaims 3 bits: one for the anchor-mappable flag and
two left spare.  This is a BREAKING format change -- every field below bit 63
shifted -- but the decoder ships inside the same bundle as the data, so there is
no version-skew window.  Saturation remains the documented behaviour for
out-of-range values rather than an error, and `encode_locus_junctions` returns
the saturation count so a build can assert it stays at zero.

Why the anchor flag is pm151 and not Umap
-----------------------------------------
The flag answers "could a split read be placed here at all", which is what
Panmask's easy regions are calibrated for at 151 bp -- the read length of the
srav3h junction set.  Umap k100 single-read uniqueness is stored per LOCUS in
the catalog instead; it is a per-base uniqueness question, and being a shorter
k-mer it understates mappability relative to real 151 bp reads.  Measured over
all 1,425,453 distinct junctions, zero have an anchor window on a contig pm151
does not cover, so a CLEAR bit unambiguously means "neither anchor is easy" and
never "no data".

Why coverage is log-scaled rather than uint8
--------------------------------------------
cov_per_sample runs 1.0 -> 49,053 on the kept set with a median of 1.5.  A
plain uint8 both clips the deep tail and rounds away the fractional low end
where most junctions sit, collapsing 1.0-1.4 onto 1.  The same 3+5 log10 fixed
point used for the sample count covers the full range at ~2% relative error and
shares one decode path with it.

The class of a junction (spanning / both ends inside / donor only / acceptor
only) is NOT packed: it is recoverable from the two offsets and the locus
length, which the renderer already has.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

SC_FLOOR = 10          # samples_count floor, from the threshold calibration
PAD = 1000             # window around the locus, in bp
CAP_INSIDE = 15        # top-N by sample count with donor and/or acceptor inside
CAP_SPAN = 1           # top-N by sample count spanning the whole locus
ANCHOR_WIN = 50        # bp of exonic flank tested for mappability at each end

LOG_SCALE = 32         # 3+5 fixed point: value = round(log10(x) * 32)
DONOR_BITS = 22
ACCEPTOR_BITS = 22
_SHIFT_DONOR = 42
_SHIFT_STRAND = 41
_SHIFT_ACC = 19
_SHIFT_ANCHOR = 18
_DONOR_MIN, _DONOR_MAX = -(2 ** (DONOR_BITS - 1)), 2 ** (DONOR_BITS - 1) - 1
_ACC_MIN, _ACC_MAX = -(2 ** (ACCEPTOR_BITS - 1)), 2 ** (ACCEPTOR_BITS - 1) - 1
_DONOR_MASK = np.uint64((1 << DONOR_BITS) - 1)
_ACC_MASK = np.uint64((1 << ACCEPTOR_BITS) - 1)


def log_code(x) -> np.ndarray:
    """round(log10(x) * 32), clamped to a byte.  x <= 0 encodes as 0."""
    v = np.asarray(x, dtype=np.float64)
    out = np.zeros(v.shape, dtype=np.float64)
    pos = v > 0
    out[pos] = np.log10(v[pos]) * LOG_SCALE
    return np.clip(np.round(out), 0, 255).astype(np.uint8)


def log_decode(code) -> np.ndarray:
    """Inverse of log_code, for verification and for the JS reference test."""
    return 10.0 ** (np.asarray(code, dtype=np.float64) / LOG_SCALE)


def pack(off_d, off_a, antisense, sc_code, cov_code, anchor_ok=0) -> np.ndarray:
    """Pack parallel arrays into uint64 words.  No clipping is done here.

    `anchor_ok` defaults to 0 so a caller that has no mappability resource still
    produces well-formed words -- but a zero there is indistinguishable from a
    measured "neither anchor easy", so `build` always passes a real value.
    """
    d = (np.asarray(off_d, np.int64) & np.int64(_DONOR_MASK)).astype(np.uint64)
    a = (np.asarray(off_a, np.int64) & np.int64(_ACC_MASK)).astype(np.uint64)
    f = np.asarray(antisense, np.uint64) & np.uint64(1)
    m = np.asarray(anchor_ok, np.uint64) & np.uint64(1)
    return ((d << np.uint64(_SHIFT_DONOR)) | (f << np.uint64(_SHIFT_STRAND))
            | (a << np.uint64(_SHIFT_ACC)) | (m << np.uint64(_SHIFT_ANCHOR))
            | (np.asarray(sc_code, np.uint64) << np.uint64(8))
            | np.asarray(cov_code, np.uint64))


def unpack(words) -> dict:
    """Inverse of pack.  Returns signed offsets, flag, and decoded values."""
    w = np.asarray(words, np.uint64)
    d = (w >> np.uint64(_SHIFT_DONOR)) & _DONOR_MASK
    a = (w >> np.uint64(_SHIFT_ACC)) & _ACC_MASK
    d = d.astype(np.int64)
    a = a.astype(np.int64)
    d = np.where(d > _DONOR_MAX, d - (1 << DONOR_BITS), d)
    a = np.where(a > _ACC_MAX, a - (1 << ACCEPTOR_BITS), a)
    sc_code = ((w >> np.uint64(8)) & np.uint64(0xFF)).astype(np.uint8)
    cov_code = (w & np.uint64(0xFF)).astype(np.uint8)
    return {"off_d": d, "off_a": a,
            "antisense": ((w >> np.uint64(_SHIFT_STRAND)) & np.uint64(1)).astype(np.uint8),
            "anchor_ok": ((w >> np.uint64(_SHIFT_ANCHOR)) & np.uint64(1)).astype(np.uint8),
            "sc_code": sc_code, "cov_code": cov_code,
            "sc": log_decode(sc_code), "cov": log_decode(cov_code)}


def _assign(jx: pd.DataFrame, loci: pd.DataFrame) -> pd.DataFrame:
    """Every (locus, junction) pair where the junction is in the padded window
    or spans the locus outright.  Spanning pairs are collected separately
    because a junction that starts far upstream is not found by a window
    search on its start coordinate."""
    out = []
    for ch, gl in loci.groupby("chrom", sort=False):
        gj = jx[jx.chrom == ch]
        if not len(gj):
            continue
        gj = gj.sort_values("start")
        js, je = gj.start.values, gj.end.values
        sc, cov, jst = gj.sc.values, gj.cov_per_sample.values, gj.strand.values
        for r in gl.itertuples():
            lo, hi = r.start - PAD, r.end + PAD
            i, k = np.searchsorted(js, lo), np.searchsorted(js, hi, side="right")
            idx = np.arange(i, k)[je[i:k] >= lo]
            spanning = np.where((js < r.start) & (je > r.end))[0]
            idx = np.union1d(idx, spanning)
            if not len(idx):
                continue
            out.append(pd.DataFrame({
                "chrom": ch,
                "locus_uid": r.locus_uid, "lstart": r.start, "lend": r.end,
                "lstrand": r.strand, "jstart": js[idx], "jend": je[idx],
                "sc": sc[idx], "cov": cov[idx], "jstrand": jst[idx]}))
    if not out:
        return pd.DataFrame()
    return pd.concat(out, ignore_index=True)


def select(pairs: pd.DataFrame) -> pd.DataFrame:
    """Apply the per-locus caps: top CAP_SPAN spanning junctions by sample
    count, plus top CAP_INSIDE with donor and/or acceptor inside the locus.
    Junctions with both ends outside and not spanning are dropped."""
    p = pairs.copy()
    # Donor/acceptor are assigned from the JUNCTION's own strand, not from
    # genomic order and not from the locus strand: the donor is the junction's
    # 5' end, which on a minus-strand junction is the HIGHER coordinate.
    #
    # This previously tested jstart as the donor unconditionally. That is a
    # coordinate label masquerading as a biological one, and it disagreed with
    # the strand-aware assignment on 46,486 of 172,678 kept junctions (27%) --
    # entirely within the donor-vs-acceptor distinction, so `both` and `span`
    # were unaffected and nothing downstream looked wrong.
    lo_in = (p.jstart >= p.lstart) & (p.jstart <= p.lend)
    hi_in = (p.jend >= p.lstart) & (p.jend <= p.lend)
    jminus = (p.jstrand.values == "-")
    p["d_in"] = np.where(jminus, hi_in, lo_in)
    p["a_in"] = np.where(jminus, lo_in, hi_in)
    p["spanning"] = (p.jstart < p.lstart) & (p.jend > p.lend)
    p["cls"] = np.where(
        p.spanning, "span",
        np.where(p.d_in & p.a_in, "both",
                 np.where(p.d_in, "donor", np.where(p.a_in, "acceptor", "outside"))))
    inside = p[p.cls.isin(["donor", "acceptor", "both"])]
    span = p[p.cls == "span"]
    keep = pd.concat([
        inside.sort_values("sc", ascending=False).groupby("locus_uid").head(CAP_INSIDE),
        span.sort_values("sc", ascending=False).groupby("locus_uid").head(CAP_SPAN),
    ], ignore_index=True)
    return keep.sort_values(["locus_uid", "sc"], ascending=[True, False])


def anchor_index(intervals: pd.DataFrame) -> dict:
    """Build a per-chrom (starts, ends) lookup from a mappability interval set.

    `intervals` needs chrom/start/end and MUST already be merged and
    non-overlapping -- `_anchor_flags` relies on that to test containment with a
    single searchsorted rather than a scan.
    """
    out = {}
    for c, g in intervals.groupby("chrom", observed=True):
        g = g.sort_values("start")
        out[str(c)] = (np.ascontiguousarray(g.start.to_numpy(np.int64)),
                       np.ascontiguousarray(g.end.to_numpy(np.int64)))
    return out


def _anchor_flags(chrom, jstart, jend, index, term=ANCHOR_WIN):
    """1 where EITHER outer anchor window is wholly inside one easy interval.

    The windows are the exonic flanks -- [jstart-term, jstart) upstream of the
    donor and [jend, jend+term) downstream of the acceptor -- because that is
    where a split read's anchoring segments actually sit.  Junction coordinates
    are 0-based half-open intron bounds (verified: 40/40 top-support junctions
    read GT..AG at zero offset), so no +/-1 correction belongs here.

    A window on a contig absent from `index` yields 0 and is counted separately;
    the caller decides whether that is acceptable.
    """
    n = len(jstart)
    flags = np.zeros(n, np.uint8)
    nodata = 0
    for i in range(n):
        t = index.get(str(chrom[i]))
        if t is None:
            nodata += 1
            continue
        st, en = t
        hit = 0
        for ws, we in ((jstart[i] - term, jstart[i]), (jend[i], jend[i] + term)):
            k = int(np.searchsorted(en, ws, "right"))
            if k < len(st) and st[k] <= ws and en[k] >= we:
                hit = 1
                break
        flags[i] = hit
    return flags, nodata


def encode(keep: pd.DataFrame, anchor_idx: dict | None = None) -> tuple[dict, dict]:
    """(per-locus packed words, stats).  Offsets are 5'-relative and signed."""
    five = np.where(keep.lstrand == "+", keep.lstart, keep.lend)
    sgn = np.where(keep.lstrand == "+", 1, -1)
    # Resolve each junction END to a genomic position FIRST, using the
    # junction's own strand (donor = the junction's 5' end = the higher
    # coordinate when the junction is on the minus strand). Only then convert
    # to locus-5'-relative signed offsets.
    #
    # This previously flipped on the LOCUS strand, which is a different
    # question: it made "donor offset" mean "lower coordinate" on plus-strand
    # loci and "higher coordinate" on minus-strand loci, neither of which is
    # the donor for an antisense junction. See select() for the 27% figure.
    jminus = (keep.jstrand.values == "-")
    donor_pos = np.where(jminus, keep.jend.values, keep.jstart.values)
    acc_pos = np.where(jminus, keep.jstart.values, keep.jend.values)
    od = (donor_pos - five) * sgn
    oa = (acc_pos - five) * sgn
    sat = int(((od < _DONOR_MIN) | (od > _DONOR_MAX)
               | (oa < _ACC_MIN) | (oa > _ACC_MAX)).sum())
    od = np.clip(od, _DONOR_MIN, _DONOR_MAX)
    oa = np.clip(oa, _ACC_MIN, _ACC_MAX)
    anti = (keep.lstrand.values != keep.jstrand.values).astype(np.uint8)
    if anchor_idx:
        amap, a_nodata = _anchor_flags(keep.chrom.values, keep.jstart.values,
                                       keep.jend.values, anchor_idx)
    else:
        amap, a_nodata = np.zeros(len(keep), np.uint8), len(keep)
    words = pack(od, oa, anti, log_code(keep.sc.values), log_code(keep["cov"].values),
                 anchor_ok=amap)
    out = {}
    for uid, w in zip(keep.locus_uid.values, words):
        out.setdefault(uid, []).append(int(w))
    stats = {"junctions": int(len(keep)), "loci": len(out), "saturated": sat,
             "antisense": int(anti.sum()), "anchor_ok": int(amap.sum()),
             "anchor_nodata": int(a_nodata),
             "max_per_locus": max((len(v) for v in out.values()), default=0),
             "bytes_raw": int(len(keep) * 8)}
    return out, stats


def build(jx_path: str, loci: pd.DataFrame, canonical_only: bool = True,
          anchor_intervals: pd.DataFrame | None = None):
    """Full pipeline: load, filter, assign, cap, pack.  `loci` needs
    locus_uid/chrom/start/end/strand for one assembly.

    `anchor_intervals` is a MERGED mappability interval set (chrom/start/end);
    when given, the anchor-mappable bit is computed from it.  When omitted the
    bit is 0 everywhere and `stats['anchor_nodata']` equals the junction count,
    which is how a caller can tell "not measured" from "measured as unmappable".
    """
    jx = pd.read_parquet(jx_path)
    jx = jx[jx.sc >= SC_FLOOR]
    if canonical_only:
        jx = jx[jx.canonical]
    pairs = _assign(jx, loci)
    keep = select(pairs)
    aidx = anchor_index(anchor_intervals) if anchor_intervals is not None else None
    words, stats = encode(keep, aidx)
    stats["pairs_considered"] = int(len(pairs))
    return words, stats, keep
