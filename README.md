# HERV catalog browser — static dashboard

A dependency-free browser over `herv_catalog.db` v1.1 (39,733 loci). No server process,
no build step, no JS framework. Two source files, one data directory.

## Run it

```bash
cd dash && python -m http.server 8000
# open http://localhost:8000
```

Any static host works (GitHub Pages, S3, `file://` will **not** — see Caveats).

## Layout

```
dash/
  index.html              search UI + keyspace dropdown + result list
  detail.js               detail panels + SVG locus graphic
  data/
    search_index.json.gz  1.4 MB  alias -> locus_uid, all 9 keyspaces
    shard_meta.json       bucket count + hash name
    loci/0..399.json.gz    13 MB  per-locus detail, djb2(uid) % 400
    gene_models.json.gz     ~3 MB  GENCODE V50 transcripts per locus window
```

Total ~17 MB. Only the index (1.4 MB) plus one 32 KB shard load for a given lookup.

## Search — two dropdown sections, because the keyspaces behave differently

**Identifier types** resolve to a single locus and autocomplete:

| type | distinct strings | ambiguous |
|---|---|---|
| `combined_id` | 46,957 | **11,475** |
| `versioned_id` | 79,463 | 412 |
| `telescope_id` | 14,968 | 0 |
| `hervd_id` | 20,228 | 1,307 |
| `ervmap_id` | 3,205 | 975 |
| `ervmap_alt_name` | 191 | 38 |

**Classifier types** resolve to a locus *list* — the whole menu is browsable with an
empty text box, since there are too few values to need typing:

`dfam_int_model` (118 values, median 182 loci each) · `dfam_accession` (252) · `repbase_name` (466)

### Ambiguity is surfaced, never silently resolved

11,475 `combined_id` strings map to **more than one locus**, because retired positional
strings get re-minted onto different loci across assignments. Example:
`ERV316A3_10p11.21c` → `HERVL000007` *and* `HERVL000008`.

The result list badges these (`2 loci`, `retired`) and clicking one opens a
disambiguation table rather than guessing. This is the same hazard documented in
`herv-catalog-compare`: a naive join on `combined_id` sends 4,778 loci to the wrong
place without erroring.

The detail header states which ID to cite (`versioned_id`) and which is the immutable
key (`locus_uid`).

## The locus graphic

SVG, lane-based, over the locus ±1 kb. Dashed vertical guides mark the locus boundaries.

| lane | source | coverage |
|---|---|---|
| segments (LTR blue / internal rust) | `locus_segment` | 14,790 loci (37.2%) — telescope-origin only |
| gEVE ORFs (strand arrows) | `hit_geve` | 3,917 loci (9.9%) |
| HERVarium domains | `hit_hervarium_domain` | 4,972 loci (12.5%) |
| gene models (thick exon / thin intron) | GENCODE V50, fetched | see below |

Loci without segments show the merged locus extent in grey with an explicit note —
that is a real property of the catalog (only telescope-origin loci carry RepeatMasker
segments), not a rendering failure.

## Gene models: why V50 specifically

`hit_gene` stores overlap metrics and gene *names* but **no gene coordinates**, so exon
structures had to be fetched. Choosing the track was not arbitrary — the first attempt
used GENCODE Basic V44 and left 8,080 loci with recorded overlaps but zero transcripts
in window, which looked like a catalog defect. It was not:

> Locus `HERVL000012` (chr10:37,930,137-37,934,250) has an `intronic` ZNF25 overlap with
> `overlap_bp` 4,112 — exactly the full locus span. ZNF25 spans 37,949,572-37,976,647 in
> GENCODE V20 **through V49** (15 kb away, no overlap) but 37,916,978-37,976,658 in
> **V50**, which contains the locus.

So `hit_gene` was computed against GENCODE V50 and the catalog is right. Pinning the
graphic to V50 keeps the drawing consistent with the overlap table it accompanies.

**Regenerating against a different version:** change `TRACK` in the fetch cell to
another `wgEncodeGencodeComp{V20..V50}` and rebuild `gene_models.json.gz`. Expect
disagreement with `hit_gene` if you do.

## Known gaps

- **T2T gene models are not rendered.** 16,608 of 57,927 `hit_gene` rows are `t2t`, and
  3,977 are hg38 RefSeq. The graphic is hg38/GENCODE only; the *table* shows all rows
  including T2T and RefSeq, so nothing is hidden — it is only undrawn.
- **The graphic is hg38-only.** T2T coordinates appear in the table with a UCSC hs1 link.
- Gene symbol drift between the catalog's annotation snapshot and V50 (e.g. `CCNY-AS1` →
  `LINC02634`) means a few table names won't match a drawn transcript label.

## Coordinates

Stored coordinates are **0-based half-open** (UCSC/BED). The UI displays 1-based
inclusive and labels it as such. This matters: the v1.1 correction fixed telescope-origin
starts that were 1-based in storage — anything exported before v1.1 has starts 1 bp high.

## Rebuilding

`build_dashboard.py` regenerates the whole `data/` directory from `herv_catalog.db`.
`index.html` and `detail.js` are static and are not written by it.

```bash
python build_dashboard.py --db herv_db/herv_catalog.db --out dash   # full, ~100 s
python build_dashboard.py --out dash --skip-genes                   # data only, ~40 s
python build_dashboard.py --out dash --validate-only                # ~3 s
```

Four stages, each skippable: search index → locus shards → gene models → validate.
Stage 3 reuses `gencode_v50_models.parquet` if present, so the 10-minute UCSC fetch
happens once; the per-tile cache in `genecache50/` makes an interrupted fetch resumable.

The validator is the point of the script, and it exits non-zero on failure:

| check | catches |
|---|---|
| strict JSON over all 402 files | the `NaN` trap below |
| all 280,499 identifier resolutions | missing loci, wrong shard hash |
| locus count vs index | dropped loci |
| `hit_gene` agreement, warns below 90% | wrong GENCODE version |

All three failure modes were verified by injecting them deliberately. The resolution
check is exhaustive rather than sampled because a strided sample missed a deleted
locus during that test.

Verified reproduction: rebuilding into a clean directory produces per-locus content
identical to the shipped bundle for all 39,733 loci, and byte-identical
`search_index.json.gz` and `gene_models.json.gz`.

## Repeat annotation and the two assemblies

`locus_repeat` holds RepeatMasker calls falling in each locus window (locus ±1 kb),
for **both** assemblies, parsed by one identical method (`parse_rmsk_out.py`) from:

| assembly | file | `rmsk_library` recorded |
|---|---|---|
| hg38 | `hg38.fa.out.gz` | `UCSC hg38.fa.out` — UCSC publishes no version file for this run |
| t2t (hs1) | `hs1.repeatMasker.out.gz` | `dc20181026-rb20181026` (from `hs1.repeatMasker.version.txt`) |

The two runs therefore use **different, and for hg38 unstated, library versions**. That is
why they are not merged, and why the cross-check below matters.

This is deliberately **separate** from `locus_segment`, which is hg38-only and carries
Telescope's own classification. Keeping them apart is what let us cross-check: 98.3 % of
`locus_segment` rows reproduce exactly (same interval, same name) in the hg38 rmsk parse,
and **100 %** have a same-name overlapping call — the 1.7 % differ only in fragment
boundaries (median 11 bp / 3 bp), as expected between library versions. Nothing is absent.
On screen the two lanes are labelled separately so provenance stays visible.

`.out` is **1-based inclusive**; the catalog is **0-based half-open**. The parser converts
(`start-1`) and its `--self-test` asserts the round-trip plus the abut/outside/spanning
boundary cases. Run it before trusting any re-parse.

A repeat between two nearby loci belongs to both windows, so `locus_repeat` is a
locus↔repeat *relation*: `COUNT(*)` is **not** a repeat count. `n_loci` records how many
windows each call falls in — divide by it, or `COUNT(DISTINCT chrom||start||end)`.

Every class is stored, including simple repeats, low complexity and satellites. The page
draws LTR/LINE/SINE/DNA by default and the **show all classes** checkbox adds the rest.

### Assembly fallback

The graphic prefers hg38 and falls back to t2t, which is what gives the 2,945 T2T-only
loci a graphic at all (previously: none). On the t2t path, lanes whose data exist only in
hg38 coordinates — gene models, gEVE ORFs, HERVarium domains, `locus_segment` — are
**omitted rather than mis-drawn**, and the panel heading says so. The validator fails if
any locus carries repeats for an assembly it has no coordinate on.

## Gene models: two assemblies, two annotation sources

`data/gene_models.json.gz` is **assembly-keyed**: `{"hg38": {uid: [tx]}, "t2t": {uid: [tx]}}`.
The two halves come from different annotation sets, and this is not a choice you can
flip without consequence:

| assembly | source | why |
|---|---|---|
| hg38 | GENCODE `wgEncodeGencodeCompV50` (UCSC API) | `hit_gene`'s hg38 rows are `source='gencode'` |
| t2t (hs1) | `hs1.ncbiRefSeq.gtf.gz`, stamped `ncbiRefSeq.2023-05-29` | `hit_gene`'s t2t rows are `source='refseq'`; UCSC publishes no GENCODE for hs1 |

Each is built in its **own** assembly's coordinates, so `detail.js` indexes the bundle
by assembly — it must never assume hg38. The build fails loudly if the bundle is the
old flat `{uid: [...]}` shape, because the page would then find no genes at all.

**Agreement against `hit_gene`** is the regression guard on annotation choice:

| assembly | loci checked | all names recovered |
|---|---|---|
| hg38 / GENCODE V50 | 20,740 | 97.6% |
| t2t / hs1 RefSeq | 14,757 | **100.0%** |

The t2t side is exact because it is literally the file `hit_gene` was computed from.

If you remember hg38 as **98.3%**: same data, stricter denominator. The old guard looped
over loci that *got* models, which silently excused 151 loci that `hit_gene` says have
genes but which the build produced nothing for. The guard now loops over expected loci,
so those 151 count against it (0.8% "none"). 97.6% is the honest figure; it is not a
regression.
The build warns if either drops below 90% — that is the signal you changed annotation
version, not that the join broke. (For hg38 the known 2.4% is unnamed GENCODE genes
and version drift; see "The GENCODE version trap" above.)

Coordinate conversion is where a gene lane goes silently wrong: GTF is **1-based
inclusive**, the UCSC API and this bundle are **0-based half-open**. `parse_refseq_gtf.py`
decrements every start, leaves every end, and its `--self-test` asserts that plus the
degenerate single-base feature (must come out 1 bp, not 0). The parser refuses to touch
real data unless the self-test passes.

Coverage: hg38 20,894 loci carry >=1 transcript; t2t 15,256. Of the 2,945 **t2t-only**
loci, 893 get a gene lane — the rest correctly show "no RefSeq transcript in this
window". That is a property of the annotation, not a gap in the build.

## Lane labels

The label gutter is `L=110` px. It was 62, which held ~9 monospace characters and
clipped **81%** of gene labels — `LOC124905662` rendered as `OC124905662`, which reads
as a different identifier rather than as a truncation. Labels longer than the gutter are
ellipsised by `laneLabel()` with the full string in a `<title>` tooltip, so nothing is
silently shortened.

## The NaN trap

**`json.dump` writes bare `NaN` for float NaN, and `NaN` is not legal JSON.** Python's
`json.load` accepts it; `JSON.parse` in every browser rejects it. A contaminated shard
fetches with HTTP 200 and then throws inside an `async` function, which surfaces as
"the click does nothing" with no console error unless you look for the rejected promise.
This shipped in the first build: 332 of 400 shards carried `NaN` in `band` (765) and
`origin` (4), from loci with no cytoband assignment.

Always serialize with non-finite values coerced to `null`:

```python
def clean(o):
    if isinstance(o, float): return None if (math.isnan(o) or math.isinf(o)) else o
    if isinstance(o, dict):  return {k: clean(v) for k, v in o.items()}
    if isinstance(o, list):  return [clean(v) for v in o]
    return o
json.dumps(clean(obj), separators=(",", ":"), allow_nan=False)   # allow_nan=False will raise
```

`allow_nan=False` turns the silent corruption into a build-time exception. Verify a
rebuild with `python jscheck.py` (syntax) and `node simall.mjs` (renders all 39,733 loci).

The page now fails loudly: shard parse errors, missing loci, and graphic exceptions all
print to the panel naming the shard and the cause, and `unhandledrejection`/`error` are
both trapped.

## Rebuilding

Shards are keyed by `djb2(locus_uid) % 400`, implemented identically in the build script
and in `index.html`. Do **not** substitute Python's `hash()` — it is salted per process,
so the page would request the wrong shard.
