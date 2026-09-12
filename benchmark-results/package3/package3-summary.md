# Package 3 — preprint evidence sweep

**Scope:** the cross-dataset sweep and efficiency experiments that form the
evidence base of the NAR GAB manuscript skeleton (`manuscript/manuscript.md`).
Repo revision recorded per artifact; working tree dirty by design.

## 1. Dataset matrix

Five datasets registered (`benchmarks/_fixtures.py`), three technologies, two
species, plus a filtered/unfiltered condition pair. All fixtures are
git-ignored; downloads are size-checked `benchmarks/_fetch_fixtures.py --extra`.

| Dataset | Species | Tissue | Technology / condition | Cells or spots | Source (license) |
|---|---|---|---|---|---|
| pbmc-parotid | human | parotid gland | 10x 3' v3 scRNA-seq, filtered | 20,015 | Zenodo 7470666 (CC-BY-4.0) |
| pbmc-parotid-raw | human | parotid gland | 10x 3' v3, unfiltered; barcodes 3001-4500 (elbow-crossing) | 1,500 (of 2,205,653) | Zenodo 7470666 (CC-BY-4.0) |
| mouse-brain-5k | mouse | adult brain | 10x 3' single-nucleus | 7,377 | Zenodo 7993711 |
| carotid-visium | human | carotid artery wall | 10x Visium v1 spatial | 1,453 spots | Zenodo 22024050 |
| mouse-olfactory-bulb | mouse | olfactory bulb | 10x Visium v1 spatial | 1,185 spots | Zenodo 15024747 |

**Sourcing note (species axis, per the goal's substitution clause):** the
planned third species could not be sourced. Attempted, after the 10x CDN WAF
block: Zenodo searches (`zebrafish 10x h5`, `"Danio rerio" 10x`, `rat 10x
chromium`, `rat snRNA-seq 10x h5ad`, `rat brain nuclei h5ad`, `worm C. elegans
10x`), figshare API search (returned unrelated records), and the full Zenodo
record space for raw-count files. Every reachable rat/zebrafish record was
either processed (no raw counts), >800 MB with an unclear counts layout, or
license-obscure. The sweep therefore covers **2 species (human, mouse)**, not
3, and this is stated in the manuscript's limitations. A third species is one
verified registry entry away once a source is found.

## 2. Sweep results (3 reps, profile S, seed 0; all CV <= 1.8%, all deterministic, zero warnings)

### QC head-to-head vs Seurat-style defaults (`benchmark-results/package3/h2h/`)

| Dataset | kept validrops | kept Seurat-style | kept-set Jaccard | validrops mito cap | QC validrops (s) | QC baseline (s) |
|---|---:|---:|---:|---:|---:|---:|
| pbmc-parotid | 51.1% | 98.6% | 0.518 | 10.7% (gmm-uik) | 24.2 | 13.4 |
| pbmc-parotid-raw | 1.6% | 0.8% | 0.500 | 10.3% | 18.7 | 9.5 |
| mouse-brain-5k | 84.3% | 94.6% | 0.891 | 8.8% | 23.6 | 12.7 |
| carotid-visium | 35.2% | 29.3% | 0.687 | 13.9% | 20.8 | 10.4 |
| mouse-olfactory-bulb | 94.1% | 77.0% | 0.761 | 28.1% | 23.4 | 11.7 |

The adaptive cap spans 8.8%-28.1%: stricter than the fixed 20% on parotid
scRNA, more permissive on olfactory-bulb Visium. Kept-set differences equal
the number of barcodes in the 10.7-28.1% mito band (verifier records
`kept_set_difference` per artifact).

### Stage-3 transfer (`benchmark-results/package3/real/`)

| Dataset | median (s) | clusters | n_obs | species | deterministic |
|---|---:|---:|---:|---|---|
| pbmc-parotid | 31.2 | 6 | 1,500 | human | yes |
| pbmc-parotid-raw | 193.4 | 2 | 1,500 | human | yes |
| mouse-brain-5k | 37.0 | 19 | 1,500 | mouse | yes |
| carotid-visium | 23.4 | 5 | 1,453 | human | yes |
| mouse-olfactory-bulb | 28.4 | 7 | 1,185 | mouse | yes |

Note (honest): Stage-3 on the unfiltered elbow-crossing slice is 6.2x slower
and yields 2 clusters — running Stage-3 on ambient-dominated barcodes is not a
meaningful workload; the meaningful raw-dataset measurement is the QC family
above. The number documents the cost of skipping Stage-1/2 filtering.

## 3. Efficiency experiments

### A — deterministic sparse SVD (implemented, exact)

Stage-3 `_embed` threaded `random_state` but never passed it to `svds`;
ARPACK's v0 came from OS entropy, so embeddings downstream were
run-to-run nondeterministic on real matrices (the Package-1 harness recorded
this as a pre-existing warning). Fix: seed the call, as Stage-4's copy already
did. Evidence: `test_embed_is_deterministic_for_fixed_seed` (bit-identical
per seed); every package-3 artifact reports `deterministic_warmup_timed: true`
**with zero warnings** (previously 2 svds-entropy warnings per PBMC4K sample).
Rerun of tier-m-real smoke on parotid after the fix: warnings 2 -> 0.

### B — Stage-4 ridge threads (measured, no change)

`benchmarks/scripts/measure_stage4_threads.py`, pbmc-parotid cap 800, 1 rep
per config, thread vars set before imports (artifacts:
`benchmark-results/package3/stage4-threads/`):

| threads | median (s) | dead callouts | flag |
|---|---:|---|---|
| 1 | 1,044.0 | 1 | Failed |
| 2 | 1,052.8 | 2 | Failed |
| 4 | 1,074.7 | 2 | Failed |

Flat (+3.0% from 1 to 4 threads): no threading benefit; the consensus loop is
`joblib.Parallel(n_jobs=1)` by default, so no joblib oversubscription exists
to fix. No runtime dependency added. The measurement surfaced a real crash
(single-class ridge calibration -> numpy zero-size reduction; see below),
fixed and re-run: the sweep now completes with the R-consistent `Failed`
soft-fail. The +-1 dead-callout difference across thread limits is
[INFERENCE] sklearn `GaussianMixture` EM float-order sensitivity in
`quality_filter` — SEEDED but not bitwise-stable across BLAS thread counts, a
property of sklearn, not validrops code.

### C — hotspot optimization (implemented, exact, 4.15x family-level)

Profile of tier-h2h-qc (parotid 1,500): `quality_metrics` = 9.25 s of ~11 s
(88%); of that, 8.4 s in 6 `numpy.isin` calls on object-dtype gene names
against fixed-width `<U` set arrays (protein_coding alone 4.3 s; numpy string
padding semantics also DIVERGE from object membership — the fixed-width path
is not merely slow but wrong-typed for this port). Replaced with pandas
`isin(set)` hashing (python-object equality, bit-identical to the R port's
`%in%` and to the previous object path).

Gate evidence (artifacts `benchmark-results/package3/c-efficiency/`):
- snapshotted all five metric columns before/after on identical input —
  `np.array_equal(equal_nan=True)` — exact equality;
- `tier-h2h-qc current` before/after (2 reps, same seed): **23.15 s -> 5.58 s
  median (4.15x)**, kept-set Jaccard 0.5179 both, deterministic, zero
  warnings;
- R-reference quality-metrics tests still pass.

Above the 1.25x merge bar; no clustering/neighbor/SVD semantics touched.

### A robustness fix surfaced by the sweep (not an efficiency change)

Stage-4 `_one_run` crashed when an iteration's QC-passing subset was
single-class: sklearn's `roc_curve` returns empty arrays and
`roc_best_threshold` hit a numpy zero-size reduction (reproduced on
pbmc-parotid under BLAS thread limit 2). Fixed: degenerate response
soft-fails with `"Failed"` (R's pROC errors on single-class y), and the
helper raises a clear `ValueError`. Regression test:
`test_roc_best_threshold_rejects_single_class_input`.

## 4. Unreported changes keep port fidelity

- `_annotation.gene_sets` selection rewritten (set-hash membership) — outputs
  bit-identical, gene-set sizes unchanged (13/102/18,935 human parotid).
- `quality_metrics` selector + var masks rewritten the same way — five metric
  columns bit-identical.
- Fixture loaders added (h5ad with raw.X counts, mtx dir); mouse-brain-5k
  detection previously mis-fired as `human` because the h5ad's raw slot lost
  its var names — loader now takes names from the main adata (alignment
  verified: same shape). Not a port bug; a fixture-format quirk documented in
  `benchmarks/_fixtures.py`.

## 5. Reproducibility

- 15 sweep artifacts + 2 c-efficiency pairs + 3 stage4-threads artifacts —
  all pass `python benchmarks/scripts/validate_artifacts.py` (round-trip,
  schema, required keys).
- Figures regenerated deterministically: `manuscript/render_figures.py`.
- Fixture downloads documented + size-checked; licenses recorded per dataset.
