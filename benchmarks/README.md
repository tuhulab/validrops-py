# valiDrops Performance Package 1 — benchmark harness

Deterministic, spawn-isolated micro-benchmarks for the exact-CPU optimization
work. Everything here lives **outside the installed package** and writes
artifacts only under the gitignored `benchmark-results/`.

## Quick start

```bash
# Schema/smoke check (one repetition, allowed only for --profile smoke)
rtk uv run python -m benchmarks.run --benchmark cluster-means --implementation baseline \
    --profile smoke --output-dir benchmark-results/package1/smoke --repetitions 1

# Tier-S baseline (this is the "before" artifact)
rtk uv run python -m benchmarks.run --benchmark cluster-means --implementation baseline \
    --profile S --output-dir benchmark-results/package1/baseline
rtk uv run python -m benchmarks.run --benchmark graph-conversion --implementation baseline \
    --profile S --output-dir benchmark-results/package1/baseline
rtk uv run python -m benchmarks.run --benchmark tier-s --implementation baseline \
    --profile S --output-dir benchmark-results/package1/baseline --repetitions 3
```

Each command writes `<benchmark>-<implementation>.json` (schema version 1) plus a
Markdown summary to `--output-dir`. Repetitions default to 7 for the two
kernels and 3 for the pipeline tiers; `--repetitions 1` is only allowed for
`smoke` (schema checks).

## Families

| Family | What is timed | Baseline implementation | Current implementation |
|---|---|---|---|
| `cluster-means` | `expm1` column means of a sparse matrix | dense float64 oracle (cast to float64 before `expm1`, float64 mean) | `_sparse_expm1_col_mean` (sparse reduction, no `todense()`) |
| `graph-conversion` | deep-clustering resolution search (`louvain` per resolution) | `louvain(adjacency, res)` — igraph built per call | `prepare_louvain_graph` once + `louvain_prepared` per call |
| `tier-s` | full Stage-3 `expression_metrics` on deterministic synthetic data | current production code | same entry point, run after the optimization |
| `tier-m-pbmc` | full Stage-3 on a bounded real PBMC4K subset | current production code | same, with a configured run |
| `tier-m-real` | full Stage-3 on a bounded subset of any registered dataset (`--dataset`) | current production code | same, with a configured run |
| `tier-h2h-qc` | full QC (Stage 1 + 2) vs a Seurat-default baseline on an identical bounded subset | Seurat-style default thresholds (`n_features >= 200`, `n_counts >= 500`, `mito <= 0.20`, `min_cells = 3`), re-implemented deterministically in Python | `pp.rank_barcodes` + `tl.quality_metrics` + `pp.quality_filter` |

`tier-m-pbmc` detects the fixture at `tests/data/pbmc4k/raw.h5` (rooted at the
repository, never the CWD). With `--pbmc-override no`, or when the fixture is
absent, the artifact reports `available: false`, the expected path, and
correctness `skipped` — it never fabricates a run.

`tier-m-real` and `tier-h2h-qc` resolve fixtures through the registry in
`benchmarks/_fixtures.py` and behave the same way when a dataset is absent:
artifacts record `fixture.available: false`, the repo-rooted expected path, and
correctness `skipped`.

Honest naming: `tier-s`, `tier-m-pbmc` and `tier-m-real` are **Stage-3
workloads** (barcode QC stages excluded); thread experiments on them say
nothing about Stage-4 ridge oversubscription. `tier-h2h-qc` is a QC-stage
agreement benchmark: the timed region includes fixture loading (identical for
baseline and current), so runtimes are directly comparable and the pipeline
delta is the difference.

## Real-data fixtures (Package 2)

Fixtures are **never committed** (git-ignored `tests/data/`); downloading is a
documented, reproducible step under the `benchmarks/_fetch_fixtures.py` script
with size checks. Registry entries live in `benchmarks/_fixtures.py`; adding a
dataset for the preprint (different tissue/condition/technology) is one entry.

| Dataset | `--dataset` | Fixture path | Tissue / species | Chemistry | License / source |
|---|---|---|---|---|---|
| PBMC 4K | (tier-m-pbmc) | `tests/data/pbmc4k/raw.h5` | PBMC / human | 10x v2 | 10x Genomics demo; 17 MB |
| Parotid gland | `pbmc-parotid` | `tests/data/pbmc-parotid/raw.h5` | salivary gland / human | 10x v3 (Cell Ranger 7.0.1) | CC-BY-4.0, Zenodo 7470666 (GSE188478); 45 MB (plus 75 MB raw) |

```bash
python3 benchmarks/_fetch_fixtures.py             # downloads both parotid H5 files
python3 benchmarks/_fetch_fixtures.py --export-mtx  # exports tests/data/pbmc-parotid/mtx/ for the R cross-check
```

The parotid dataset is deliberately not PBMC: different tissue and a newer
chemistry than the PBMC4K reference, so the harness can show the pipeline
generalises beyond the reference fixture (run `tier-m-real --dataset
pbmc-parotid` and read `species_detected`, `n_clusters`). Its `raw_all.h5`
(unfiltered, with empty droplets) is retained for a future ambient/empty-drop
comparison family.

## Head-to-head vs Seurat (tier-h2h-qc)

`tier-h2h-qc` compares validrops' Stage 1+2 callouts against a deterministic
re-implementation of the classic **Seurat default QC subset** on the identical
bounded fixture. The baseline thresholds are part of the artifact parameters
(`seurat_mito_cap`, `seurat_min_features`, `seurat_min_counts`,
`seurat_min_cells`) and are cross-validated against **real Seurat** in
`benchmarks/R/h2h_seurat_crosscheck.R` (R only needed for the cross-check, never
for the benchmark itself):

```bash
Rscript benchmarks/R/h2h_seurat_crosscheck.R tests/data/pbmc-parotid/mtx
```

On the full parotid matrix (20,015 cells) the Python baseline reproduces Seurat
5.5.1 exactly: 19,811 kept (203 mito-breaches, 18 low-feature, 0 low-count).
The agreement metrics (kept-set Jaccard, per-tool kept fractions, validrops'
adaptive mito threshold vs Seurat's fixed cap) are the preprint-facing numbers;
`compare_artifacts` on the baseline/current pair gives the runtime comparison.

## Package 3 — preprint evidence sweep

The preprint sweep (`benchmark-results/package3/`) runs both real-data families
on every registered dataset:

```bash
python3 benchmarks/_fetch_fixtures.py --extra   # mouse brain nuclei, carotid Visium, olfactory bulb, parotid raw subset
uv run --group test --group dev python benchmarks/scripts/run_sweep.py --reps 3 --profile S
uv run --group test --group dev python benchmarks/scripts/validate_artifacts.py
uv run --group test --group dev python manuscript/render_figures.py
```

Registered datasets (see `benchmarks/_fixtures.py` for loaders and metadata):

| `--dataset` | Species / tissue | Technology | Fixture |
|---|---|---|---|
| `pbmc-parotid` | human parotid | 10x 3' v3 scRNA-seq | `tests/data/pbmc-parotid/raw.h5` |
| `pbmc-parotid-raw` | human parotid | 10x 3' v3, unfiltered (ambient) | `tests/data/pbmc-parotid/raw_subset.h5ad` (barcodes 3001-4500) |
| `mouse-brain-5k` | mouse brain | 10x 3' single-nucleus | `tests/data/mouse-brain-5k/adata.h5ad` (counts in raw.X) |
| `carotid-visium` | human carotid artery | 10x Visium v1 spatial | `tests/data/carotid-visium/raw.h5` |
| `mouse-olfactory-bulb` | mouse olfactory bulb | 10x Visium v1 spatial | `tests/data/mouse-olfactory-bulb/mtx/` (gzipped MTX) |

The manuscript skeleton (`manuscript/manuscript.md`), the deterministic figure
renderer (`manuscript/render_figures.py`), artifact validation
(`benchmarks/scripts/validate_artifacts.py`) and the Stage-4 ridge thread
measurement (`benchmarks/scripts/measure_stage4_threads.py`) all read from
these committed artifacts only — no network, no live package state.

Package-3 Efficiency C (exact-output hotspot optimization): `tl.quality_metrics`
gene-set membership 9.25 s -> 0.78 s on 1,500 parotid cells (11.8x; bit-identical
output, before/after artifacts in `benchmark-results/package3/c-efficiency/`).
Efficiency A (deterministic Stage-3 SVD) and the Stage-4 thread measurement
(criterion B) are documented in `package3-summary.md`.

## Measurement protocol

- Every measured repetition runs in a fresh `multiprocessing("spawn")` child.
- Warm-up runs inside the child, untimed, before the measured run.
- Timing: `time.perf_counter()`; RSS: `resource.getrusage(RUSAGE_SELF)` high
  water read inside the child (bytes on macOS, KiB normalized to bytes on
  Linux), so native numpy/scipy/igraph allocations are included.
- Report median, minimum, median-absolute deviation, raw samples, and
  coefficient of variation. A CV above 10% is flagged inconclusive, never
  hidden.
- Inputs are regenerated deterministically from `--seed` in each child; a
  warm-up/timed fingerprint equality check doubles as the determinism gate.
- Artifacts record git revision/dirty state, OS/architecture, logical and
  physical CPU counts (when discoverable), available memory, Python/package
  versions, NumPy BLAS config, effective native-thread environment variables,
  input fingerprints and exact parameters, source and harness file hashes.
- Incompatible comparisons (different benchmark, profile, or parameters) are
  rejected by `compare_artifacts`.

## Thread-control experiments (Package 1)

Decision gate: do **not** add `threadpoolctl` or alter package runtime behavior
in Package 1. These experiments only **document** the effect of native thread
limits on the Stage-3 workload so a future Stage-4 plan can decide whether a
runtime dependency is justified.

First start a fresh interpreter with the native thread env vars set **before**
NumPy/SciPy import (artifact `environment.thread_environment` records what the
child inherited):

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  rtk uv run python -m benchmarks.run --benchmark tier-s --implementation current \
  --profile S --output-dir benchmark-results/package1/threads/1 --repetitions 3
```

Repeat with limits 2 and 4 (`threads/2`, `threads/4`). If Tier-S is too slow on
the machine, run `smoke` first and label `S` as unavailable rather than
extrapolating. Conclusions from these Stage-3 measurements must not be
extrapolated to Stage-4 ridge oversubscription.
