# Package 2 — Real-data generalization + head-to-head QC benchmark

**Scope:** extends the Package-1 harness (which measured only the PBMC4K
reference and deterministic synthetic data) toward the preprint's dataset
matrix: a registry-driven real-data tier (`tier-m-real`) and a head-to-head QC
family (`tier-h2h-qc`) that pits validrops' Stage 1+2 against the classic
Seurat-default QC subset on identical inputs. First non-PBMC fixture: human
parotid gland, 10x v3 (Cell Ranger 7.0.1), CC-BY-4.0 (Zenodo 7470666 /
GSE188478, Chen et al. 2021).

**Repo revision:** `HEAD` at run time (recorded per artifact); working tree
dirty by design (Package-2 additions + untracked fixture download).
**Hardware:** this machine (Linux x86_64, AMD EPYC-Genoa).
**Measurement:** same protocol as Package 1 — fresh spawned child per
repetition, `time.perf_counter`, `resource.getrusage` high-water inside the
child, median/min/MAD/CV, warm-up == timed fingerprint gate. The timed window
includes deterministic fixture loading (identical for baseline and current), so
the pipeline delta is the difference.

## 1. Results

### tier-m-real — Stage-3 expression_metrics on a bounded parotid subset

| Workload | Cells | Median (s) | CV | Correctness |
|---|---:|---:|---:|---|
| `tier-m-real` / current / S | 1,500 (of 20,015) | 30.96 | 0.73% | pass |

Verifier metrics: all 1,500 cells clustered, 6 clusters, stats columns match
`STAT_COLUMNS`, species detected `human`. The family proves the Stage-3
workload generalizes beyond PBMC4K on a different tissue + chemistry without
harness changes — adding a dataset for the preprint is one registry entry in
`benchmarks/_fixtures.py`.

### tier-h2h-qc — validrops vs Seurat-default QC on the identical 1,500-cell subset

| Implementation | Median (s) | CV | Correctness |
|---|---:|---:|---|
| baseline (Seurat-style defaults) | 14.09 | 3.0% | pass |
| current (validrops Stage 1+2) | 24.20 | 1.2% | pass |
| compare (current vs baseline) | speedup 0.58x | — | compatible |

Agreement on 1,500 cells: validrops keeps **766 (51.1%)**; Seurat-style keeps
**1,479 (98.6%)**; kept-set **Jaccard 0.518**; difference 713 cells. validrops'
adaptive mitochondrial threshold landed at **10.7%** (gmm-uik method) vs the
Seurat baseline's fixed 20% cap — the entire disagreement is the 713 cells with
mitochondrial fraction between 10.7% and 20%, i.e. the tools disagree exactly
on what validrops is designed to do differently (data-adaptive mito cutoff
instead of a universal cap). Both gates pass and warm-up == timed
(deterministic callouts, seed 0).

### R cross-check of the baseline (evidence, not a CI gate)

The Python "Seurat-style" baseline is a re-implementation of the published
defaults (`n_features >= 200`, `n_counts >= 500`, `mito <= 0.20`,
`min_cells = 3`, gene drop before any metric, exactly like
`CreateSeuratObject`). It is validated against real Seurat 5.5.1 on the same
full matrix (20,015 cells) exported to MTX
(`benchmarks/R/h2h_seurat_crosscheck.R`):

| Metric | R Seurat 5.5.1 | Python baseline |
|---|---|---|
| n_kept | 19,811 | 19,811 |
| n_fail_mito | 203 | 203 |
| n_fail_features | 18 | 18 |
| n_fail_counts | 0 | 0 |

Exact equality. The benchmark family itself needs no R; R appears only in this
optional cross-check.

## 2. What changed

- `benchmarks/_fixtures.py` — dataset registry (path, loader, tissue/species/
  chemistry/license metadata, recommended cap) + repo-rooted `resolve_fixture`
  + 10x H5 loader mirroring the PBMC4K conventions.
- `benchmarks/_measure.py` — two new families wired into the existing
  dispatch/verifier/skip machinery: `tier-m-real` (generic Stage-3 on a
  registry fixture) and `tier-h2h-qc` (validrops QC vs Seurat-style baseline,
  agreement metrics in the correctness record). Skip behavior for absent
  fixtures matches Package 1 (artifact with `available: false`, named expected
  path, correctness `skipped`).
- `benchmarks/run.py` — `--dataset` CLI option.
- `benchmarks/_fetch_fixtures.py` — documented fixture download (size-checked,
  re-runnable) + `--export-mtx` for the R cross-check.
- `benchmarks/R/h2h_seurat_crosscheck.R` — Seurat cross-check script.
- `tests/test_benchmarks.py` — registry-rooting test, absent-fixture skip tests
  for both families, and fixture-gated agreement/clustering tests (skip in CI,
  run when fixtures are present). `tests/data/**` stays git-ignored.

Unchanged: the PBMC4K family and all Package-1 artifacts, the port algorithms,
`_constants.py`, and everything outside the benchmark harness + tests.

## 3. Honest limitations

- validrops' Stage-1+2 is **slower** than the trivial threshold baseline
  (0.58x); the preprint framing is correctness/adaptivity (adaptive mito cap,
  residual bands), never runtime for this family.
- One dataset (parotid v3) is benchmarked so far; PBMC4K is absent from this
  worktree, so no cross-dataset agreement table yet. The registry makes the
  next datasets (mouse brain, Visium, other conditions) additive.
- The 1,500-cell cap bounds Louvain sweep cost; larger caps and full-matrix
  runs are deferred to the preprint benchmark sweep.
- `raw_all.h5` (unfiltered parotid) is staged for a future empty-drop/
  ambient comparison family; not consumed yet.
