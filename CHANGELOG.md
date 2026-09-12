# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][],
and this project adheres to [Semantic Versioning][].

[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html

## [Unreleased]

### Added

- Stage 1 `pp.rank_barcodes`: barcode ranking and knee-point detection separating cells from ambient droplets
- Stage 2a `tl.quality_metrics`: per-barcode mito/ribo/coding fractions and log UMI/feature counts
- Stage 2b `pp.quality_filter`: mito-fraction cap, feature↔UMI residual band, coding-fraction band
- Stage 3a `tl.expression_metrics`: deviance HVG selection, SVD embedding, SNN + Louvain clustering, per-cluster marker statistics
- Stage 3b `pp.expression_filter`: expression-based cluster outlier filtering
- Stage 4 `pp.label_dead`: heuristic dead-cell score and consensus ridge-training loop (off by default)
- `validrops()` end-to-end pipeline orchestrator annotating an `AnnData` in place
- `validrops.pl`: QC diagnostic plots (barcode rank, mito threshold, UMI vs features, coding fraction, dead score)
- Validation suite against R-generated fixtures on PBMC 4K, with a bounded end-to-end barcode concordance of 0.9193
- Package scaffolding: CI (test/build/release), Sphinx docs on readthedocs, pre-commit, codecov
- Benchmark harness (Package 2): `tier-m-real` registry-driven real-data Stage-3 family and `tier-h2h-qc` head-to-head QC family (validrops Stage 1+2 vs a Seurat-default baseline cross-validated against Seurat 5.5.1), with the first non-PBMC fixture — human parotid gland, 10x v3 (CC-BY-4.0) — documented in `benchmarks/README.md`
- Benchmark harness (Package 3): four more registry datasets (mouse brain nuclei, carotid Visium, mouse olfactory bulb Visium, unfiltered parotid), a cross-dataset sweep runner, artifact validator and a NAR GAB manuscript skeleton with artifact-regenerated figures (`manuscript/`)
- `tl.quality_metrics` hotspot optimization: gene-set membership now uses python-set hashing instead of `np.isin` on object-dtype arrays (`_annotation.gene_sets`, metric selectors, `adata.var` masks) — quality metrics 9.25 s -> 0.78 s on 1,500-cell parotid (11.8x), full QC benchmark family 23.15 s -> 5.58 s (4.15x), with bit-identical output (verified column-wise with `equal_nan=True` and kept-set Jaccard 0.5179 before/after)

### Fixed

- `test_pbmc_detection_rooted_at_repo` asserted the git-ignored PBMC4K fixture is present, which could never hold on CI; availability now mirrors on-disk presence (fixture-gated tests skip when absent)
- Sphinx `-W` build failed on `docs/internal/**` (untracked working notes never added to a toctree, with xrefs into excluded `docs/superpowers/**`); internal notes are now excluded from the published build like `docs/superpowers/**`
- Stage-3 sparse SVD (`expression_metrics._embed`) drew its ARPACK starting vector from OS entropy because the seeded `random_state` parameter was never passed to `svds`; embeddings (and clustering downstream) were run-to-run nondeterministic on real data, flagged by the benchmark harness as a pre-existing warning. The call is now seeded; output is bit-identical for a fixed seed and every real-data tier artifact reports `deterministic_warmup_timed: true` with zero warnings
- Stage-4 ridge consensus loop crashed with a numpy zero-size reduction when an iteration's QC-passing subset had a single class (surfaced on pbmc-parotid under a different BLAS thread order); the iteration now soft-fails with the R-consistent `"Failed"` flag, and `roc_best_threshold` raises a clear `ValueError` on single-class input
