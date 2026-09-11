#!/usr/bin/env Rscript
# Cross-check the harness's Seurat-style QC baseline against real Seurat.
#
# The Python baseline in benchmarks/_measure.py (`_seurat_style_qc`) re-implements
# the classic Seurat QC defaults deterministically so the benchmark family runs
# without an R runtime. This script validates that baseline: it applies the SAME
# thresholds with Seurat 5.x on the SAME 10x H5 and prints the per-condition and
# kept-set counts. The numbers feed the package-2 summary; they are evidence that
# the Python baseline equals what a Seurat user would measure, not a CI gate.
#
# Usage:
#   Rscript benchmarks/R/h2h_seurat_crosscheck.R tests/data/pbmc-parotid/raw.h5
#
# Output (CSV to stdout):
#   metric,value
#   n_total,<cells>
#   n_kept_seurat,<cells passing all three thresholds>
#   n_fail_mito,<cells with percent.mt > mito_cap>
#   n_fail_features,<cells with nFeature_RNA < min_features>
#   n_fail_counts,<cells with nCount_RNA < min_counts>
#   mito_cap,0.2
#   min_features,200
#   min_counts,500
#   min_cells,3
#   seurat_version,<version>
#
# Threshold constants MUST match the harness parameters
# (seurat_mito_cap=0.20, seurat_min_features=200, seurat_min_counts=500,
#  seurat_min_cells=3 in `_family_params("tier-h2h-qc", ...)`).

MITO_CAP <- 0.20
MIN_FEATURES <- 200
MIN_COUNTS <- 500
MIN_CELLS <- 3

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("usage: Rscript h2h_seurat_crosscheck.R <10x_mtx_dir>")
}
path <- args[[1]]
if (!dir.exists(path)) {
    stop(sprintf("fixture directory not found: %s", path))
}

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
})

# Read the 10x MTX directory (matrix.mtx, features.tsv, barcodes.tsv). The
# matrix is exported from the same H5 the harness uses (see
# benchmarks/_fetch_fixtures.py --export-mtx), so R and the Python baseline
# see bit-identical counts. Column 2 of features.tsv carries gene symbols (the
# MT- pattern match needs them).
mat <- ReadMtx(
    mtx = file.path(path, "matrix.mtx"),
    cells = file.path(path, "barcodes.tsv"),
    features = file.path(path, "features.tsv"),
    feature.column = 2
)
obj <- CreateSeuratObject(counts = mat, min.cells = MIN_CELLS, min.features = 0)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^Mt-|^mt-")

n_feature <- obj$nFeature_RNA
n_count <- obj$nCount_RNA
pct_mt <- obj$percent.mt

keep <- (pct_mt <= 100 * MITO_CAP) & (n_feature >= MIN_FEATURES) & (n_count >= MIN_COUNTS)
rows <- c(
    n_total = as.integer(ncol(mat)),
    n_kept_seurat = as.integer(sum(keep)),
    n_fail_mito = as.integer(sum(pct_mt > 100 * MITO_CAP)),
    n_fail_features = as.integer(sum(n_feature < MIN_FEATURES)),
    n_fail_counts = as.integer(sum(n_count < MIN_COUNTS)),
    mito_cap = MITO_CAP,
    min_features = MIN_FEATURES,
    min_counts = MIN_COUNTS,
    min_cells = MIN_CELLS,
    seurat_version = as.character(packageVersion("Seurat"))
)
cat("metric,value\n")
for (name in names(rows)) {
    cat(sprintf("%s,%s\n", name, rows[[name]]))
}
