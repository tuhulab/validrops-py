"""Real-dataset fixture registry for the benchmark harness (Package 2+).

Each entry describes a public raw count matrix (10x H5, 10x MTX, or AnnData
H5AD with counts) placed under the git-ignored ``tests/data/<name>/``
directory. Fixtures are never committed (see ``.gitignore``); the harness
reports ``available: false`` and skips the family when a fixture is absent,
exactly like the Package-1 PBMC4K gate. Downloading is a documented,
size-checked step in ``benchmarks/_fetch_fixtures.py`` (see also
``benchmarks/README.md``).

The package-3 sweep covers five datasets across three technologies (10x v3
scRNA-seq, 10x single-nucleus, 10x Visium v1 spatial) and two species (human,
mouse), plus a filtered/unfiltered condition pair on the same tissue:
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]

#: Loader names understood by :func:`load_fixture_adata`.
LOADER_10X_H5 = "10x_h5"
LOADER_H5AD = "h5ad"
LOADER_10X_MTX = "10x_mtx"

#: One entry per benchmark dataset, keyed by the ``--dataset`` CLI value.
FIXTURES = {
    "pbmc-parotid": {
        "path": "tests/data/pbmc-parotid/raw.h5",
        "loader": LOADER_10X_H5,
        "tissue": "parotid (salivary) gland",
        "species": "human",
        "chemistry": "10x 3' v3 feature-barcode (Cell Ranger 7.0.1)",
        "condition": "filtered (Cell Ranger cell call)",
        "license": "CC-BY-4.0 (Zenodo record 7470666; GSE188478, Chen et al. 2021)",
        "n_cells": 20_015,
        "n_genes": 36_601,
        "recommended_cap": 1_500,
        "fetch": "python3 benchmarks/_fetch_fixtures.py",
    },
    "pbmc-parotid-raw": {
        "path": "tests/data/pbmc-parotid/raw_subset.h5ad",
        "loader": LOADER_H5AD,
        "tissue": "parotid (salivary) gland",
        "species": "human",
        "chemistry": "10x 3' v3 feature-barcode (Cell Ranger 7.0.1)",
        "condition": "unfiltered; barcodes 3001-4500 by descending UMI (crosses the rank elbow: real cells + ambient droplets)",
        "license": "CC-BY-4.0 (Zenodo record 7470666; GSE188478, Chen et al. 2021)",
        "n_cells": 1_500,
        "n_genes": 36_601,
        "recommended_cap": 1_500,
        "fetch": "python3 benchmarks/_fetch_fixtures.py --extra",
    },
    "mouse-brain-5k": {
        "path": "tests/data/mouse-brain-5k/adata.h5ad",
        "loader": LOADER_H5AD,
        "tissue": "adult mouse brain",
        "species": "mouse",
        "chemistry": "10x 3' single-nucleus (5k Adult Mouse Brain Nuclei, Cell Ranger 7.0.0)",
        "condition": "filtered; counts in raw.X (X is normalized)",
        "license": "Zenodo record 7993711 (scDGD training data; permissive reuse)",
        "n_cells": 7_377,
        "n_genes": 32_285,
        "recommended_cap": 1_500,
        "fetch": "python3 benchmarks/_fetch_fixtures.py --extra",
    },
    "carotid-visium": {
        "path": "tests/data/carotid-visium/raw.h5",
        "loader": LOADER_10X_H5,
        "tissue": "carotid artery wall (aorta)",
        "species": "human",
        "chemistry": "10x Visium v1 spatial (Space Ranger filtered_feature_bc_matrix)",
        "condition": "filtered; 1,453 tissue spots",
        "license": "Zenodo record 22024050 (no explicit license in record)",
        "n_cells": 1_453,
        "n_genes": 36_601,
        "recommended_cap": 1_453,
        "fetch": "python3 benchmarks/_fetch_fixtures.py --extra",
    },
    "mouse-olfactory-bulb": {
        "path": "tests/data/mouse-olfactory-bulb/mtx",
        "loader": LOADER_10X_MTX,
        "tissue": "adult mouse olfactory bulb",
        "species": "mouse",
        "chemistry": "10x Visium v1 spatial (V1, 10x demo dataset)",
        "condition": "tissue spots with count matrix exported from the record's list.rds (one-time conversion, no R at test time)",
        "license": "Zenodo record 15024747 (10x demo data redistribution)",
        "n_cells": 1_185,
        "n_genes": 32_285,
        "recommended_cap": 1_185,
        "fetch": "python3 benchmarks/_fetch_fixtures.py --extra",
    },
}

#: Default dataset for families that take ``--dataset`` (tier-m-real, tier-h2h-qc).
DEFAULT_DATASET = "pbmc-parotid"


def known_datasets() -> tuple[str, ...]:
    return tuple(sorted(FIXTURES))


def fixture_expected_path(dataset: str) -> str:
    """Repo-rooted fixture path for a dataset (never depends on the CWD)."""
    if dataset not in FIXTURES:
        raise ValueError(f"unknown dataset {dataset!r}; choose from {known_datasets()}")
    return str(REPO_ROOT / FIXTURES[dataset]["path"])


def resolve_fixture(dataset: str, override: str = "auto", path: str | None = None) -> dict:
    """Detect whether a dataset fixture is present; rooted at the repo.

    Mirrors the Package-1 ``resolve_pbmc_state`` contract: ``override='yes'``
    raises if the fixture is missing, ``'no'`` forces unavailability, and
    ``auto`` reports on-disk presence.
    """
    requested = path or fixture_expected_path(dataset)
    exists = Path(requested).exists()  # files and mtx-dir fixtures alike
    if override == "yes" and not exists:
        raise FileNotFoundError(f"{dataset} fixture requested (override yes) but missing: {requested}")
    available = exists and override != "no"
    return {
        "available": available,
        "path": requested,
        "override": override,
        "dataset": dataset,
        "meta": dict(FIXTURES[dataset]),
    }


def load_fixture_adata(path: str, n_cells_cap: int):
    """Load a fixture as an AnnData (cells x genes, counts), bounded to ``cap``.

    Dispatch on the registry loader; every loader applies the shared
    conventions: unique var names, positional ``cell_1 .. cell_N`` obs names
    (R positional convention), and the cap applied by slicing after load.
    """
    import warnings

    import scanpy as sc

    spec = FIXTURES[_dataset_for_path(path)]
    loader = spec["loader"]
    if loader == LOADER_10X_H5:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Variable names are not unique")
            adata = sc.read_10x_h5(path)
        adata.var_names_make_unique()
    elif loader == LOADER_H5AD:
        adata = sc.read_h5ad(path)
        counts_src = adata.raw if (adata.raw is not None and adata.raw.X is not None) else adata
        # Some upstream h5ads store their counts in raw.X but lose the raw slot's
        # var names (degenerate '0'..'N' index); the main adata keeps the real
        # gene symbols and both matrices are column-aligned (same shape).
        names = counts_src.var_names
        if len(names) and all(str(n).isdigit() for n in names[:64]):
            names = adata.var_names
        cap = min(int(n_cells_cap), adata.n_obs)
        matrix = counts_src.X[:cap]  # csr slice, backed-safe
        if hasattr(matrix, "toarray"):
            matrix = matrix
        import anndata as ad

        adata = ad.AnnData(X=matrix, var=pd.DataFrame(index=pd.Index(names, dtype=str)))
        adata.var_names_make_unique()
    elif loader == LOADER_10X_MTX:
        import gzip

        import numpy as np
        import scipy.io
        import scipy.sparse as sp

        mtx_dir = Path(path)

        def _read_lines(name: str) -> list[str]:
            p = mtx_dir / name
            if not p.exists():
                p = Path(str(p) + ".gz")
            opener = gzip.open if str(p).endswith(".gz") else open
            with opener(p, "rt") as fh:
                return [line.rstrip("\n") for line in fh]

        matrix = scipy.io.mmread(mtx_dir / "matrix.mtx.gz")
        features = _read_lines("features.tsv")
        barcodes = _read_lines("barcodes.tsv")
        import anndata as ad

        adata = ad.AnnData(
            X=sp.csr_matrix(matrix.T, dtype=np.float32),
            var=pd.DataFrame(index=pd.Index(features, dtype=str)),
            obs=pd.DataFrame(index=pd.Index(barcodes, dtype=str)),
        )
        adata.var_names_make_unique()
    else:  # pragma: no cover - registry is internal
        raise ValueError(f"unknown loader {loader!r} for {path}")
    adata.obs_names = [f"cell_{i + 1}" for i in range(adata.n_obs)]
    cap = min(int(n_cells_cap), adata.n_obs)
    return adata[:cap].copy()


def _dataset_for_path(path: str) -> str:
    for name, spec in FIXTURES.items():
        if str(REPO_ROOT / spec["path"]) == str(path):
            return name
    raise ValueError(f"fixture path {path!r} is not registered in FIXTURES")
