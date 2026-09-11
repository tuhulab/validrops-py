"""Real-dataset fixture registry for the benchmark harness (Package 2+).

Each entry describes a public 10x-style raw count matrix that can be placed
under the git-ignored ``tests/data/<name>/`` directory. Fixtures are never
committed (see ``.gitignore``); the harness reports ``available: false`` and
skips the family when a fixture is absent, exactly like the Package-1 PBMC4K
gate. Downloading is a one-line documented command (see
``benchmarks/_fetch_fixtures.py`` and ``benchmarks/README.md``).
"""

from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

#: Loader names understood by :func:`load_fixture_adata`.
LOADER_10X_H5 = "10x_h5"

#: One entry per benchmark dataset, keyed by the ``--dataset`` CLI value.
FIXTURES = {
    "pbmc-parotid": {
        "path": "tests/data/pbmc-parotid/raw.h5",
        "loader": LOADER_10X_H5,
        "tissue": "parotid (salivary) gland",
        "species": "human",
        "chemistry": "10x v3 feature-barcode (Cell Ranger 7.0.1)",
        "license": "CC-BY-4.0 (Zenodo record 7470666; GSE188478, Chen et al. 2021)",
        "n_cells": 20_015,
        "n_genes": 36_601,
        "recommended_cap": 2_000,
        "fetch": "python3 benchmarks/_fetch_fixtures.py",
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
    exists = os.path.isfile(requested)
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
    """Load a 10x H5 fixture as an AnnData, bounded to the first ``cap`` cells.

    Matches the PBMC4K loader conventions: unique var names, positional
    ``cell_1 .. cell_N`` obs names (R positional convention), float32 CSR kept
    as-is, capped copy after loading.
    """
    import warnings

    import scanpy as sc

    with warnings.catch_warnings():
        # read_10x_h5 warns about duplicated gene symbols; var_names_make_unique
        # on the next line is precisely the remedy, so the warning is noise.
        warnings.filterwarnings("ignore", message="Variable names are not unique")
        adata = sc.read_10x_h5(path)
    adata.var_names_make_unique()
    adata.obs_names = [f"cell_{i + 1}" for i in range(adata.n_obs)]
    cap = min(int(n_cells_cap), adata.n_obs)
    return adata[:cap].copy()
