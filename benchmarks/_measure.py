"""Spawn-isolated measurement harness for the Package-1 benchmarks.

Every timed repetition runs in a fresh ``multiprocessing`` spawn child so that
elapsed time and process peak RSS (``resource.getrusage``) include native
NumPy/SciPy/igraph allocations without cross-repetition high-water
contamination. Warm-up happens *inside* the child before timing, mutable
inputs are regenerated deterministically from the seed (never shared between
samples), and the candidate production helpers are imported lazily inside the
kernels so the baseline harness runs before any optimization exists.

Artifacts (JSON + Markdown) live only under the gitignored
``benchmark-results/`` directory.
"""

from __future__ import annotations

import hashlib
import importlib
import importlib.metadata
import json
import multiprocessing
import os
import platform
import statistics
import subprocess
import sys
import time
import traceback
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import numpy as np

from benchmarks import _fixtures

SCHEMA_VERSION = 1

REPO_ROOT = Path(__file__).resolve().parents[1]
PBMC_EXPECTED = "tests/data/pbmc4k/raw.h5"

DEFAULT_REPETITIONS = {
    "cluster-means": 7,
    "graph-conversion": 7,
    "tier-s": 3,
    "tier-m-pbmc": 3,
    "tier-m-real": 3,
    "tier-h2h-qc": 3,
}
KERNEL_NAMES = tuple(DEFAULT_REPETITIONS)

THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)

_CORRECTNESS_TOL = {"rtol": 1e-12, "atol": 1e-12}


# ---------------------------------------------------------------------------
# Git / environment / fingerprint helpers
# ---------------------------------------------------------------------------


def _git_state() -> dict:
    def _run(*args: str) -> str | None:
        try:
            out = subprocess.run(args, cwd=REPO_ROOT, capture_output=True, text=True, timeout=10)
            return out.stdout.strip() if out.returncode == 0 else None
        except Exception:  # noqa: BLE001
            return None

    sha = _run("git", "rev-parse", "HEAD")
    dirty = sha is not None and bool(_run("git", "status", "--porcelain"))
    return {"revision": sha or "unavailable", "dirty": dirty if sha is not None else True}


def _jsonable(value):
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return str(value)


def _physical_cpu_count() -> int | None:
    try:
        import psutil  # optional, not a project dependency

        return psutil.cpu_count(logical=False)
    except Exception:  # noqa: BLE001
        pass
    if sys.platform == "darwin":
        try:
            out = subprocess.run(["sysctl", "-n", "hw.physicalcpu"], capture_output=True, text=True, timeout=5)
            return int(out.stdout.strip()) if out.returncode == 0 and out.stdout.strip() else None
        except Exception:  # noqa: BLE001
            return None
    if sys.platform.startswith("linux"):
        try:
            ids = set()
            for line in Path("/proc/cpuinfo").read_text().splitlines():
                if line.startswith("physical id"):
                    ids.add(line.split(":", 1)[1].strip())
            return len(ids) or None
        except Exception:  # noqa: BLE001
            return None
    return None


def _available_memory_bytes() -> int | None:
    try:
        import psutil

        return int(psutil.virtual_memory().total)
    except Exception:  # noqa: BLE001
        pass
    if sys.platform == "darwin":
        try:
            out = subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True, timeout=5)
            return int(out.stdout.strip()) if out.returncode == 0 and out.stdout.strip() else None
        except Exception:  # noqa: BLE001
            return None
    if sys.platform.startswith("linux"):
        try:
            for line in Path("/proc/meminfo").read_text().splitlines():
                if line.startswith("MemTotal:"):
                    return int(line.split()[1]) * 1024
        except Exception:  # noqa: BLE001
            return None
    return None


def _capture_call(func, *args, **kwargs) -> str:
    buf = StringIO()
    try:
        with redirect_stdout(buf):
            result = func(*args, **kwargs)
        text = buf.getvalue() or (result if isinstance(result, str) else "")
    except Exception:  # noqa: BLE001
        return "unavailable"
    return text.strip()[:800] or "unavailable"


def _blas_config() -> dict:
    try:
        cfg = np.show_config(mode="dicts")
    except Exception:  # noqa: BLE001
        return {"source": "unavailable"}
    out = {}
    for section in ("blas", "blas_opt", "lapack", "lapack_opt"):
        if section in cfg:
            out[section] = _jsonable(cfg[section])
    return out


def _threadpool_config() -> str:
    try:
        return _capture_call(np.show_runtime)
    except Exception:  # noqa: BLE001
        return "unavailable"


def _collect_environment() -> dict:
    env = {
        "os": platform.system(),
        "os_release": platform.release(),
        "architecture": platform.machine(),
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "logical_cpus": os.cpu_count(),
        "physical_cpus": _physical_cpu_count(),
        "available_memory_bytes": _available_memory_bytes(),
    }
    for name in ("numpy", "scipy", "sklearn", "anndata", "igraph"):
        try:
            env[f"{name}_version"] = importlib.metadata.version(name)
        except Exception:  # noqa: BLE001
            env[f"{name}_version"] = "unavailable"
    env["blas_config"] = _blas_config()
    env["blas_threadpool"] = _threadpool_config()
    env["thread_environment"] = {k: os.environ[k] for k in THREAD_ENV_VARS if k in os.environ} or None
    return env


def _sha256(*parts: bytes) -> str:
    digest = hashlib.sha256()
    for part in parts:
        digest.update(part)
    return digest.hexdigest()[:16]


def fingerprint_matrix(matrix) -> str:
    m = matrix.tocsr()
    return _sha256(
        repr(m.shape).encode(),
        str(m.dtype).encode(),
        np.asarray(m.indptr).tobytes(),
        np.asarray(m.indices).tobytes(),
        np.asarray(m.data).tobytes(),
    )


def fingerprint_array(array: np.ndarray) -> str:
    return _sha256(repr(array.shape).encode(), str(array.dtype).encode(), np.ascontiguousarray(array).tobytes())


def _file_hashes(paths) -> dict:
    out = {}
    for path in sorted(paths):
        p = Path(path)
        if p.exists():
            out[str(p.relative_to(REPO_ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()[:16]
    return out


def _implementation_record(benchmark: str, implementation: str) -> dict:
    tl = REPO_ROOT / "src" / "validrops" / "tl"
    pp = REPO_ROOT / "src" / "validrops" / "pp"
    sources = [REPO_ROOT / "src" / "validrops" / "_constants.py"]
    if benchmark in ("cluster-means", "tier-s", "tier-m-pbmc", "tier-m-real", "tier-h2h-qc"):
        sources.append(tl / "expression_metrics.py")
    if benchmark in ("graph-conversion", "tier-s", "tier-m-pbmc", "tier-m-real", "tier-h2h-qc"):
        sources.append(tl / "_snn.py")
    if benchmark in ("tier-m-real", "tier-h2h-qc"):
        sources.append(tl / "quality_metrics.py")
        sources.append(tl / "_deviance.py")
    if benchmark == "tier-h2h-qc":
        sources.extend(
            [
                pp / "rank_barcodes.py",
                pp / "quality_filter.py",
                tl / "_segmented.py",
                tl / "_stats.py",
                tl / "_uik.py",
                REPO_ROOT / "src" / "validrops" / "_pipeline.py",
            ]
        )
    sources.extend([REPO_ROOT / "benchmarks" / f"{name}.py" for name in ("_datasets", "_fixtures", "_measure", "run")])
    return {
        "label": implementation,
        "source_hashes": _file_hashes(sources),
        "harness_hashes": _file_hashes(sources[-4:]),
    }


def resolve_pbmc_state(params: dict) -> dict:
    """PBMC4K fixture detection rooted at the repository, not the CWD."""
    requested = params.get("pbmc_path") or str(REPO_ROOT / PBMC_EXPECTED)
    override = params.get("pbmc_override", "auto")
    exists = os.path.isfile(requested)
    if override == "yes" and not exists:
        raise FileNotFoundError(f"PBMC4K fixture requested (--pbmc-override yes) but missing: {requested}")
    available = exists and override != "no"
    return {"available": available, "path": requested, "override": override}


# ---------------------------------------------------------------------------
# Benchmark kernels (production imports are lazy so baseline runs pre-refactor)
# ---------------------------------------------------------------------------


def _sweep_list_for_search(k_min: int, coarse_mins: dict[float, int]) -> list[float]:
    """Fine-resolution sweep exactly as ``expression_metrics._deep_clustering``."""
    closest = min(abs(v - k_min) for v in coarse_mins.values())
    close_res = [r for r, v in coarse_mins.items() if abs(v - k_min) == closest]
    fine: list[float] = []
    for r in close_res:
        for delta in np.arange(-0.9, 0.9 + 1e-9, 0.1):
            value = r + delta
            if value not in fine:
                fine.append(value)
    return fine


def _deep_search(louvain_call, k_min: int, random_state: int) -> tuple[np.ndarray, float, list[float]]:
    """Mirror of ``expression_metrics._deep_clustering``.

    The production copy in ``expression_metrics.py`` is authoritative; this
    exists so the benchmark can time the same resolution search under both the
    wrapper and prepared Louvain paths. Returns
    ``(labels, chosen_resolution, call_trace)`` with the trace in the exact call
    order production makes them.
    """

    def smallest(resolution: float) -> tuple[int, np.ndarray]:
        labels = louvain_call(resolution)
        trace.append(resolution)
        return int(np.bincount(labels).min()), labels

    trace: list[float] = []
    coarse = {float(r): smallest(float(r))[0] for r in range(1, 21)}
    fine = _sweep_list_for_search(k_min, coarse)
    fine_pairs = [smallest(res) for res in fine]
    fine_mins = {r: m for r, (m, _) in zip(fine, fine_pairs, strict=True)}
    exact = [r for r, v in fine_mins.items() if v == k_min]
    if not exact:
        target = min(fine_mins.values(), key=lambda v: abs(v - k_min))
        exact = [r for r, v in fine_mins.items() if v == target]
    _, labels = fine_pairs[fine.index(max(exact))]
    return labels, float(max(exact)), trace


def _maxrss_bytes() -> int | None:
    """Process high-water RSS of the *current* process (called inside the child).

    macOS ``ru_maxrss`` is in bytes; Linux (and other BSDs) report KiB.
    """
    try:
        import resource
    except ImportError:
        return None
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value) if sys.platform == "darwin" else int(value) * 1024


def _cluster_means_param(payload: dict) -> dict:
    from benchmarks._datasets import make_sparse_counts

    matrix, meta = make_sparse_counts(payload["profile"], seed=payload["seed"])
    return {
        "data": matrix,
        "meta": meta,
        "fingerprint": fingerprint_matrix(matrix),
        "params": {"profile": payload["profile"], "seed": payload["seed"]},
        "impl": payload["implementation"],
    }


def _run_cluster_means_impl(input_: dict, impl: str) -> np.ndarray:
    matrix = input_["data"]
    if impl == "baseline":
        # Dense float64 oracle: cast to float64 BEFORE expm1, accumulate mean in
        # float64 (reviewer-mandated); replaces production's former float32 mean.
        return np.expm1(matrix.toarray().astype(np.float64)).mean(axis=0, dtype=np.float64)
    from validrops.tl.expression_metrics import _sparse_expm1_col_mean

    return _sparse_expm1_col_mean(matrix)


def _verify_cluster_means(input_: dict) -> dict:
    matrix = input_["data"]
    snapshot = (
        np.asarray(matrix.indptr).copy(),
        np.asarray(matrix.indices).copy(),
        np.asarray(matrix.data).copy(),
        matrix.shape,
        matrix.dtype,
    )
    dense_mean = np.expm1(matrix.toarray().astype(np.float64)).mean(axis=0, dtype=np.float64)
    # Independent float64 accumulation order (sum / n vs mean) validates the oracle
    # itself; this is the gate that can run before the sparse helper exists.
    dense_sum = (
        np.asarray(np.expm1(matrix.toarray().astype(np.float64)).sum(axis=0, dtype=np.float64)).ravel()
        / matrix.shape[0]
    )
    oracle_ok = bool(np.allclose(dense_mean, dense_sum, rtol=_CORRECTNESS_TOL["rtol"], atol=_CORRECTNESS_TOL["atol"]))
    metrics = {
        "oracle_self_consistent": oracle_ok,
        "dense_oracle_dtype": str(dense_mean.dtype),
        "n_genes": int(dense_mean.size),
        "allclose_rtol": _CORRECTNESS_TOL["rtol"],
        "allclose_atol": _CORRECTNESS_TOL["atol"],
    }
    try:
        from validrops.tl.expression_metrics import _sparse_expm1_col_mean

        sparse = _sparse_expm1_col_mean(matrix)
        sparse_ok = bool(np.allclose(sparse, dense_mean, rtol=_CORRECTNESS_TOL["rtol"], atol=_CORRECTNESS_TOL["atol"]))
        unchanged = (
            np.array_equal(snapshot[0], matrix.indptr)
            and np.array_equal(snapshot[1], matrix.indices)
            and np.array_equal(snapshot[2], matrix.data)
            and snapshot[3] == matrix.shape
            and snapshot[4] == matrix.dtype
        )
        metrics.update(
            {
                "sparse_candidate": "present",
                "sparse_matches_dense": sparse_ok,
                "max_abs_diff": float(np.max(np.abs(sparse - dense_mean))) if dense_mean.size else 0.0,
                "sparse_mean_dtype": str(sparse.dtype),
                "input_sparse_structure_unchanged": bool(unchanged),
            }
        )
        ok = oracle_ok and sparse_ok and unchanged
    except (ImportError, AttributeError) as exc:
        metrics["sparse_candidate"] = "absent"
        metrics["reason"] = f"candidate helper absent ({exc})"
        ok = oracle_ok
    return {"status": "pass" if ok else "fail", "metrics": metrics}


def _cluster_means_sample(input_: dict) -> dict:
    result = _run_cluster_means_impl(input_, input_["impl"])
    return {
        "elapsed": None,  # filled by the timed wrapper
        "rss_bytes": None,
        "output_fingerprint": fingerprint_array(np.asarray(result)),
        "correctness": _verify_cluster_means(input_),
        "warnings": [],
    }


def _graph_conversion_param(payload: dict) -> dict:
    from benchmarks._datasets import make_clustered_embedding

    p = payload["params"]
    emb, labels, meta = make_clustered_embedding(p["n_cells"], p["n_pcs"], p["n_clusters"], seed=payload["seed"])
    return {
        "data": emb,
        "cluster_truth": labels,
        "meta": meta,
        "fingerprint": fingerprint_array(emb),
        "params": payload["params"],
        "impl": payload["implementation"],
    }


def _build_snn_adjacency(input_: dict):
    from validrops.tl._snn import snn_graph

    p = input_["params"]
    return snn_graph(input_["data"], k=p["knn"], prune=p["prune"])


def _run_graph_search(input_: dict, impl: str) -> tuple[np.ndarray, float, list[float]]:
    p = input_["params"]
    if impl == "current":
        from validrops.tl._snn import louvain_prepared, prepare_louvain_graph

        prepared = prepare_louvain_graph(input_["adjacency"])
        call = lambda res: louvain_prepared(prepared, res, random_state=p["seed"])
    else:
        from validrops.tl._snn import louvain

        call = lambda res: louvain(input_["adjacency"], res, random_state=p["seed"])
    return _deep_search(call, p["k_min"], p["seed"])


def _graph_state(graph) -> dict:
    return {"vcount": graph.vcount(), "ecount": graph.ecount(), "weights": list(graph.es["weight"])}


def _verify_graph_conversion(input_: dict) -> dict:
    """Gate: prepared path == wrapper path over the actual production sweep.

    Also verifies the prepared igraph object is not mutated by the sweep.
    Before the prepared helpers exist (baseline capture), the gate checks
    wrapper-search determinism and sane output instead, and records the
    candidate as absent.
    """
    p = input_["params"]
    from validrops.tl._snn import louvain  # the wrapper always exists

    try:
        from validrops.tl._snn import louvain_prepared, prepare_louvain_graph
    except (ImportError, AttributeError) as exc:
        wrapper_a, res_a, trace_a = _deep_search(
            lambda res: louvain(input_["adjacency"], res, random_state=p["seed"]), p["k_min"], p["seed"]
        )
        wrapper_b, res_b, trace_b = _deep_search(
            lambda res: louvain(input_["adjacency"], res, random_state=p["seed"]), p["k_min"], p["seed"]
        )
        deterministic = bool(np.array_equal(wrapper_a, wrapper_b) and res_a == res_b and trace_a == trace_b)
        metrics = {
            "prepared_candidate": "absent",
            "reason": f"prepared helpers absent ({exc})",
            "n_resolutions_actual_sweep": int(len(trace_a)),
            "wrapper_search_deterministic": deterministic,
            "labels_cover_all_cells": bool(wrapper_a.size == input_["adjacency"].shape[0]),
        }
        ok = deterministic and wrapper_a.size == input_["adjacency"].shape[0]
        return {"status": "pass" if ok else "fail", "metrics": metrics}

    prepared = prepare_louvain_graph(input_["adjacency"])
    before = _graph_state(prepared)
    wrapper_labels, wrapper_res, wrapper_trace = _deep_search(
        lambda res: louvain(input_["adjacency"], res, random_state=p["seed"]), p["k_min"], p["seed"]
    )
    prepared_labels, prepared_res, prepared_trace = _deep_search(
        lambda res: louvain_prepared(prepared, res, random_state=p["seed"]), p["k_min"], p["seed"]
    )
    after = _graph_state(prepared)
    labels_match = bool(np.array_equal(wrapper_labels, prepared_labels)) and wrapper_res == prepared_res
    state_unchanged = before == after and prepared.vcount() == input_["adjacency"].shape[0]
    metrics = {
        "prepared_candidate": "present",
        "n_resolutions_actual_sweep": int(len(wrapper_trace)),
        "wrapper_trace_len": int(len(wrapper_trace)),
        "prepared_trace_len": int(len(prepared_trace)),
        "labels_equal_after_relabel": bool(labels_match),
        "chosen_resolution_equal": bool(wrapper_res == prepared_res),
        "graph_vcount_ecount_unchanged": bool(state_unchanged),
        "trace_lists_equal": bool(wrapper_trace == prepared_trace),
    }
    ok = labels_match and state_unchanged and wrapper_trace == prepared_trace
    return {"status": "pass" if ok else "fail", "metrics": metrics}


def _graph_conversion_sample(input_: dict) -> dict:
    labels, resolution, trace = _run_graph_search(input_, input_["impl"])
    fingerprint = _sha256(
        np.asarray(labels).tobytes(),
        repr(resolution).encode(),
        repr(len(trace)).encode(),
    )
    return {
        "elapsed": None,
        "rss_bytes": None,
        "output_fingerprint": fingerprint,
        "correctness": _verify_graph_conversion(input_),
        "warnings": [],
    }


def _tier_s_param(payload: dict) -> dict:
    from benchmarks._datasets import make_sparse_counts

    matrix, meta = make_sparse_counts(payload["profile"], seed=payload["seed"])
    return {
        "data": matrix,
        "meta": meta,
        "fingerprint": fingerprint_matrix(matrix),
        "params": payload["params"],
        "impl": payload["implementation"],
    }


def _build_tier_s_adata(input_: dict) -> object:
    import anndata as ad

    meta = input_["meta"]
    matrix = input_["data"].copy()
    var_names = list(meta["mitochondrial"]) + list(meta["ribosomal"])
    # Every remaining position is a unique base gene name; rebuild names so the
    # mito/ribo sets are real columns of the matrix.
    names = [f"gene_{i:06d}" for i in range(matrix.shape[1])]
    names[: len(var_names)] = var_names
    adata = ad.AnnData(X=matrix, var={"gene_names": names})
    adata.var_names = names
    adata.obs["qc_pass"] = True
    adata.uns["validrops"] = {
        "gene_sets": {
            "mitochondrial": list(meta["mitochondrial"]),
            "ribosomal": list(meta["ribosomal"]),
        }
    }
    return adata


def _run_tier_s(input_: dict) -> dict:
    import validrops

    p = input_["params"]
    adata = _build_tier_s_adata(input_)
    validrops.tl.expression_metrics(
        adata,
        nfeats=p["nfeats"],
        npcs=p["npcs"],
        k_min=p["k_min"],
        res_shallow=p["res_shallow"],
        top_n=p["top_n"],
        random_state=p["seed"],
    )
    return _tier_s_summary(adata)


def _verify_tier_s(input_: dict) -> dict:
    p = input_["params"]
    adata = _build_tier_s_adata(input_)
    import validrops

    validrops.tl.expression_metrics(
        adata,
        nfeats=p["nfeats"],
        npcs=p["npcs"],
        k_min=p["k_min"],
        res_shallow=p["res_shallow"],
        top_n=p["top_n"],
        random_state=p["seed"],
    )
    obs = adata.obs
    stats = adata.uns["validrops"]["cluster_stats"]
    assigned = (~obs["cluster"].isna()).all()
    non_empty = stats.shape[0] >= 1
    numeric = stats[[c for c in stats.columns if c not in ("cluster",) if c in stats.columns]]
    finite = bool(np.isfinite(numeric.to_numpy(dtype=float)).all()) if numeric.shape[0] else True
    metrics = {
        "all_cells_clustered": bool(assigned),
        "stats_non_empty": bool(non_empty),
        "stats_finite": finite,
        "n_clusters": int(obs["cluster"].nunique()),
        "n_shallow_clusters": int(obs["cluster_shallow"].nunique()),
        "stats_shape": [int(stats.shape[0]), int(stats.shape[1])],
        "workload": "Stage-3 expression_metrics on deterministic synthetic data"
        + ("" if p["profile"] != "tier-m-pbmc" else "; real PBMC4K subset"),
    }
    ok = assigned and non_empty
    return {"status": "pass" if ok else "fail", "metrics": metrics}


def _tier_s_sample(input_: dict) -> dict:
    out = _run_tier_s(input_)
    return {
        "elapsed": None,
        "rss_bytes": None,
        "output_fingerprint": out["cluster_fingerprint"] + f"-{out['n_clusters']}",
        "correctness": _verify_tier_s(input_),
        "warnings": [],
    }


def _pbmc_param(payload: dict) -> dict:
    state = resolve_pbmc_state(payload)
    return {
        "data": state,
        "meta": {"cells_cap": payload["params"]["n_cells_cap"], "source": "real PBMC4K raw.h5"},
        "fingerprint": _sha256(state["path"].encode()),
        "params": payload["params"],
        "impl": payload["implementation"],
    }


def _build_pbmc_adata(input_: dict):
    import warnings

    import scanpy as sc

    p = input_["params"]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Variable names are not unique")
        adata = sc.read_10x_h5(input_["data"]["path"])
    adata.var_names_make_unique()
    adata.obs_names = [f"cell_{i + 1}" for i in range(adata.n_obs)]
    adata = adata[: min(p["n_cells_cap"], adata.n_obs)].copy()
    adata.obs["qc_pass"] = True
    return adata


def _run_pbmc_tier(input_: dict) -> dict:
    import validrops

    p = input_["params"]
    adata = _build_pbmc_adata(input_)
    # real gene sets; annotation detection runs on real symbols
    validrops.tl.quality_metrics(adata)
    validrops.tl.expression_metrics(adata, nfeats=p["nfeats"], npcs=p["npcs"], random_state=p["seed"])
    return _tier_s_summary(adata)


def _tier_s_summary(adata) -> dict:
    obs = adata.obs
    stats = adata.uns["validrops"]["cluster_stats"]
    cluster_codes = obs["cluster"].cat.codes.to_numpy()
    return {
        "cluster_fingerprint": _sha256(cluster_codes.astype(np.int64).tobytes()),
        "n_clusters": int(cluster_codes.max()) + 1,
        "n_cluster_shallow": int(obs["cluster_shallow"].nunique()),
        "stats_shape": [int(stats.shape[0]), int(stats.shape[1])],
        "n_obs": int(obs.shape[0]),
    }


def _verify_pbmc_tier(input_: dict) -> dict:
    import validrops

    p = input_["params"]
    adata = _build_pbmc_adata(input_)
    validrops.tl.quality_metrics(adata)
    validrops.tl.expression_metrics(adata, nfeats=p["nfeats"], npcs=p["npcs"], k_min=p["k_min"], random_state=p["seed"])
    summary = _tier_s_summary(adata)
    stats = adata.uns["validrops"]["cluster_stats"]
    import importlib as _importlib

    em_module = _importlib.import_module("validrops.tl.expression_metrics")
    metrics = {
        "workload": "Stage-3 expression_metrics on real PBMC4K (bounded subset)",
        "all_cells_clustered": bool((~adata.obs["cluster"].isna()).all()),
        "stats_frame_columns": list(stats.columns) == list(em_module.STAT_COLUMNS),
        "n_cluster_rows": int(stats.shape[0]),
        "n_clusters": summary["n_clusters"],
        "stats_shape": summary["stats_shape"],
        "n_obs": summary["n_obs"],
        "species_detected": adata.uns["validrops"].get("species"),
    }
    ok = metrics["all_cells_clustered"] and metrics["n_clusters"] >= 1
    return {"status": "pass" if ok else "fail", "metrics": metrics}


def _pbmc_sample(input_: dict) -> dict:
    out = _run_pbmc_tier(input_)
    return {
        "elapsed": None,
        "rss_bytes": None,
        "output_fingerprint": out["cluster_fingerprint"] + f"-{out['n_clusters']}",
        "correctness": _verify_pbmc_tier(input_),
        "warnings": [],
    }


# ---------------------------------------------------------------------------
# tier-m-real: generic real-dataset Stage-3 workload (fixture registry)
# ---------------------------------------------------------------------------


def _real_param(payload: dict) -> dict:
    state = _fixtures.resolve_fixture(payload["dataset"], payload["pbmc_override"], payload["pbmc_path"])
    p = payload["params"]
    return {
        "data": state,
        "meta": {
            "cells_cap": p.get("n_cells_cap"),
            "source": state["meta"]["tissue"],
            "species": state["meta"]["species"],
            "chemistry": state["meta"]["chemistry"],
            "license": state["meta"]["license"],
            "fixture_cells": state["meta"]["n_cells"],
            "fixture_genes": state["meta"]["n_genes"],
        },
        "fingerprint": _sha256(state["path"].encode()),
        "params": p,
        "impl": payload["implementation"],
    }


def _build_real_adata(input_: dict):
    p = input_["params"]
    adata = _fixtures.load_fixture_adata(input_["data"]["path"], p["n_cells_cap"])
    adata.obs["qc_pass"] = True
    return adata


def _run_real_tier(input_: dict) -> dict:
    import validrops

    p = input_["params"]
    adata = _build_real_adata(input_)
    validrops.tl.quality_metrics(adata)
    validrops.tl.expression_metrics(adata, nfeats=p["nfeats"], npcs=p["npcs"], random_state=p["seed"])
    return _tier_s_summary(adata)


def _verify_real_tier(input_: dict) -> dict:
    import importlib

    import validrops

    p = input_["params"]
    adata = _build_real_adata(input_)
    validrops.tl.quality_metrics(adata)
    validrops.tl.expression_metrics(adata, nfeats=p["nfeats"], npcs=p["npcs"], k_min=p["k_min"], random_state=p["seed"])
    summary = _tier_s_summary(adata)
    stats = adata.uns["validrops"]["cluster_stats"]
    em_module = importlib.import_module("validrops.tl.expression_metrics")
    metrics = {
        "workload": f"Stage-3 expression_metrics on {p['dataset']} (bounded subset)",
        "all_cells_clustered": bool((~adata.obs["cluster"].isna()).all()),
        "stats_frame_columns": list(stats.columns) == list(em_module.STAT_COLUMNS),
        "n_clusters": summary["n_clusters"],
        "stats_shape": summary["stats_shape"],
        "n_obs": summary["n_obs"],
        "species_detected": adata.uns["validrops"].get("species"),
        "dataset": p["dataset"],
    }
    ok = metrics["all_cells_clustered"] and metrics["n_clusters"] >= 1
    return {"status": "pass" if ok else "fail", "metrics": metrics}


def _real_tier_sample(input_: dict) -> dict:
    out = _run_real_tier(input_)
    return {
        "elapsed": None,
        "rss_bytes": None,
        "output_fingerprint": out["cluster_fingerprint"] + f"-{out['n_clusters']}",
        "correctness": _verify_real_tier(input_),
        "warnings": [],
    }


# ---------------------------------------------------------------------------
# tier-h2h-qc: validrops vs Seurat-style default QC on an identical fixture
# ---------------------------------------------------------------------------
#
# The reference side re-implements the classic, published Seurat QC defaults
# (nFeature/nCount minima and a mitochondrial-fraction cap) deterministically
# in Python so the family runs anywhere without an R runtime. The thresholds
# are recorded in every artifact; benchmarks/R/h2h_seurat_crosscheck.R
# cross-validates this baseline against real Seurat 5.5.1 on the same matrix.


def _seurat_style_qc(
    adata,
    mito_cap: float,
    min_features: int,
    min_counts: int,
    min_cells: int = 3,
) -> np.ndarray:
    """Deterministic re-implementation of Seurat's default QC subset.

    Mirrors what ``CreateSeuratObject`` + ``subset()`` achieve with the
    classic defaults: genes present in fewer than ``min_cells`` cells are
    dropped before ANY cell-level metric is computed (Seurat's ``min.cells``
    default is 3; ``nCount_RNA``/``nFeature_RNA``/``percent.mt`` are all
    calculated on the gene-filtered object). A barcode passes iff the
    mitochondrial fraction, feature count, and UMI count meet the caps. Mito
    symbols are matched by the standard ``MT-/Mt-/mt-`` prefixes
    (human/mouse), the same families validrops detects in
    ``quality_metrics``. Cross-checked against Seurat 5.x in
    ``benchmarks/R/h2h_seurat_crosscheck.R``.
    """
    X = adata.X
    gene_presence = X > 0
    keep_genes = np.asarray(gene_presence.sum(axis=0)).ravel() >= min_cells
    mito_mask = np.asarray(adata.var_names.str.match(r"^(MT|Mt|mt)-"))
    if keep_genes.any() and not keep_genes.all():
        X = X[:, keep_genes]
        mito_mask = mito_mask[keep_genes]
    totals = np.asarray(X.sum(axis=1, dtype=np.float64)).ravel()
    features = np.asarray((X > 0).sum(axis=1)).ravel().astype(np.float64)
    mito_frac = np.zeros(adata.n_obs, dtype=np.float64)
    if mito_mask.any():
        mito_totals = np.asarray(X[:, mito_mask].sum(axis=1, dtype=np.float64)).ravel()
        mito_frac = mito_totals / np.where(totals > 0, totals, 1.0)
    return ((mito_frac <= mito_cap) & (features >= min_features) & (totals >= min_counts)).astype(bool)


def _validrops_qc(adata, random_state: int) -> np.ndarray:
    """Stage 1 + 2 QC callouts exactly as ``validrops()`` runs them."""
    import validrops

    validrops.pp.rank_barcodes(adata, random_state=random_state)
    validrops.tl.quality_metrics(adata)
    validrops.pp.quality_filter(adata, random_state=random_state)
    return adata.obs["qc_pass"].to_numpy(dtype=bool)


def _h2h_masks(input_: dict) -> tuple[np.ndarray, np.ndarray, dict]:
    """Run both QC pipelines on independent copies; returns (validrops, seurat, info)."""
    base = _build_real_adata(input_)
    v_adata = base.copy()
    v_mask = _validrops_qc(v_adata, input_["params"]["seed"])
    s_mask = _seurat_style_qc(
        base.copy(),
        input_["params"]["seurat_mito_cap"],
        input_["params"]["seurat_min_features"],
        input_["params"]["seurat_min_counts"],
    )
    uns = v_adata.uns.get("validrops", {})
    info = {
        "validrops_mito_threshold": float(uns.get("mitochondrial_threshold", float("nan"))),
        "validrops_mito_method": uns.get("mito_threshold_method"),
    }
    return v_mask, s_mask, info


def _agreement_metrics(v_mask: np.ndarray, s_mask: np.ndarray, n_total: int) -> dict:
    kept_v = int(v_mask.sum())
    kept_s = int(s_mask.sum())
    jaccard = float(
        np.intersect1d(np.flatnonzero(v_mask), np.flatnonzero(s_mask)).size
        / max(1, np.union1d(np.flatnonzero(v_mask), np.flatnonzero(s_mask)).size)
    )
    return {
        "n_total": n_total,
        "n_kept_validrops": kept_v,
        "n_kept_seurat_style": kept_s,
        "fraction_kept_validrops": kept_v / n_total,
        "fraction_kept_seurat_style": kept_s / n_total,
        "jaccard_kept_sets": jaccard,
        "kept_set_difference": int(np.logical_xor(v_mask, s_mask).sum()),
    }


def _run_h2h_impl(input_: dict, impl: str) -> np.ndarray:
    adata = _build_real_adata(input_)
    if impl == "baseline":
        p = input_["params"]
        return _seurat_style_qc(adata, p["seurat_mito_cap"], p["seurat_min_features"], p["seurat_min_counts"])
    return _validrops_qc(adata, input_["params"]["seed"])


def _verify_h2h(input_: dict) -> dict:
    p = input_["params"]
    v_mask, s_mask, info = _h2h_masks(input_)
    metrics = {
        "workload": f"QC head-to-head on {p['dataset']} (bounded subset)",
        "dataset": p["dataset"],
        "baseline": "Seurat-style default thresholds (n_features>=200, n_counts>=500, mito<=0.20, min_cells=3)",
        **info,
        **_agreement_metrics(v_mask, s_mask, int(v_mask.size)),
    }
    ok = (
        metrics["n_kept_validrops"] >= 1
        and metrics["n_kept_seurat_style"] >= 1
        and np.isfinite(metrics["jaccard_kept_sets"])
        and 0.0 <= metrics["jaccard_kept_sets"] <= 1.0
    )
    return {"status": "pass" if ok else "fail", "metrics": metrics}


def _h2h_sample(input_: dict) -> dict:
    mask = _run_h2h_impl(input_, input_["impl"])
    return {
        "elapsed": None,
        "rss_bytes": None,
        "output_fingerprint": _sha256(np.asarray(mask, dtype=np.uint8).tobytes(), repr(int(mask.sum())).encode()),
        "correctness": _verify_h2h(input_),
        "warnings": [],
    }


# Dispatch tables -----------------------------------------------------------

_INPUT_BUILDERS = {
    "cluster-means": _cluster_means_param,
    "graph-conversion": _graph_conversion_param,
    "tier-s": _tier_s_param,
    "tier-m-pbmc": _pbmc_param,
    "tier-m-real": _real_param,
    "tier-h2h-qc": _real_param,
}

_SAMPLE_RUNNERS = {
    "cluster-means": _cluster_means_sample,
    "graph-conversion": _graph_conversion_sample,
    "tier-s": _tier_s_sample,
    "tier-m-pbmc": _pbmc_sample,
    "tier-m-real": _real_tier_sample,
    "tier-h2h-qc": _h2h_sample,
}

_VERIFIERS = {
    "cluster-means": _verify_cluster_means,
    "graph-conversion": _verify_graph_conversion,
    "tier-s": _verify_tier_s,
    "tier-m-pbmc": _verify_pbmc_tier,
    "tier-m-real": _verify_real_tier,
    "tier-h2h-qc": _verify_h2h,
}


def _child_prepare(payload: dict) -> dict:
    """Reproduce the deterministic input inside the spawned child."""
    builder = _INPUT_BUILDERS[payload["benchmark"]]
    input_ = builder(payload)
    if payload["benchmark"] == "graph-conversion":
        input_["adjacency"] = _build_snn_adjacency(input_)
    return input_


def _execute_sample(payload: dict) -> dict:
    """One warm-up plus one timed run of the target implementation, in process."""
    input_ = _child_prepare(payload)
    runner = _SAMPLE_RUNNERS[payload["benchmark"]]
    verifier = _VERIFIERS[payload["benchmark"]]

    verify = verifier(input_)  # untimed correctness gate (may be 'skipped')
    out_warm = runner(input_)  # untimed warm-up inside the child
    t0 = time.perf_counter()
    out = runner(input_)
    elapsed = time.perf_counter() - t0
    rss = _maxrss_bytes()

    deterministic = bool(out_warm["output_fingerprint"] == out["output_fingerprint"])
    warnings = list(out["warnings"])
    if not deterministic:
        # Pre-existing production behavior, NOT a Package-1 regression: scipy's
        # svds draws its ARPACK starting vector from OS entropy when v0 is
        # unset (scipy.sparse.linalg.svds), so real (ill-conditioned) inputs
        # can land in v0-dependent near-degenerate subspaces run to run.
        # Synthetic Tier-S inputs are well-conditioned and stay deterministic.
        # Dense SVD is explicitly out of Package-1 scope and unchanged here.
        warnings.append(
            "warm-up vs timed output fingerprints differ: pre-existing scipy svds "
            "OS-entropy starting vector (expression_metrics._embed); dense SVD is "
            "unchanged and out of Package-1 scope"
        )
    return {
        "elapsed": float(elapsed),
        "rss_bytes": rss,
        "output_fingerprint": out["output_fingerprint"],
        "deterministic_warmup_timed": deterministic,
        "correctness": verify,
        "warnings": warnings,
    }


def _child_entry(payload: dict, queue):
    """Spawn-child entry point (top-level so it pickles)."""
    try:
        queue.put({"ok": True, "result": _execute_sample(payload)})
    except Exception:  # noqa: BLE001
        queue.put({"ok": False, "error": traceback.format_exc()})


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def _bounded_stage_params(profile: str) -> dict:
    from benchmarks._datasets import PROFILES

    genes = PROFILES[profile]["genes"]
    nfeats = min(300, genes)
    if nfeats < 20:  # smoke is tiny; keep the embedding meaningful
        nfeats = max(nfeats, 20)
    return {
        "nfeats": nfeats,
        "npcs": 10,
        "k_min": 5,
        "res_shallow": 0.1,
        "top_n": 10,
    }


def _family_params(benchmark: str, profile: str, dataset: str | None = None) -> dict:
    from benchmarks._datasets import PROFILES

    if benchmark == "cluster-means":
        return {"profile": profile}
    if benchmark == "graph-conversion":
        cells = PROFILES[profile]["cells"]
        return {
            "profile": profile,
            "n_cells": cells,
            "n_pcs": 10,
            "n_clusters": 4 if cells >= 20 else 2,
            "knn": 20,
            "prune": 1 / 15,
            "k_min": 5,
        }
    if benchmark == "tier-s":
        return {"profile": profile, **_bounded_stage_params(profile), "workload": "Stage-3 synthetic"}
    if benchmark == "tier-m-pbmc":
        # Bounded real-data Stage-3 workload. k_min=150: on a 1200-cell PBMC4K
        # subset the SNN graph fragments to singletons above res ~3, so small
        # k_min targets are unachievable and the faithful resolution search
        # explodes into ~1000 micro-clusters and hundreds of Louvain calls. A
        # larger target keeps the sweep bounded and the stage meaningful.
        return {
            "profile": profile,
            **_bounded_stage_params(profile),
            "k_min": 150,
            "n_cells_cap": 1200,
            "workload": "Stage-3 PBMC4K",
        }
    if benchmark == "tier-m-real":
        name = dataset or _fixtures.DEFAULT_DATASET
        return {
            "profile": profile,
            **_bounded_stage_params(profile),
            "k_min": 150,
            "n_cells_cap": 1500,
            "dataset": name,
            "workload": f"Stage-3 {name}",
        }
    if benchmark == "tier-h2h-qc":
        name = dataset or _fixtures.DEFAULT_DATASET
        # Same bounded real subset as tier-m-real; the Seurat-style defaults are
        # recorded as parameters so every artifact is self-describing and the
        # baseline is reproducible (cross-checked against Seurat 5.5.1 in
        # benchmarks/R/h2h_seurat_crosscheck.R).
        return {
            "profile": profile,
            "n_cells_cap": 1500,
            "dataset": name,
            "seurat_mito_cap": 0.20,
            "seurat_min_features": 200,
            "seurat_min_counts": 500,
            "seurat_min_cells": 3,
            "workload": f"QC head-to-head {name}",
        }
    raise ValueError(f"unknown benchmark {benchmark!r}")


def run_benchmark(
    benchmark: str,
    implementation: str,
    profile: str,
    seed: int = 0,
    repetitions: int | None = None,
    timeout: float = 120.0,
    output_dir: str | None = None,
    pbmc_override: str = "auto",
    pbmc_path: str | None = None,
    dataset: str | None = None,
) -> dict:
    """Run one benchmark family and return the schema-version-1 artifact dict.

    Each repetition executes in a fresh spawn child (warm-up inside the child,
    then one timed run). ``repetitions=1`` is only allowed for the ``smoke``
    profile (schema/smoke checks). ``dataset`` selects the real fixture for
    ``tier-m-real`` and ``tier-h2h-qc`` (see ``benchmarks/_fixtures.py``).
    """
    if benchmark not in KERNEL_NAMES:
        raise ValueError(f"unknown benchmark {benchmark!r}; choose from {sorted(KERNEL_NAMES)}")
    if implementation not in ("baseline", "current"):
        raise ValueError(f"implementation must be 'baseline' or 'current', got {implementation!r}")
    from benchmarks._datasets import PROFILES

    if profile not in PROFILES:
        raise ValueError(f"unknown profile {profile!r}; choose from {sorted(PROFILES)}")
    if dataset is not None and dataset not in _fixtures.FIXTURES:
        raise ValueError(f"unknown dataset {dataset!r}; choose from {_fixtures.known_datasets()}")

    n_reps = DEFAULT_REPETITIONS[benchmark] if repetitions is None else repetitions
    if n_reps < 1:
        raise ValueError("repetitions must be >= 1")
    if n_reps == 1 and profile != "smoke":
        raise ValueError("--repetitions 1 is only allowed for the smoke profile (schema/smoke checks)")

    params = _family_params(benchmark, profile, dataset)
    params["seed"] = seed  # part of the exact-parameter record; also gates comparisons
    payload = {
        "benchmark": benchmark,
        "implementation": implementation,
        "profile": profile,
        "seed": seed,
        "params": params,
        "pbmc_override": pbmc_override,
        "pbmc_path": pbmc_path,
        "dataset": params.get("dataset"),
    }

    git = _git_state()
    pbmc = resolve_pbmc_state(payload)
    fixture = None
    fixture_expected = None
    if benchmark == "tier-m-pbmc":
        # PBMC4K family keeps its historical hardwired resolution.
        fixture = pbmc
        fixture_expected = PBMC_EXPECTED
    elif params.get("dataset"):
        # Registry-based real-data families resolve their own fixture.
        fixture = _fixtures.resolve_fixture(params["dataset"], pbmc_override, pbmc_path)
        fixture_expected = fixture["path"]
    input_spec = _INPUT_BUILDERS[benchmark](payload)
    input_record = {
        "profile": profile,
        "params": params,
        "fingerprint": input_spec["fingerprint"],
        "meta": input_spec["meta"],
    }

    warnings: list[str] = []
    if fixture is not None and not fixture["available"]:
        label = "PBMC4K" if benchmark == "tier-m-pbmc" else f"dataset {params['dataset']}"
        artifact = {
            "schema_version": SCHEMA_VERSION,
            "benchmark": benchmark,
            "revision": git["revision"],
            "dirty": git["dirty"],
            "environment": _collect_environment(),
            "input": input_record,
            "implementation": _implementation_record(benchmark, implementation),
            "timing_seconds": {
                "median": None,
                "minimum": None,
                "mad": None,
                "repetitions": n_reps,
                "raw": [],
                "coefficient_of_variation": None,
            },
            "peak_rss_bytes": {
                "worker_max": None,
                "worker_min": None,
                "worker_mean": None,
                "per_sample": [],
                "method": "resource.getrusage(RUSAGE_SELF).ru_maxrss inside each spawned child",
                "scope": "per-repetition spawned child process (native allocations included)",
                "unit": "bytes",
                "unavailable_reason": f"benchmark skipped ({label} fixture absent)",
            },
            "correctness": {"status": "skipped", "metrics": {}},
            "warnings": [f"{label} fixture not available; expected at {fixture_expected}"],
            "pbmc": pbmc if benchmark == "tier-m-pbmc" else None,
            "fixture": fixture,
        }
        if output_dir:
            write_artifacts(artifact, output_dir=output_dir)
        return artifact

    ctx = multiprocessing.get_context("spawn")
    samples = []
    for i in range(n_reps):
        queue = ctx.Queue()
        proc = ctx.Process(
            target=_child_entry,
            args=(payload, queue),
            name=f"bench-{benchmark}-{implementation}-{i}",
            daemon=False,
        )
        proc.start()
        proc.join(timeout=timeout)
        if proc.is_alive():
            proc.terminate()
            proc.join(timeout=5)
            if proc.is_alive():
                proc.kill()
                proc.join(timeout=5)
            raise RuntimeError(
                f"benchmark child {i} exceeded the {timeout}s timeout (benchmark={benchmark}, "
                f"implementation={implementation}, profile={profile}); rerun with a larger --timeout"
            )
        try:
            message = queue.get(timeout=30)
        except Exception as exc:  # queue.Empty on py<3.13
            raise RuntimeError(f"benchmark child {i} died without a result: {exc}") from exc
        if not message.get("ok"):
            raise RuntimeError(f"benchmark child {i} failed:\n{message.get('error', 'unknown error')}")
        samples.append(message["result"])

    raw = [s["elapsed"] for s in samples]
    rss_vals = [s["rss_bytes"] for s in samples]
    median = float(statistics.median(raw))
    minimum = float(min(raw))
    deviations = [abs(x - median) for x in raw]
    mad = float(statistics.median(deviations)) if raw else 0.0
    mean = float(np.mean(raw))
    cv = float(np.std(raw, ddof=1) / mean) if (mean > 0 and len(raw) > 1) else None

    correctness_metrics = {}
    combined_status = "pass"
    for i, s in enumerate(samples):
        c = s["correctness"]
        for k, v in c.get("metrics", {}).items():
            correctness_metrics.setdefault(k, v)
        if c["status"] == "fail":
            combined_status = "fail"
            warnings.append(f"correctness gate failed in sample {i}")
        elif c["status"] == "skipped" and combined_status == "pass":
            combined_status = "skipped"
            warnings.append(f"correctness gate skipped in sample {i}: {c.get('reason', 'reason unknown')}")
        warnings.extend(s["warnings"])
    correctness_metrics["deterministic_warmup_timed"] = bool(all(s["deterministic_warmup_timed"] for s in samples))
    correctness_metrics["deterministic_sample_fraction"] = float(
        sum(s["deterministic_warmup_timed"] for s in samples) / n_reps
    )
    if cv is not None and cv > 0.10:
        warnings.append(f"coefficient of variation {cv:.1%} > 10%; treat runtime comparisons as inconclusive")

    artifact = {
        "schema_version": SCHEMA_VERSION,
        "benchmark": benchmark,
        "revision": git["revision"],
        "dirty": git["dirty"],
        "environment": _collect_environment(),
        "input": input_record,
        "implementation": _implementation_record(benchmark, implementation),
        "timing_seconds": {
            "median": median,
            "minimum": minimum,
            "mad": mad,
            "repetitions": n_reps,
            "raw": raw,
            "coefficient_of_variation": cv,
        },
        "peak_rss_bytes": {
            "worker_max": int(max(rss_vals)) if any(v is not None for v in rss_vals) else None,
            "worker_min": int(min(rss_vals)) if any(v is not None for v in rss_vals) else None,
            "worker_mean": float(np.mean([v for v in rss_vals if v is not None]))
            if any(v is not None for v in rss_vals)
            else None,
            "per_sample": rss_vals,
            "method": "resource.getrusage(RUSAGE_SELF).ru_maxrss inside each spawned child",
            "scope": "process high-water including deterministic input generation and warm-up "
            "in the same child; identical inputs across baseline/current isolate the kernel delta",
            "unit": "bytes",
        },
        "correctness": {"status": combined_status, "metrics": correctness_metrics},
        "warnings": warnings,
        "pbmc": pbmc if benchmark == "tier-m-pbmc" else None,
        "fixture": fixture,
    }
    if output_dir:
        write_artifacts(artifact, output_dir=output_dir)
    return artifact


def _fmt(value, spec: str = "") -> str:
    if value is None:
        return "unavailable"
    return format(value, spec)


def _render_markdown(artifact: dict) -> str:
    t = artifact["timing_seconds"]
    rss = artifact["peak_rss_bytes"]
    lines = [
        f"# {artifact['benchmark']} / {artifact['implementation']}",
        "",
        f"- revision: {artifact['revision']} (dirty={artifact['dirty']})",
        f"- profile: {artifact['input']['profile']}",
        f"- input fingerprint: {artifact['input']['fingerprint']}",
        f"- parameters: `{json.dumps(artifact['input']['params'])}`",
        f"- timing seconds: median={_fmt(t['median'], '.6f')}, min={_fmt(t['minimum'], '.6f')}, "
        f"mad={_fmt(t['mad'], '.6f')}, reps={_fmt(t['repetitions'])}, CV={_fmt(t['coefficient_of_variation'], '.4f')}",
        f"- peak RSS bytes: worker_max={_fmt(rss['worker_max'])}, unit={rss.get('unit')}",
        f"- correctness: {artifact['correctness']['status']}",
        f"- environment: {json.dumps(artifact['environment'], default=str)[:400]}",
        "",
    ]
    if artifact["warnings"]:
        lines.append("## Warnings")
        lines.extend(f"- {w}" for w in artifact["warnings"])
    return "\n".join(lines) + "\n"


def write_artifacts(artifact: dict, output_dir: str) -> None:
    """Write `<benchmark>-<implementation>.json` and `.md` under output_dir."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    impl = artifact["implementation"]
    label = impl["label"] if isinstance(impl, dict) else str(impl)
    stem = f"{artifact['benchmark']}-{label}"
    (out / f"{stem}.json").write_text(json.dumps(artifact, indent=2, default=str) + "\n")
    (out / f"{stem}.md").write_text(_render_markdown(artifact))


def compare_artifacts(a: dict, b: dict) -> dict:
    """Compare two artifacts; raises ``ValueError`` on incompatible comparisons."""
    problems = []
    if a["benchmark"] != b["benchmark"]:
        problems.append(f"benchmark {a['benchmark']} vs {b['benchmark']}")
    pa, pb = a["input"]["params"], b["input"]["params"]
    for key in sorted(set(pa) | set(pb)):
        va, vb = pa.get(key), pb.get(key)
        if va != vb:
            problems.append(f"input param {key}: {va!r} vs {vb!r}")
    if a["input"]["profile"] != b["input"]["profile"]:
        problems.append("input profile differs")
    if a["implementation"]["label"] == b["implementation"]["label"] and a is not b:
        problems.append("same implementation label without shared session")
    if a["timing_seconds"]["median"] is None or b["timing_seconds"]["median"] is None:
        problems.append("one artifact was skipped (no timing)")
    if problems:
        raise ValueError("incompatible artifacts: " + "; ".join(problems))
    return {
        "compatible": True,
        "median_speedup": float(a["timing_seconds"]["median"] / b["timing_seconds"]["median"]),
        "candidate": a["implementation"]["label"],
        "reference": b["implementation"]["label"],
        "workspace": "comparison is meaningful for kernels that reject incompatible merges",
    }
