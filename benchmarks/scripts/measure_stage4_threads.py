"""Efficiency B: Stage-4 ridge thread measurement (measurement-only).

The Package-1 ledger deferred a Stage-4 thread question: the consensus ridge
training loop is a ``joblib.Parallel(n_jobs=1)`` port of the R BiocParallel
call, so by DEFAULT there is no joblib parallelism to oversubscribe — each
ridge fit runs in-process. What remains measurable is whether BLAS thread
limits (OMP/OPENBLAS/MKL/VECLIB) change Stage-4 wall time or dead-cell output
for a single worker.

This script runs the full Stage 1-4 pipeline (rank -> quality metrics ->
quality filter -> label_dead training) on a bounded corpus under
OMP_NUM_THREADS=1/2/4, in fresh subprocesses with the vars set BEFORE imports,
and writes JSON artifacts under benchmark-results/package3/stage4-threads/.

Usage: python3 benchmarks/scripts/measure_stage4_threads.py [--reps 2] [--cap 1200] [--dataset pbmc-parotid]
"""

import argparse
import json
import os
import statistics
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from benchmarks import _fixtures  # noqa: E402

THREAD_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)

CHILD = """
import json, os, sys, time
sys.path.insert(0, {repo!r})
from benchmarks import _fixtures
import validrops
import numpy as np

path = {path!r}
cap = {cap}
adata = _fixtures.load_fixture_adata(path, cap)
adata.obs["qc_pass"] = True
t0 = time.perf_counter()
validrops.pp.rank_barcodes(adata, random_state=0)
validrops.tl.quality_metrics(adata)
validrops.pp.quality_filter(adata, random_state=0)
validrops.pp.label_dead(adata, random_state=0)
elapsed = time.perf_counter() - t0
labels = adata.obs["label"].astype(str)
out = {{
    "elapsed": elapsed,
    "n_dead": int((labels == "dead").sum()),
    "n_live": int((labels == "live").sum()),
    "n_uncertain": int((labels == "uncertain").sum()),
    "flag": adata.uns["validrops"].get("label_flag"),
    "threshold": adata.uns["validrops"].get("label_threshold"),
    "qc_pass": int(adata.obs["qc_pass"].sum()),
}}
print(json.dumps(out))
"""


def run_child(threads: int, dataset: str, cap: int, seed_env: dict) -> dict:
    """Run the Stage-1-4 pipeline in a subprocess under a thread limit."""
    env = dict(os.environ)
    for var in THREAD_VARS:
        env[var] = str(threads)
    env.update(seed_env)
    path = _fixtures.fixture_expected_path(dataset)
    code = CHILD.format(repo=str(REPO), path=path, cap=cap)
    proc = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        env=env,
        cwd=str(REPO),
        timeout=1800,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"child failed (threads={threads}):\n{proc.stderr[-2000:]}")
    return json.loads(proc.stdout.strip().splitlines()[-1])


def main() -> int:
    """Run the thread sweep and write per-config artifacts."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=2)
    parser.add_argument("--cap", type=int, default=1200)
    parser.add_argument("--dataset", default="pbmc-parotid")
    args = parser.parse_args()

    if not _fixtures.resolve_fixture(args.dataset, "auto")["available"]:
        print(f"fixture absent; nothing to measure ({_fixtures.fixture_expected_path(args.dataset)})")
        return 1

    out_root = REPO / "benchmark-results" / "package3" / "stage4-threads"
    out_root.mkdir(parents=True, exist_ok=True)
    results = {}
    for threads in (1, 2, 4):
        samples = []
        for rep in range(args.reps):
            print(f"[stage4] threads={threads} rep={rep}", flush=True)
            samples.append(run_child(threads, args.dataset, args.cap, {"VALIDROPS_REP_SEED": str(rep)}))
        elapsed = [s["elapsed"] for s in samples]
        sample = samples[0]
        artifact = {
            "benchmark": "stage4-ridge-threads",
            "dataset": args.dataset,
            "n_cells_cap": args.cap,
            "threads": threads,
            "thread_variables": list(THREAD_VARS),
            "n_jobs": 1,  # Stage-4 consensus loop default (joblib sequential)
            "timing_seconds": {
                "median": statistics.median(elapsed),
                "raw": elapsed,
                "repetitions": args.reps,
            },
            "output": {
                "n_dead": [s["n_dead"] for s in samples],
                "n_live": [s["n_live"] for s in samples],
                "n_uncertain": [s["n_uncertain"] for s in samples],
                "flag": sample["flag"],
                "threshold": sample["threshold"],
                "qc_pass": sample["qc_pass"],
            },
        }
        results[str(threads)] = artifact
        (out_root / f"threads-{threads}.json").write_text(json.dumps(artifact, indent=2) + "\n")
        print(f"[stage4] median {statistics.median(elapsed):.2f}s dead={sample['n_dead']} flag={sample['flag']}")

    med = {t: a["timing_seconds"]["median"] for t, a in results.items()}
    print("\nsummary medians:", json.dumps(med, indent=2))
    dead_sets = {t: set(a["output"]["n_dead"]) for t, a in results.items()}
    print("dead counts per threads:", dead_sets)
    return 0


if __name__ == "__main__":
    sys.exit(main())
