"""Run the Package-3 cross-dataset sweep (tier-m-real + tier-h2h-qc per dataset).

Sequential runner writing artifacts under benchmark-results/package3/.
Usage: python3 benchmarks/scripts/run_sweep.py [--reps N] [--profile P]
"""

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from benchmarks import _fixtures, _measure  # noqa: E402


def main() -> int:
    """Run the sweep and report failures."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--profile", default="S")
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--datasets", default=",".join(_fixtures.known_datasets()))
    args = parser.parse_args()
    out_root = REPO / "benchmark-results" / "package3"
    datasets = [d for d in args.datasets.split(",") if d in _fixtures.FIXTURES]

    failures = []
    for dataset in datasets:
        for spec, benchmark, impls in (
            ("real", "tier-m-real", ("current",)),
            ("h2h", "tier-h2h-qc", ("baseline", "current")),
        ):
            for impl in impls:
                out = out_root / spec / dataset
                out.mkdir(parents=True, exist_ok=True)
                print(f"[sweep] {benchmark} {impl} dataset={dataset}", flush=True)
                try:
                    artifact = _measure.run_benchmark(
                        benchmark=benchmark,
                        implementation=impl,
                        profile=args.profile,
                        seed=0,
                        repetitions=args.reps,
                        timeout=args.timeout,
                        output_dir=str(out),
                        dataset=dataset,
                    )
                    status = artifact["correctness"]["status"]
                    med = artifact["timing_seconds"]["median"]
                    print(f"[sweep]   -> {status} median={med:.3f}s warnings={len(artifact['warnings'])}", flush=True)
                    if status != "pass" or artifact["warnings"]:
                        failures.append(
                            f"{benchmark}/{impl}/{dataset}: status={status} warnings={artifact['warnings']}"
                        )
                except Exception as exc:  # noqa: BLE001
                    print(f"[sweep]   -> FAILED: {type(exc).__name__}: {exc}", flush=True)
                    failures.append(f"{benchmark}/{impl}/{dataset}: {type(exc).__name__}: {exc}")
    if failures:
        print("\nSWEEP FAILURES:")
        for f in failures:
            print(" -", f)
        return 1
    print("\nsweep complete: all families pass with zero warnings")
    return 0


if __name__ == "__main__":
    sys.exit(main())
