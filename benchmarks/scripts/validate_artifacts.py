"""Validate every committed benchmark artifact under benchmark-results/.

Checks per JSON artifact: required top-level keys, identical JSON round-trip
(via the harness writer's default=str serializer), non-negative timing when
present, and schema_version. Exits 1 with a report if anything fails.

Usage: python3 benchmarks/scripts/validate_artifacts.py
"""

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from benchmarks import _measure  # noqa: E402

REQUIRED = [
    "schema_version",
    "benchmark",
    "revision",
    "dirty",
    "environment",
    "input",
    "implementation",
    "timing_seconds",
    "peak_rss_bytes",
    "correctness",
    "warnings",
    "pbmc",
]
#: Package-2+ artifacts add the registry fixture record; Package-1 artifacts
#: legitimately predate it, so the key is validated only when present.
OPTIONAL_KEYS = ("fixture",)
ROOT = REPO / "benchmark-results"


def main() -> int:
    """Validate all committed artifacts; nonzero exit on failure."""
    artifacts = sorted(ROOT.rglob("*.json"))
    if not artifacts:
        print("no artifacts found under benchmark-results/")
        return 1
    failures = []
    checked = 0
    aux = 0
    for path in artifacts:
        data = json.loads(path.read_text())
        checked += 1
        if "schema_version" not in data:
            # Auxiliary measurement records (e.g. stage4-threads JSONs) are
            # validated only for round-trip stability, not the harness schema.
            aux += 1
            rt = json.loads(json.dumps(data, default=str)) == data
            if not rt:
                failures.append(f"{path.relative_to(REPO)}: JSON round-trip differs")
            continue
        missing = [k for k in REQUIRED if k not in data]
        if missing:
            failures.append(f"{path.relative_to(REPO)}: missing keys {missing}")
        for key in OPTIONAL_KEYS:
            if key in data and not (data[key] is None or isinstance(data[key], dict)):
                failures.append(f"{path.relative_to(REPO)}: {key} must be dict or null")
        if data.get("schema_version") != _measure.SCHEMA_VERSION:
            failures.append(f"{path.relative_to(REPO)}: schema_version {data.get('schema_version')}")
        rt = json.loads(json.dumps(data, default=str)) == data
        if not rt:
            failures.append(f"{path.relative_to(REPO)}: JSON round-trip differs")
        med = (data.get("timing_seconds") or {}).get("median")
        if med is not None and (not isinstance(med, (int, float)) or med < 0):
            failures.append(f"{path.relative_to(REPO)}: invalid median {med!r}")
    if failures:
        print(f"FAIL {len(failures)} of {checked} records:")
        for f in failures:
            print(" -", f)
        return 1
    print(
        f"ok: {checked} records validate ({checked - aux} harness artifacts, {aux} auxiliary); round-trip, schema, required keys"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
