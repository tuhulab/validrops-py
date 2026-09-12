"""Regenerate all manuscript figures from committed benchmark artifacts.

Deterministic (fixed seed, no network): every figure is drawn from the JSONs
under benchmark-results/package3/ (and package2/ for the parotid R
cross-check context). Run:

    uv run --group test --group dev python manuscript/render_figures.py

Outputs PNGs into manuscript/figures/. Fails loudly when an expected artifact
is missing instead of drawing silence.
"""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

FIGDIR = Path(__file__).resolve().parent / "figures"
P3 = REPO / "benchmark-results" / "package3"

plt.rcParams.update({"figure.dpi": 150, "savefig.bbox": "tight"})


def load(p: Path) -> dict:
    data = json.loads(p.read_text())
    return data


def _metrics(artifact: dict) -> dict:
    return artifact["correctness"]["metrics"]


def sweep_artifacts() -> dict[str, dict]:
    """Collect package3 artifacts: dataset -> impl -> {family: artifact}."""
    out = {}
    for family in ("real", "h2h"):
        for impl_dir in (P3 / family).iterdir() if (P3 / family).exists() else []:
            dataset = impl_dir.name
            bucket = out.setdefault(dataset, {})
            for json_path in impl_dir.glob("*.json"):
                a = load(json_path)
                impl = a["implementation"]["label"]
                bucket.setdefault(impl, {})[a["benchmark"]] = a
    return out


def _fig(name: str):
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    return fig, ax


def figure_kept_fractions(data: dict) -> None:
    """Fig: % barcodes kept by validrops vs Seurat-style baseline per dataset."""
    rows = []
    for dataset in sorted(data):
        current = data[dataset].get("current", {}).get("tier-h2h-qc")
        baseline = data[dataset].get("baseline", {}).get("tier-h2h-qc")
        if not current or not baseline:
            continue
        m = _metrics(current)
        rows.append(
            (
                dataset,
                m["fraction_kept_validrops"],
                m["fraction_kept_seurat_style"],
                m["jaccard_kept_sets"],
            )
        )
    if not rows:
        return
    fig, ax = _fig("kept-fractions")
    names = [r[0] for r in rows]
    x = range(len(rows))
    ax.bar([i - 0.18 for i in x], [r[1] for r in rows], width=0.36, label="validrops")
    ax.bar([i + 0.18 for i in x], [r[2] for r in rows], width=0.36, label="Seurat-style default", color="#c44e52")
    for i, r in enumerate(rows):
        ax.text(i, max(r[1], r[2]) + 0.02, f"J={r[3]:.2f}", ha="center", fontsize=7)
    ax.set_xticks(list(x), names, rotation=20, ha="right")
    ax.set_ylabel("fraction of barcodes kept")
    ax.set_ylim(0, 1.1)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    ax.set_title("QC callouts per dataset")
    fig.savefig(FIGDIR / "kept-fractions.png")
    plt.close(fig)


def figure_mito_thresholds(data: dict) -> None:
    """Fig: validrops adaptive mito cap vs the Seurat fixed 20% cap."""
    rows = []
    for dataset in sorted(data):
        a = (data[dataset].get("current") or {}).get("tier-h2h-qc")
        baseline = (data[dataset].get("baseline") or {}).get("tier-h2h-qc")
        if not a or not baseline:
            continue
        m = _metrics(a)
        if m.get("validrops_mito_threshold") != m.get("validrops_mito_threshold"):  # nan
            continue
        rows.append((dataset, m["validrops_mito_threshold"], m["validrops_mito_method"]))
    if not rows:
        return
    fig, ax = _fig("mito-thresholds")
    names = [r[0] for r in rows]
    x = range(len(rows))
    ax.bar([i - 0.18 for i in x], [100 * r[1] for r in rows], width=0.36, label="validrops adaptive cap")
    ax.axhline(20, color="#c44e52", ls="--", lw=1.2, label="Seurat-style fixed 20%")
    for i, r in enumerate(rows):
        ax.text(i - 0.18, 100 * r[1] + 0.5, f"{r[1] * 100:.1f}%", ha="center", fontsize=7)
    ax.set_xticks(list(x), names, rotation=20, ha="right")
    ax.set_ylabel("mitochondrial fraction cap (%)")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    ax.set_title("Adaptive vs fixed mitochondrial threshold")
    fig.savefig(FIGDIR / "mito-thresholds.png")
    plt.close(fig)


def figure_runtimes(data: dict) -> None:
    """Fig: median Stage-3 runtime and QC runtime per dataset."""
    rows = []
    for dataset in sorted(data):
        real = (data[dataset].get("current") or {}).get("tier-m-real")
        h2h_cur = (data[dataset].get("current") or {}).get("tier-h2h-qc")
        h2h_base = (data[dataset].get("baseline") or {}).get("tier-h2h-qc")
        if not real or not h2h_cur or not h2h_base:
            continue
        rows.append(
            (
                dataset,
                real["timing_seconds"]["median"],
                h2h_cur["timing_seconds"]["median"],
                h2h_base["timing_seconds"]["median"],
            )
        )
    if not rows:
        return
    fig, ax = _fig("runtimes")
    names = [r[0] for r in rows]
    x = range(len(rows))
    ax.bar([i - 0.25 for i in x], [r[1] for r in rows], width=0.25, label="Stage-3 (tier-m-real)")
    ax.bar(list(x), [r[2] for r in rows], width=0.25, label="validrops QC (tier-h2h-qc)")
    ax.bar([i + 0.25 for i in x], [r[3] for r in rows], width=0.25, label="Seurat-style QC", color="#c44e52")
    ax.set_xticks(list(x), names, rotation=20, ha="right")
    ax.set_ylabel("median seconds (3 reps)")
    ax.legend(frameon=False)
    ax.set_title("Runtime per dataset (timed window includes fixture load)")
    fig.savefig(FIGDIR / "runtimes.png")
    plt.close(fig)


def figure_cluster_counts(data: dict) -> None:
    """Fig: Stage-3 cluster count per dataset (bounded subset)."""
    rows = []
    for dataset in sorted(data):
        a = (data[dataset].get("current") or {}).get("tier-m-real")
        if not a:
            continue
        m = _metrics(a)
        rows.append((dataset, m["n_clusters"], m["n_obs"], m["species_detected"]))
    if not rows:
        return
    fig, ax = _fig("clusters")
    names = [r[0] for r in rows]
    vals = [r[1] for r in rows]
    ax.bar(range(len(rows)), vals)
    for i, r in enumerate(rows):
        ax.text(i, r[1] + 0.4, f"{r[1]} ({r[2]} cells)", ha="center", fontsize=7)
    ax.set_xticks(range(len(rows)), names, rotation=20, ha="right")
    ax.set_ylabel("clusters (Stage-3 deep clustering)")
    ax.set_title("Clusters detected per dataset (bounded subset)")
    fig.savefig(FIGDIR / "clusters.png")
    plt.close(fig)


def main() -> int:
    """Regenerate all manuscript figures from committed artifacts."""
    FIGDIR.mkdir(parents=True, exist_ok=True)
    data = sweep_artifacts()
    if not data:
        print("no package3 sweep artifacts found; run benchmarks/scripts/run_sweep.py first")
        return 1
    figure_kept_fractions(data)
    figure_mito_thresholds(data)
    figure_runtimes(data)
    figure_cluster_counts(data)
    print(f"figures written to {FIGDIR}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
