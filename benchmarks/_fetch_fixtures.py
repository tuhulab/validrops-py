"""Download the parotid-gland 10x fixtures (Zenodo 7470666, CC-BY-4.0).

Writes to tests/data/pbmc-parotid/ (git-ignored). Resumable: an existing
file of the expected size is kept; partial files are re-downloaded.
"""

import sys
import urllib.request
from pathlib import Path

BASE = "https://zenodo.org/api/records/7470666/files/{name}/content"
TARGETS = {
    "raw_all.h5": ("raw_feature_bc_matrix.h5", 74_632_310),
    "raw.h5": ("filtered_feature_bc_matrix.h5", 45_000_867),
}
OUT = Path(__file__).resolve().parents[1] / "tests" / "data" / "pbmc-parotid"


def export_mtx() -> None:
    """Export the filtered H5 as a 10x MTX directory for the R cross-check.

    R's Seurat reads the same counts through base ``ReadMtx`` (the conda R
    toolchain here cannot compile hdf5r, and the harness itself stays R-free).
    Output: ``tests/data/pbmc-parotid/mtx/`` with matrix.mtx, features.tsv
    (gene id + symbol columns), barcodes.tsv.
    """
    import warnings

    import scanpy as sc
    from scipy.io import mmwrite

    src = OUT / "raw.h5"
    if not src.exists():
        raise SystemExit(f"run the download first; missing {src}")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Variable names are not unique")
        adata = sc.read_10x_h5(src)
    adata.var_names_make_unique()
    mtx = OUT / "mtx"
    mtx.mkdir(parents=True, exist_ok=True)
    mmwrite(mtx / "matrix.mtx", adata.X.T)  # 10x convention: genes x cells
    features = adata.var["gene_ids"].astype(str).to_frame("gene_id")
    features.insert(1, "symbol", adata.var_names.astype(str))
    features.to_csv(mtx / "features.tsv", sep="\t", header=False, index=False)
    adata.obs_names.astype(str).to_series().to_csv(mtx / "barcodes.tsv", sep="\t", header=False, index=False)
    print(f"ok mtx export -> {mtx} ({adata.shape[0]} cells x {adata.shape[1]} genes)")


def download(url: str, dest: Path, expected: int) -> None:
    if dest.exists() and dest.stat().st_size == expected:
        print(f"keep {dest.name} ({dest.stat().st_size} bytes)")
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "validrops-bench/0.1"})
    with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as fh:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    if tmp.stat().st_size != expected:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"{dest.name}: size mismatch ({tmp.stat().st_size} != {expected})")
    tmp.rename(dest)
    print(f"ok {dest.name} ({dest.stat().st_size} bytes)")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    for local, (remote, size) in TARGETS.items():
        download(BASE.format(name=urllib.parse.quote(remote)), OUT / local, size)
    return 0


if __name__ == "__main__":
    import sys
    import urllib.parse

    if len(sys.argv) > 1 and sys.argv[1] == "--export-mtx":
        export_mtx()
        sys.exit(0)
    sys.exit(main())
