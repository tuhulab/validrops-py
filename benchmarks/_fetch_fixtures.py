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

# Additional cross-dataset fixtures (Package 3 sweep).
EXTRA_TARGETS = {
    "mouse-brain-5k": {
        "record": "https://zenodo.org/api/records/7993711",
        "files": [("adata.h5ad", 267_996_169)],
        "out": "adata.h5ad",
    },
    "carotid-visium": {
        "record": "https://zenodo.org/api/records/22024050",
        "files": [("filtered_feature_bc_matrix.h5", 1_114_377)],
        "out": "raw.h5",
    },
    "mouse-olfactory-bulb": {
        "record": "https://zenodo.org/api/records/15024747",
        "files": [("10X_V1 Mouse Olfactory Bulb.zip", 115_403_244)],
        "out": "source.zip",
        "extract": True,
    },
}
EXTRA_OUT = Path(__file__).resolve().parents[1] / "tests" / "data"


def _record_files(record_url: str) -> dict[str, str]:
    import json

    with urllib.request.urlopen(
        urllib.request.Request(record_url, headers={"User-Agent": "validrops-bench/0.1"}), timeout=120
    ) as r:
        record = json.load(r)
    return {f["key"]: f["links"]["self"] for f in record["files"]}


def download_extra() -> None:
    """Download the Package-3 sweep fixtures (mouse brain nuclei h5ad, Visium h5).

    Same contract as the parotid download: size-checked, re-runnable, skipped
    when a file is already present at the expected size.
    """
    for name, spec in EXTRA_TARGETS.items():
        dest_dir = EXTRA_OUT / name
        dest_dir.mkdir(parents=True, exist_ok=True)
        files = _record_files(spec["record"])
        for remote, size in spec["files"]:
            dest = dest_dir / spec["out"]
            if dest.exists() and dest.stat().st_size == size:
                print(f"keep {name}/{spec['out']} ({size} bytes)")
                continue
            if remote not in files:
                print(f"SKIP {name}: file {remote!r} not in record {spec['record']}")
                continue
            tmp = dest.with_suffix(dest.suffix + ".part")
            req = urllib.request.Request(files[remote], headers={"User-Agent": "validrops-bench/0.1"})
            with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            if tmp.stat().st_size != size:
                tmp.unlink(missing_ok=True)
                raise RuntimeError(f"{name}/{spec['out']}: size mismatch ({tmp.stat().st_size} != {size})")
            tmp.rename(dest)
            print(f"ok {name}/{spec['out']} ({dest.stat().st_size} bytes)")
            if spec.get("extract"):
                import zipfile

                with zipfile.ZipFile(dest) as zf:
                    zf.extractall(dest_dir / "extracted")
                    print(f"extracted {name}: {[i.filename for i in zf.infolist()][:12]}")


def build_parotid_raw_subset() -> None:
    """Build a bounded elbow-crossing subset of the unfiltered parotid matrix.

    The raw matrix is 2.2M barcodes (93M UMI); loading it per benchmark child
    costs ~15 s and >1 GB. For the sweep the registry uses a deterministic
    subset that crosses the barcode-rank elbow: barcode positions 3001-4500 in
    Cell Ranger's descending-UMI order contain a mix of real cells and ambient
    droplets — the operating regime of Stage-1 knee detection. Output is a
    standard 10x H5 next to the source (git-ignored).
    """
    import warnings

    import scanpy as sc

    src = OUT / "raw_all.h5"
    dest = OUT / "raw_subset.h5ad"
    if dest.exists():
        print(f"keep pbmc-parotid/{dest.name} ({dest.stat().st_size} bytes)")
        return
    if not src.exists():
        raise SystemExit(f"missing {src}; run the parotid download first")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Variable names are not unique")
        full = sc.read_10x_h5(src)
    full.var_names_make_unique()
    subset = full[3000:4500].copy()
    subset.write_h5ad(dest)
    print(f"ok pbmc-parotid/{dest.name} (barcodes 3001-4500, {subset.shape})")


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
    if len(sys.argv) > 1 and sys.argv[1] == "--all":
        main()
        download_extra()
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == "--extra":
        download_extra()
        build_parotid_raw_subset()
        sys.exit(0)
    sys.exit(main())
