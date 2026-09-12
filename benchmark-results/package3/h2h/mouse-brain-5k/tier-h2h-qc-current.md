# tier-h2h-qc / {'label': 'current', 'source_hashes': {'benchmarks/_datasets.py': '500d54e5eb650f36', 'benchmarks/_fixtures.py': '314d5c2621efc977', 'benchmarks/_measure.py': '6b66366b407d04f1', 'benchmarks/run.py': '67fbcdfcd2d72d31', 'src/validrops/_constants.py': 'd5b85f7e71247b1d', 'src/validrops/_pipeline.py': '46190a5c959e9860', 'src/validrops/pp/quality_filter.py': '33ad7dd7b6b9b95f', 'src/validrops/pp/rank_barcodes.py': '233fb13c244c283a', 'src/validrops/tl/_deviance.py': 'b3bb4d9afa657104', 'src/validrops/tl/_segmented.py': '188361e99803d228', 'src/validrops/tl/_snn.py': 'c0536038a02c3c1b', 'src/validrops/tl/_stats.py': '3c8e6e279f78e741', 'src/validrops/tl/_uik.py': '24f430dcd6bd663e', 'src/validrops/tl/expression_metrics.py': '74dd2b785e2a94f7', 'src/validrops/tl/quality_metrics.py': 'e0331993952771b3'}, 'harness_hashes': {'benchmarks/_datasets.py': '500d54e5eb650f36', 'benchmarks/_fixtures.py': '314d5c2621efc977', 'benchmarks/_measure.py': '6b66366b407d04f1', 'benchmarks/run.py': '67fbcdfcd2d72d31'}}

- revision: a905519f0beac866a88a7977bb350edca9ed2fff (dirty=True)
- profile: S
- input fingerprint: d62e267b222d61da
- parameters: `{"profile": "S", "n_cells_cap": 1500, "dataset": "mouse-brain-5k", "seurat_mito_cap": 0.2, "seurat_min_features": 200, "seurat_min_counts": 500, "seurat_min_cells": 3, "workload": "QC head-to-head mouse-brain-5k", "seed": 0}`
- timing seconds: median=23.613588, min=23.497295, mad=0.116293, reps=3, CV=0.0077
- peak RSS bytes: worker_max=1490763776, unit=bytes
- correctness: pass
- environment: {"os": "Linux", "os_release": "7.0.0-31-generic", "architecture": "x86_64", "python_version": "3.12.14", "python_implementation": "CPython", "logical_cpus": 8, "physical_cpus": 4, "available_memory_bytes": 23450210304, "numpy_version": "2.4.4", "scipy_version": "1.17.1", "sklearn_version": "unavailable", "anndata_version": "0.12.10", "igraph_version": "1.0.0", "blas_config": {}, "blas_threadpool":
