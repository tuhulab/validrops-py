# tier-m-real / {'label': 'current', 'source_hashes': {'benchmarks/_datasets.py': '500d54e5eb650f36', 'benchmarks/_fixtures.py': '314d5c2621efc977', 'benchmarks/_measure.py': '6b66366b407d04f1', 'benchmarks/run.py': '67fbcdfcd2d72d31', 'src/validrops/_constants.py': 'd5b85f7e71247b1d', 'src/validrops/tl/_deviance.py': 'b3bb4d9afa657104', 'src/validrops/tl/_snn.py': 'c0536038a02c3c1b', 'src/validrops/tl/expression_metrics.py': '74dd2b785e2a94f7', 'src/validrops/tl/quality_metrics.py': 'e0331993952771b3'}, 'harness_hashes': {'benchmarks/_datasets.py': '500d54e5eb650f36', 'benchmarks/_fixtures.py': '314d5c2621efc977', 'benchmarks/_measure.py': '6b66366b407d04f1', 'benchmarks/run.py': '67fbcdfcd2d72d31'}}

- revision: a905519f0beac866a88a7977bb350edca9ed2fff (dirty=True)
- profile: S
- input fingerprint: cc8ae08f9cce0fc5
- parameters: `{"profile": "S", "nfeats": 300, "npcs": 10, "k_min": 150, "res_shallow": 0.1, "top_n": 10, "n_cells_cap": 1500, "dataset": "mouse-olfactory-bulb", "workload": "Stage-3 mouse-olfactory-bulb", "seed": 0}`
- timing seconds: median=28.375134, min=28.336712, mad=0.038422, reps=3, CV=0.0177
- peak RSS bytes: worker_max=1720778752, unit=bytes
- correctness: pass
- environment: {"os": "Linux", "os_release": "7.0.0-31-generic", "architecture": "x86_64", "python_version": "3.12.14", "python_implementation": "CPython", "logical_cpus": 8, "physical_cpus": 4, "available_memory_bytes": 23450210304, "numpy_version": "2.4.4", "scipy_version": "1.17.1", "sklearn_version": "unavailable", "anndata_version": "0.12.10", "igraph_version": "1.0.0", "blas_config": {}, "blas_threadpool":
