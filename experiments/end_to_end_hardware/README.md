# End-to-End Hardware

The full QPA fidelity-decay pipeline on **real IBM Quantum hardware**: search for
good transpilations, run the λ noise sweep, analyze QPU runtime cost, and plot
the combined results.

Common parameters: `K = 2`, `T = 1`. Each λ point is run as an ensemble of
Pauli-twirled circuit instances (`N_RANDOM`), batched into Qiskit Runtime jobs.

## Notebooks

### Transpilation search — run first
| Notebook | Purpose |
|---|---|
| `find_best_transpilation_unrolled_boston_t1.ipynb` | For each `N ∈ {3,5,7,9,11}` and each unrolled path, transpile many seeds against `ibm_boston` (opt level 3), keep the seed with lowest 2Q depth, and cache the result (`.qpy` + diagram + per-N `summary.json`). |
| `find_best_transpilation_unrolled_pittsburgh_t1.ipynb` | Same search for `N ∈ {3,5,7}` against `ibm_pittsburgh`. |

### λ sweeps
| Notebook | Purpose |
|---|---|
| `end_to_end_unrolled_boston_t1.ipynb` | Unrolled-strategy fidelity decay on `ibm_boston`, `N ∈ {3,5,7,9,11}`. Loads the cached transpilations above; λ sweep split into even/odd passes; checkpoints to CSV+JSON. |
| `end_to_end_unrolled_pittsburgh_t1.ipynb` | Same on `ibm_pittsburgh`, `N ∈ {3,5,7}`, adding job tags, a dynamical-decoupling on/off sweep, and a larger batch size. |
| `end_to_end_unrolled_pittsburgh_t1_show_3.ipynb` | Pittsburgh run / presentation variant focused on `N = 3`. |
| `end_to_end_unrolled_pittsburgh_t1_show_5.ipynb` | Pittsburgh run / presentation variant focused on `N = 5`. |
| `end_to_end_dynamic_pittsburgh.ipynb` | Fidelity decay using **dynamic** circuits (`if_test`) on `ibm_pittsburgh`. Imports `HybridStrategy` / `utils` from [`../dynamic_transpilation/`](../dynamic_transpilation/). |

### Analysis
| Notebook | Purpose |
|---|---|
| `runtime_analysis_pittsburgh_t1.ipynb` | QPU-usage breakdown of the Pittsburgh jobs, scoped to `N = 3`. Reads the job ids from the checkpoint CSVs. |
| `runtime_analysis_pittsburgh_t1_5.ipynb` | Same QPU-usage analysis scoped to `N = 5`. |
| `plot_all_results.ipynb` | Discovers and plots every fidelity-decay CSV under `shared_data/experiment_results/end_to_end_hardware/results/`. |

## Run order

For a given backend: `find_best_transpilation_*` → `end_to_end_unrolled_*`
(or `end_to_end_dynamic_*`) → `runtime_analysis_*` → `plot_all_results`. The
λ-sweep notebooks raise a clear error if the transpilation cache is missing.

## Outputs

Written to [`shared_data/experiment_results/end_to_end_hardware/`](../../shared_data/experiment_results/end_to_end_hardware/):

```
end_to_end_hardware/
├── results/              # one folder per experiment: fidelity CSV/JSON,
│                         #   decay plots, runtime CSV/plots, checkpoints
└── transpilations/       # per-N cached transpiled circuits (.qpy),
                          #   pre/post diagrams, summary.json
```

## Credentials

These notebooks call `QiskitRuntimeService()` and submit to IBM Quantum. Supply
credentials via a saved account / environment — never commit tokens or CRNs into
a notebook cell.
