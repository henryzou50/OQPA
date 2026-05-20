# Unrolled Scaling

How the **unrolled (static) QPA circuit** — the longest-path / all-Schur-tests-pass
branch — scales with the register count `N`. These notebooks build circuits and
transpile them; none of them execute on hardware or submit jobs.

Common parameters: `K = 2` (qubits per data register), `N = 3, 5, …, 23`,
`NO_RESET = True` (fresh ancillas per trial, so no `reset` instruction).

## Notebooks

| Notebook | Purpose |
|---|---|
| `qubit_scaling.ipynb` | Qubit count of the unrolled circuit vs `N`, for trial depths `T = 1, 2, 3`. |
| `transpile_scaling_unrolled_T1.ipynb` | Post-transpilation **2Q depth** and **2Q gate count** vs `N` at `T=1`. Transpiles each circuit with `N_SEEDS = 100` seeds at `optimization_level=3` against `ibm_boston` and `ibm_miami`, keeping the best. |
| `transpile_scaling_unrolled_T2.ipynb` | Same scan at `T = 2`. |
| `transpile_scaling_unrolled_T3.ipynb` | Same scan at `T = 3`. |
| `transpile_scaling_unrolled_overview.ipynb` | Combined view of the `T = 1, 2, 3` results. Reads the per-`T` JSON files produced by the notebooks above — run those first. |
| `cost_analysis.ipynb` | Circuit and Qiskit Runtime **job count** the unrolled strategy submits per λ point as a function of `N` (`paths(N, T=1) = 2^((N-1)/2)`). Validates the formula against `core.CircuitFactory`. |

## Run order

`transpile_scaling_unrolled_T{1,2,3}` are independent and write their own data;
run any subset, then `transpile_scaling_unrolled_overview` to combine them.
`qubit_scaling` and `cost_analysis` are standalone.

## Outputs

Written to [`shared_data/experiment_results/unrolled_scaling/`](../../shared_data/experiment_results/unrolled_scaling/):

```
unrolled_scaling/
├── T1/ , T2/ , T3/        # scaling_results_T*.{csv,json}, scaling_plot_T*.png,
│                          #   and per-N pre/post circuit diagrams
├── qubit_scaling.{csv,json,png}
└── scaling_overview_best.png
```
