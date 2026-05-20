# Dynamic Transpilation

Transpiling the **dynamic** QPA circuit — a single circuit with classical
feedback via `if_test` control flow — onto IBM hardware coupling maps. This
folder also holds the circuit engine and QASM3 utilities that the other
experiment campaigns import.

## Library code

| File | Contents |
|---|---|
| `qpa_engine.py` | The QPA circuit engine. `HybridStrategy` builds one dynamic circuit with `if_test` classical feedback; `UnrolledStrategy` builds a set of static circuits, one per execution path; `HybridStrategy_n5_ntrials_1` is a flat-`if` variant. Also `QPARegisters`, the Pauli-twirling noise strategy, and the core ops (`apply_schur_test`, `apply_cyclic_rotation`). |
| `utils.py` | QASM3 manipulation and transpilation helpers. `strip_if_keep_else` / `strip_if_keep_then` convert a dynamic circuit into a static one by removing `if/else` control flow; `transpile_dynamic_best_of_n_unrolled` and friends transpile the dynamic circuit using a layout extracted from the unrolled longest-path circuit; `get_static_depth` / `get_static_gate_counts` measure a chosen branch. |

> Both files are imported by `experiments/end_to_end_hardware/end_to_end_dynamic_pittsburgh.ipynb`
> (via a `sys.path` insert) — keep them here.

## Notebooks

| Notebook | Purpose |
|---|---|
| `transpile_qpa_n3.ipynb` | Transpile the dynamic QPA circuit (`HybridStrategy`, n=3, k=2, n_trials=3) using the **unrolled-layout** strategy: build a static longest-path circuit, transpile *it* to get a good `initial_layout`, then re-transpile the dynamic circuit with that layout. Baseline vs. best-of-many-seeds. |
| `transpile_qpa_n5.ipynb` | Same procedure at n=5. |
| `transpile_qpa_n7.ipynb` | Same procedure at n=7. |
| `swapnet_switch_demo.ipynb` | Visual demo of the swap-network / cyclic-rotation construction used inside `HybridStrategy` for `T = 1, 2, 3`. |

## Running

These notebooks import `qpa_engine` / `utils` as siblings, so open them from
**inside this folder** (the default in Jupyter/VSCode). They do not import the
`core/` package and do not need the repo-root bootstrap.

## History

An earlier flat-`if` transpilation notebook (`transpile_qpa_flat_if_5.ipynb`)
was removed during the May 2026 reorganization — it was superseded by the
unrolled-layout approach used in `transpile_qpa_n5.ipynb`.
