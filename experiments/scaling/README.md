# Scaling — purity-amplification resource estimate

How the cyclic QPA circuit of [arXiv:2409.18167](https://arxiv.org/abs/2409.18167) scales under
transpilation, written to answer two reviewer questions about a NISQ-era (2026) demonstration:

1. **How do depth and gate count grow with the number of qubits in `ρ` and with `m`, and how fast
   do we reach 5k–15k gates?**
2. **Are the circuits known to the point of transpilation — can `F` and `P` (Fig. 3) be written in
   1- and 2-qubit gates?**

## Terminology (this repo ↔ the paper)

| Paper | OQPA knob | Meaning |
|---|---|---|
| qubits in `ρ` | `k` | qubits per copy of `ρ` (qudit dimension `d = 2^k`) |
| amplification order `m` | `N` (`n_registers`, odd) **and** `T` (`n_trials`) | `N` copies combined over `T` Schur-test + cyclic-rotation rounds; the minimal `m=2` demo is one swap test, i.e. `N=3, T=1` |

`F → Hadamard`, `controlled-P → controlled-SWAP (Fredkin, 7 two-qubit gates)`,
`P → SWAP (3 two-qubit gates)`. Intrinsic two-qubit cost: `N_2Q ≈ 3.5·k·T·(N−1)`.

## Notebook

| Notebook | Purpose |
|---|---|
| `transpile_scaling_purity_amplification.ipynb` | Self-contained (no `qiskit_aer` import). Rebuilds the longest-path QPA circuit (mirrors `UnrolledStrategy.build_longest_path` in `../dynamic_transpilation/qpa_engine.py`), decomposes `F`/`P`, and transpiles `k`/`N`/`T` sweeps against `FakeBrisbane` (127-qubit heavy-hex, realistic routing) and the `ibm_boston` basis with all-to-all connectivity (early-FTQC lower bound). |

Run with the **`opqa`** kernel (`~/.venvs/opqa`, qiskit 2.3.x). Both transpiler backends load
offline — no IBM Quantum account or `qiskit_aer` needed.

## Outputs

Plots are written to [`shared_data/experiment_results/scaling/`](../../shared_data/experiment_results/scaling/):
`scaling_vs_k.png`, `scaling_vs_N.png`, `scaling_vs_T.png`.
