# Experiments

Research notebooks for the Quantum Purification Algorithm (QPA) cyclic project,
organized by experiment campaign. Each subfolder is self-contained and has its
own `README.md` with per-notebook detail.

| Folder | What it studies |
|---|---|
| [`unrolled_scaling/`](unrolled_scaling/) | How the **unrolled (static)** QPA circuit scales with the register count `N` — qubit count, transpiled 2Q depth / gate count, and circuit/job cost. No hardware execution. |
| [`dynamic_transpilation/`](dynamic_transpilation/) | Transpiling the **dynamic** QPA circuit (`if_test` control flow) for hardware. Holds the circuit engine (`qpa_engine.py`) and QASM3 utilities (`utils.py`) shared by the other campaigns. |
| [`end_to_end_hardware/`](end_to_end_hardware/) | Full fidelity-decay pipeline on **real IBM hardware** — transpilation search, λ sweeps, runtime/QPU-cost analysis, and combined plotting. |

## Running the notebooks

Every notebook in `unrolled_scaling/` and `end_to_end_hardware/` begins with a
**repo-root bootstrap cell** that walks up to the repository root (the folder
containing `core/`), `chdir`s there, and puts it on `sys.path`. As a result the
notebooks import `core` / `execution` / `analysis` and resolve every
`shared_data/...` path correctly **no matter which directory they are launched
from**.

The `dynamic_transpilation/` notebooks have no bootstrap by design — they import
`qpa_engine` / `utils` as siblings and should be opened from inside that folder
(the default when you open the notebook in Jupyter/VSCode).

## Where the data lives

Experiment outputs are **not** stored next to the notebooks. They live under
[`shared_data/experiment_results/`](../shared_data/experiment_results/) so the
results can be version-controlled and shared:

```
shared_data/experiment_results/
├── end_to_end_hardware/
│   ├── results/          # fidelity-decay CSV/JSON + plots, per experiment
│   └── transpilations/   # cached transpiled circuits (.qpy) + diagrams
└── unrolled_scaling/     # T1/T2/T3 scaling results + qubit-scaling data
```

## Hardware / credentials

The `end_to_end_hardware/` notebooks and the `dynamic_transpilation` scaling
work submit to IBM Quantum and need a configured `QiskitRuntimeService`
account. Provide credentials through the environment / saved account — **do not
hard-code tokens or CRNs in the notebooks.**
