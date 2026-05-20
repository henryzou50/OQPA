# Experiment Results

Output data for the experiment campaigns in [`experiments/`](../../experiments/).
This folder is version-controlled (unlike the old gitignored `data/`) so results
can be shared.

```
experiment_results/
├── end_to_end_hardware/        <- produced by experiments/end_to_end_hardware/
│   ├── results/                #   fidelity-decay CSV/JSON, plots, runtime data,
│   │                           #   per-(N,dd) checkpoints
│   │   ├── end_to_end_unrolled_boston_t1/
│   │   ├── end_to_end_unrolled_pittsburgh_t1/
│   │   └── end_to_end_dynamic_pittsburgh/
│   └── transpilations/         #   cached transpiled circuits (.qpy) + diagrams
│       ├── end_to_end_unrolled_boston_t1/
│       └── end_to_end_unrolled_pittsburgh_t1/
└── unrolled_scaling/           <- produced by experiments/unrolled_scaling/
    ├── T1/ , T2/ , T3/         #   transpilation-scaling results per trial depth
    ├── qubit_scaling.{csv,json,png}
    └── scaling_overview_best.png
```

Notebooks reference these paths relative to the repository root; their
repo-root bootstrap cell makes that work regardless of launch directory.

## Note on removed data

Five stale result folders from earlier `ibm_marrakesh` runs and a preliminary
`ibm_pittsburgh` "small" run were removed during the May 2026 reorganization
(`end_to_end`, `end_to_end_unrolled_t1`, `end_to_end_unrolled_t3`,
`end_to_end_dynamic_t1`, `end_to_end_unrolled_small`). They had no surviving
notebook that produced or consumed them.
