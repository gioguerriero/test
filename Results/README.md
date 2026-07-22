# `Results/` — Run Outputs

## Purpose

This folder collects the outputs produced by the optimization pipeline. Every time [`main.m`](../main.m) is executed, a **new dedicated sub-folder is created** and all the results of that run are saved into it. Nothing is overwritten: successive runs accumulate side by side, so past results are always preserved and reproducible.

## Naming convention

Each run folder is named:

```
<CometName>_run<N>
```

where:

- `<CometName>` is the sanitised name of the target comet (`selected_comet.name` with dashes and spaces removed — see `safe_name` in `main.m`), e.g. `Synthetic1`, `C2023X1`.
- `<N>` is an auto-incrementing integer. On each run, `main.m` scans the existing `<CometName>_run*` folders and picks the next free number, so run indices are unique per comet.

Example: `Results/Synthetic1_run1/`, `Results/Synthetic1_run2/`, `Results/C2023X1_run1/`.

## Contents of a run folder

Each run folder contains the following files, written by `main.m` at the corresponding stage of the pipeline:

| File | Produced after | Contents |
|---|---|---|
| `search_params.mat` | CR3BP global search setup | The `search_params` struct — the exact search configuration used (sweep bounds `vinf_ub_vec` / `q_vec` / `maximum_dv_vec`, manifold resolution, matching thresholds, filter factors, maneuver spacing). Guarantees the run is reproducible. |
| `global_results.mat` | CR3BP global search | The `global_results` cell array — every converged Halo→Moon→Comet transfer found by the sweep, each stored as a solution struct (optimal variables, matched manifold data, TOF and ΔV components). |
| `pareto_data.mat` | Pareto analysis | The `pareto_data` struct returned by `compute_pareto_front`: per-solution effective ΔV and TOF, plus the non-dominated front indices. |
| `pareto_front.fig` | Pareto analysis | The ΔV–TOF Pareto-front figure (MATLAB `.fig`, re-openable and editable). |
| `pareto_front.mat` | Pareto post-processing | The `pareto_front` struct — the front solutions repackaged and sorted by increasing TOF, with their indices, TOF/ΔV vectors, and full solution structs. |
| `refine_params.mat` | Ephemeris refinement setup | The `refine_params` struct — the refinement configuration (ΔV budget, minimum flyby altitude, node counts per phase, method flags). |
| `refined_front.mat` | Ephemeris refinement | The `refined_front` struct array — the ephemeris-refined trajectories for the selected front solutions, with convergence flags and the refined ΔV/TOF values. |

## Data lineage

The files reflect the order in which data is produced along the pipeline:

```
search_params ─▶ global_results ─▶ pareto_data ─▶ pareto_front ─▶ refine_params ─▶ refined_front
                                        │
                                        └─▶ pareto_front.fig
```

This mirrors the execution flow in `main.m` and in [`../Opt Manager/`](../Opt%20Manager/): the global search produces `global_results`, `compute_pareto_front` distils the Pareto set, and `run_refinement` produces the refined trajectories.

## Reloading results

Any file can be reloaded independently for post-processing without re-running the search. For example:

```matlab
load('Results/Synthetic1_run1/global_results.mat');   % -> global_results
load('Results/Synthetic1_run1/pareto_front.mat');     % -> pareto_front
```

The saved structs are also the expected inputs of the diagnostic tools in [`../Auxiliar/Checks and tests/`](../Auxiliar/) (e.g. `vinf_escape_study`, `earth_altitude_check`), which operate on a saved result set.

## Related folders

- [`../Opt Manager/`](../Opt%20Manager/) — the orchestration layer that generates every object saved here.
- Root [`../README.md`](../README.md) — the project-level overview and full data-flow description.
