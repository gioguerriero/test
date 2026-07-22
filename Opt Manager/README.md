# `Opt Manager/` — Optimization Orchestration Layer

## Table of Contents

- [Overview](#overview)
- [Architecture](#architecture)
- [Main Functions](#main-functions)
- [Data Flow](#data-flow)
- [Execution Logic](#execution-logic)
- [Related Folders](#related-folders)

---

## Overview

`Opt Manager/` is the **orchestration layer** of the framework. It contains no dynamics of its own: instead, it *sequences* the phase-specific optimizers and post-processing steps into the complete pipeline, manages the data that flows between them, and drives the parametric sweep that produces the Pareto front.

Everything the top-level entry point [`main.m`](../main.m) does — beyond setting constants and loading data — is delegated to three functions in this folder:

1. `cr3bp_search` — the complete CR3BP trajectory search (Moon-to-Comet + Halo-to-Moon, swept parametrically).
2. `compute_pareto_front` — extraction of the non-dominated ΔV–ToF front.
3. `run_refinement` — the ephemeris-refinement fidelity ladder.

Keeping orchestration separate from physics is a deliberate design choice: the dynamical model lives once in [`../Auxiliar/`](../Auxiliar/), each optimization phase owns its subproblem, and this layer only decides **what runs, in what order, on which data**.

---

## Architecture

### `cr3bp_search` — the search driver

This is the heart of the framework. It combines the two coupled optimization phases and wraps them in the parametric sweep.

- **Responsibility:** build (or load) the halo manifold **KD-tree** once, then for every combination of sweep parameters run the M2C phase, generate the Moon manifold database, match it against the halo tree, filter the matches, refine each survivor, and collect the converged transfers.
- **Interactions:**
  - calls [`optimization_moon2comet`](../Moon2Comet/optimization_moon2comet.m) (M2C phase);
  - calls [`build_moon_manifold_db`](../Halo2Moon/build_moon_manifold_db.m), [`find_manifold_matches`](../Halo2Moon/find_manifold_matches.m), and [`refinement_halo2moon`](../Halo2Moon/refinement_halo2moon.m) (H2M phase);
  - uses [`build_halo_manifold_db`](../Halo2Moon/build_halo_manifold_db.m) + `KDTreeSearcher` to build the manifold index.
- **Execution order:** M2C **first** (fixes the lunar encounter), then H2M matching/refinement — repeated across the triple loop.
- **Why it matters:** it implements the backward-solving strategy and the sweep that turns single transfers into a full trade-off study. It is the sole producer of `global_results`.

### `compute_pareto_front` — front extraction

- **Responsibility:** re-propagate every converged solution to compute its *effective* total ΔV (crucially, the DSM1 delta-v is the velocity mismatch at the Halo/Moon rendezvous, which must be measured rather than read off a design variable), then extract the **non-dominated set** in the (ΔV, ToF) plane and plot it.
- **Interactions:** consumes `global_results`; uses `CR3BP`, `state_finder`, and `vinf_rotation` (from `../Auxiliar/` and `../Halo2Moon/`) for the re-propagation.
- **Why it matters:** it converts a heterogeneous pile of solutions into the framework's primary deliverable — the ΔV–ToF Pareto front.

### `run_refinement` — refinement driver

- **Responsibility:** take one selected front solution and climb the fidelity ladder — reconstruct its CR3BP arcs, refine it in the Sun+Earth ephemeris model (multiple shooting), then in the full-ephemeris model with Moon gravity and TCMs.
- **Interactions:** calls the whole of [`../CR3BP 2 Ephemeris Sun-Earth/`](../CR3BP%202%20Ephemeris%20Sun-Earth/) (`build_initial_guess_MS`, `optimization_ephe`, `variables_organizer`, `check_altitude`) and [`../Full_ephemeris_conversion/`](../Full_ephemeris_conversion/) (`bplane_from_vinf`, `bplane_tcm`, the pre/post-flyby refinements, `variables_organizer_refined`).
- **Execution order:** Step 1 (Sun+Earth) → B-plane targeting → Step 2 post-flyby → Step 2 pre-flyby → assembly.
- **Why it matters:** it validates that CR3BP solutions survive in realistic dynamics, producing the refined trajectory objects.

### `plot_trajectory` — thin plotting wrapper

- **Responsibility:** a thin wrapper around `plot_full_trajectory` for the final visualisation.
- **Why it matters:** keeps the plotting entry point stable and decoupled from the refinement internals.

---

## Main Functions

### `cr3bp_search(c, selected_comet, S_halo, unstable_dir, sp)`

| | |
|---|---|
| **Purpose** | Run the full CR3BP search (M2C + H2M) across the `(v_inf^UB, q, ΔV_max)` sweep. |
| **Inputs** | `c` (constants), `selected_comet` (position, epoch, synodic geometry), `S_halo`, `unstable_dir`, `sp` (all search parameters: sweep vectors, manifold resolution, matching thresholds, filter factors, spacing). |
| **Outputs** | `global_results` — cell array of solution structs (optimal M2C/H2M variables, matched manifold data, TOF, ΔV components). |
| **Dependencies** | `Moon2Comet/`, `Halo2Moon/`, `KDTreeSearcher`, `parpool` (optional). |
| **Role** | Sole producer of the raw solution set; implements the parametric sweep. |

### `compute_pareto_front(global_results, c, S_halo, unstable_dir, eps_vel_ms)`

| | |
|---|---|
| **Purpose** | Compute the effective ΔV/ToF of each solution and extract the non-dominated front. |
| **Inputs** | `global_results`, constants, halo data, injection perturbation `eps_vel_ms`. |
| **Outputs** | `pareto_data` (per-solution ΔV/ToF + front indices); optionally the front figure handle. |
| **Dependencies** | `CR3BP`, `state_finder`, `vinf_rotation`. |
| **Role** | Turns solutions into the primary ΔV–ToF trade-off deliverable. |

### `run_refinement(c, best, selected_comet, S_halo, unstable_dir, rp)`

| | |
|---|---|
| **Purpose** | Refine a selected CR3BP solution through the Sun+Earth and full-ephemeris models. |
| **Inputs** | `c`, `best` (one front solution), `selected_comet`, `S_halo`, `unstable_dir`, `rp` (refinement parameters: budget, min flyby altitude, node counts, method flags). |
| **Outputs** | `out` — the fully refined trajectory struct, with pre/post convergence flags. |
| **Dependencies** | `CR3BP 2 Ephemeris Sun-Earth/`, `Full_ephemeris_conversion/`, `Auxiliar/`. |
| **Role** | Fidelity ladder; validates solutions in realistic dynamics. |

### `plot_trajectory(out, S_halo, c, comet_name)`

| | |
|---|---|
| **Purpose** | Plot the final full-ephemeris trajectory. |
| **Inputs** | Refined `out` struct, `S_halo`, constants, comet name. |
| **Outputs** | A figure (no return value). |
| **Dependencies** | `Full_ephemeris_conversion/plot_full_trajectory`. |

---

## Data Flow

The manager transforms a small set of well-defined objects:

```
sp (search_params) ─▶ cr3bp_search ─▶ global_results
                                          │
                                          ▼
                                 compute_pareto_front ─▶ pareto_data ─▶ (pareto_front, sorted by ToF)
                                                                             │
                                                                    select rank / best
                                                                             ▼
                              rp (refine_params) ────────────▶ run_refinement ─▶ out ─▶ plot_trajectory
```

- **`global_results`** is the boundary object between the *search* and the *analysis*: it carries every converged transfer with enough information to reconstruct and re-propagate it.
- **`pareto_data` / `pareto_front`** is the boundary object between *analysis* and *refinement*: the caller selects a front solution (`best`) to hand to `run_refinement`.
- **`out`** is the boundary object between *refinement* and *visualisation/output*: a single struct holding all arcs, epochs, and manoeuvres in synodic + J2000 form.

Each of these objects is also persisted to disk by `main.m` into the run's [`../Results/`](../Results/) folder, so any stage can be resumed or post-processed independently.

---

## Execution Logic

When `main.m` runs, control passes through this layer in the following order:

1. **Setup (in `main.m`).** Constants, SPICE kernels, and halo data are loaded; the target comet is selected and enriched with the Earth encounter geometry; a `Results/<comet>_runN/` folder is created.
2. **`cr3bp_search`.** The KD-tree is built or loaded. For each `(v_inf^UB, q)` pair: M2C is optimised, the Moon manifold database is generated and matched against the halo tree, matches are filtered by ΔV feasibility and diversity, and each survivor is refined by `refinement_halo2moon`. For each `ΔV_max` the top `N_save` transfers (by total TOF) are stored. → `global_results`.
3. **`compute_pareto_front`.** Every solution is re-propagated to obtain its effective ΔV; the non-dominated front is extracted and plotted. → `pareto_data`, and (in `main.m`) the `pareto_front` struct sorted by TOF.
4. **`run_refinement`** (looped over the selected front solutions). Each solution is refined through the Sun+Earth and full-ephemeris models. → `out` (per solution), collected into `refined_front`.
5. **Output & plotting.** `main.m` saves all objects to the run folder and calls `plot_full_trajectory` (directly or via `plot_trajectory`) on the chosen solution.

---

## Related Folders

- [`../Moon2Comet/`](../Moon2Comet/) and [`../Halo2Moon/`](../Halo2Moon/) — the two optimization phases invoked by `cr3bp_search`.
- [`../CR3BP 2 Ephemeris Sun-Earth/`](../CR3BP%202%20Ephemeris%20Sun-Earth/) and [`../Full_ephemeris_conversion/`](../Full_ephemeris_conversion/) — the refinement stages invoked by `run_refinement`.
- [`../Auxiliar/`](../Auxiliar/) — shared dynamics, frame conversions, and plotting.
- [`../Results/`](../Results/) — where the produced objects are saved.
- Root [`../README.md`](../README.md) — the project-level architecture and workflow overview.
