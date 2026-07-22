# `Halo2Moon/` — Halo-to-Moon (H2M) Matching & Refinement

## Table of Contents

- [Overview](#overview)
- [Optimization Pipeline](#optimization-pipeline)
- [Design Variables](#design-variables)
- [Constraints](#constraints)
- [Cost Function](#cost-function)
- [Folder Structure](#folder-structure)

---

## Overview

This module implements the **Halo-to-Moon (H2M)** phase: the second half of the CR3BP trajectory search. Given the lunar-encounter condition fixed by the [`../Moon2Comet/`](../Moon2Comet/) phase, H2M finds a **dynamically natural departure** from the L2 halo orbit that reaches that encounter, and reconciles the residual velocity mismatch with the first deep-space manoeuvre, **DSM1**.

The core idea is **manifold matching**: the Earth-bound unstable manifold of the halo orbit already sweeps a corridor through the Earth–Moon region, so rather than blindly optimising a departure, the module *searches a precomputed database* of manifold trajectories for those that geometrically connect to the required lunar encounter. This turns an expensive optimisation into a fast nearest-neighbour query, made tractable by a **KD-tree**.

This folder also stores the precomputed halo data it depends on:

- `S_halo.mat` — states of the parking L2 halo orbit (synodic, non-dimensional).
- `unstable_dir.mat` — the unstable eigenvector at each halo point (the direction of the manifold departure).

---

## Optimization Pipeline

### 1. Initialization — build the manifold database (once)

`build_halo_manifold_db` perturbs a discretised set of halo departure points along the unstable eigendirection (with a small injection `ΔV_inj ≈ 15 m/s`) and propagates them **forward** in the CR3BP, keeping only trajectories that clear a minimum Earth altitude. The resulting dense cloud of manifold states (~10⁶–10⁷ points) is organised into a **KD-tree** (`KDTreeSearcher`), indexed by trajectory and time. This database is built **once** and reused across every optimisation run (cached to `.mat` by the orchestrator).

### 2. Trajectory generation — admissible pre-flyby states

For a given M2C lunar encounter, `build_moon_manifold_db` builds the *target side* of the match:

- The reference `v_inf` at the Moon is rotated over the **deflection cone** of half-aperture `δ_max` (the maximum bending achievable at the minimum flyby altitude), sampled on a polar grid of uniform areal density (~`M = 2000` directions).
- Each rotated `v_inf` is added to the Moon's velocity and **back-propagated** in the CR3BP, producing a database of admissible pre-flyby trajectories.

### 3. Matching — connect the two databases

`find_manifold_matches` queries the back-propagated Moon trajectories against the halo KD-tree with a **range search** (threshold `limit_dist_km`). Each spatial hit is a candidate Halo→Moon connection; the local velocity difference between the halo state and the Moon state is the delta-v that **DSM1** must supply. To avoid near-duplicate solutions, a **greedy diversity filter** keeps at most `K_per_moon` matches per Moon trajectory that differ sufficiently in TOF and phase (`div_tof_days`, `div_theta`).

> The subsequent **ΔV-feasibility** and additional diversity filtering are applied by the orchestrator [`cr3bp_search`](../Opt%20Manager/cr3bp_search.m), which discards matches whose implied DSM1 exceeds the remaining budget before refinement.

### 4. Optimization — local refinement

Each surviving match seeds a local `fmincon` optimisation, `refinement_halo2moon`, which turns the discrete match into a **continuous Halo-to-Moon trajectory**. `state_finder` locates the exact halo departure point for a given phase angle. The refinement **minimises the total time of flight** while enforcing arc closure, budget, deflection, and Earth-altitude constraints (see below).

### 5. Outputs

`refinement_halo2moon` returns the optimised design vector (`refined_traj`), the total H2M time of flight (`tof_halo2moon`), and the solver `exitflag`. Converged solutions are collected by `cr3bp_search` into the global result set.

---

## Design Variables

The refinement operates on a **5-element** design vector, defined in `refinement_halo2moon` and unpacked identically in `NC_refinement_halo2moon`:

| Var | Symbol | Physical meaning | Optimization role | Related functions |
|---|---|---|---|---|
| `x(1)` | `moon_fpa` | In-plane angle of the (H2M) flyby `v_inf` rotation | Steers the incoming lunar geometry | `vinf_rotation` |
| `x(2)` | `moon_out_of_plane` | Out-of-plane angle of the flyby `v_inf` rotation | Out-of-plane lunar geometry | `vinf_rotation` |
| `x(3)` | `halo_theta` | Departure phase angle on the halo orbit | Selects the manifold trajectory (departure point) | `state_finder` |
| `x(4)` | `tof_moon` | Time of flight on the Moon-side arc (DSM1 → flyby) | Length of the back-propagated lunar arc | `CR3BP` propagation |
| `x(5)` | `tof_halo` | Time of flight on the halo-side arc (injection → DSM1) | Length of the forward manifold arc | `CR3BP` propagation |

Two quantities are **not** free variables here:

- The **injection** is fixed at `ΔV_inj ≈ 15 m/s` along the unstable eigendirection (a modelling choice, not optimised).
- **DSM1** is *implicit*: it is the velocity discontinuity between the halo-side and Moon-side arcs at the matching point, i.e. an outcome of the design variables rather than a direct variable. Its magnitude is bounded through the constraints.

---

## Constraints

Implemented in `NC_refinement_halo2moon`, returning `[c, ceq]`:

| Constraint | Type | Physical interpretation | Implementation |
|---|---|---|---|
| **Arc position match** | `ceq` | The forward halo arc and the back-propagated Moon arc must meet in space | `ceq = final_state_halo(1:3) - final_state_moon(1:3)` |
| **ΔV budget** | `c` | The DSM1 velocity mismatch must stay within the remaining budget | `c = ‖Δv‖ - dv_max/Vstar ≤ 0` |
| **Maximum deflection angle** | `c` | The flyby cannot bend `v_inf` more than a patched-conic swingby allows at the minimum altitude | `c = cos(δ_max) - cos(θ) ≤ 0` |
| **Minimum Earth distance** | `c` | The trajectory must not pass below a minimum Earth altitude | `c = (h_min_earth + rEarth)/Lstar - d_min ≤ 0` |

Mathematically, the equality constraint enforces **arc closure** (the two-point boundary-value condition of the Halo→Moon connection); the inequalities encode the **physical feasibility** of the swingby (`δ_max` from Eq. of the maximum bending angle) and mission safety (Earth altitude). The `dv_max` on the right-hand side is the *remaining* budget `ΔV_max − ΔV_DSM2 − ΔV_inj`, so DSM1 competes with the manoeuvres already committed by M2C.

---

## Cost Function

Implemented in `OF_refinement_moon2comet`:

```matlab
tof = tof_moon + tof_halo;   % = x(4) + x(5)
```

The refinement **minimises the total Halo-to-Moon time of flight** (`x(4) + x(5)`). Propellant is deliberately **not** in the objective: instead the total delta-v is enforced as a **hard constraint** bounded by the remaining budget (see the ΔV-budget constraint above). This separation reflects the design intent — at this stage the ΔV envelope is already set by M2C and the mission budget, so the free objective to improve is transfer duration.

The global ΔV–ToF Pareto front is not built here; it emerges from the outer parametric sweep in `cr3bp_search`, which repeats the whole M2C + H2M process across the budget and weight grid.

---

## Folder Structure

| File | Role |
|---|---|
| `build_halo_manifold_db.m` | Propagates the halo unstable manifold from many departure points; produces the manifold point cloud used to build the KD-tree. |
| `build_moon_manifold_db.m` | Back-propagates admissible pre-flyby trajectories over the `δ_max` deflection cone (~`M = 2000` directions). |
| `find_manifold_matches.m` | KD-tree range search matching Moon trajectories to halo points; greedy diversity filter; returns candidate connections + DSM1 mismatch. |
| `refinement_halo2moon.m` | **Entry point** for the local refinement: `fmincon` minimising total TOF under the closure/budget/deflection/altitude constraints. |
| `NC_refinement_halo2moon.m` | Nonlinear constraints for the refinement (arc match, ΔV, deflection, Earth altitude). |
| `OF_refinement_moon2comet.m` | Objective for the refinement: total H2M time of flight. |
| `state_finder.m` | Returns the halo-orbit state (and index) closest to a given phase angle. |
| `S_halo.mat`, `unstable_dir.mat` | Precomputed halo states and unstable eigenvectors (inputs to the manifold generation). |

## Related folders

- [`../Moon2Comet/`](../Moon2Comet/) — the M2C phase that fixes the lunar encounter targeted here.
- [`../Opt Manager/`](../Opt%20Manager/) — `cr3bp_search`, which builds the KD-tree, applies the ΔV/diversity filters, and drives the outer sweep.
- [`../Auxiliar/`](../Auxiliar/) — `CR3BP`, `vinf_rotation`, `moon_state` used throughout.
