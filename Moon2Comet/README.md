# `Moon2Comet/` — Moon-to-Comet (M2C) Optimization

## Table of Contents

- [Overview](#overview)
- [Optimization Pipeline](#optimization-pipeline)
- [Design Variables](#design-variables)
- [Constraints](#constraints)
- [Cost Function](#cost-function)
- [Folder Structure](#folder-structure)

---

## Overview

This module implements the **Moon-to-Comet (M2C)** optimization: the first of the two coupled phases that make up the CR3BP trajectory search. Its job is to determine the **lunar encounter conditions** and the **post-flyby heliocentric leg** (DSM2 + coast) that deliver the spacecraft to the target comet at its prescribed rendezvous epoch.

M2C is optimised **before** the Halo-to-Moon phase ([`../Halo2Moon/`](../Halo2Moon/)) because the comet is a *rigid* boundary condition (a fixed position at a fixed epoch), whereas the halo departure offers a continuous family of natural manifold trajectories. Solving M2C first fixes the lunar encounter that the H2M phase must then reproduce.

The entry point is [`optimization_moon2comet`](optimization_moon2comet.m), invoked by [`cr3bp_search`](../Opt%20Manager/cr3bp_search.m) once per outer-loop parameter combination. The optimised lunar-encounter state (pre-flyby `v_inf` and lunar phase `θ_moon`) and `ΔV_DSM2` are then passed to the H2M phase as boundary conditions.

---

## Optimization Pipeline

The transfer modelled here starts at the **Moon** (parametrised lunar orbit), applies the flyby as an instantaneous `v_inf` rotation, coasts to the deep-space manoeuvre **DSM2**, and reaches the comet. The optimisation proceeds in **three sequential stages** for robustness on this nonconvex, multimodal problem.

### 1. Initialization

`optimization_moon2comet` receives the comet position (inertial and synodic), the encounter epoch, the `v_inf` bounds (`vinf_data`), the cost weight `q`, and the minimum maneuver-spacing vector `min_days_between`. From the spacing it derives the non-dimensional lower bounds on the DSM2 timing (`D_flyby_dsm2_adim`, `D_dsm2_comet_adim`).

### 2. Trajectory generation

For a candidate design vector, the trajectory is reconstructed as:

- `moon_state` returns the Moon's synodic state at phase angle `θ_moon`.
- `vinf_rotation` (in [`../Auxiliar/`](../Auxiliar/)) rotates the reference `v_inf` by the flight-path (`fpa`) and out-of-plane (`oop`) angles, producing the post-flyby velocity.
- The spacecraft is propagated in the CR3BP (`CR3BP`) from the Moon to the DSM2 point (a fraction `β` of the total TOF).
- The remaining leg to the comet is closed either by a Lambert arc (`lambert`, used in the GA stage via `get_DSM_info`) or directly by the DSM2 design variable (fmincon stages).

### 3. Optimization (three stages)

| Stage | Solver | Design vars | Objective function | Constraints |
|---|---|---|---|---|
| **1. Global** | `ga` (genetic algorithm) | 6 | `objective_function_moon2comet_ga` | none (bounds only) |
| **2. Local refine** | `fmincon` (SQP) | 9 | `objective_function_moon2comet_fmincon` | `nonlinear_constraints_moon2comet` |
| **3. Fixed-TOF** | `fmincon` (SQP) | 7 | `objective_function_moon2comet_tof` | `nonlinear_constraints_moon2comet_tof` |

- **Stage 1** globally explores the search space and provides an initial guess.
- **Stage 2** refines it with a gradient solver, now including the explicit DSM2 impulse as design variables (scaled to order 1 via `scaling1`, obtained from `get_DSM_info`).
- **Stage 3** fixes the **time of flight** (adjusted so the Moon phasing is physically consistent at the flyby epoch, using its real ephemeris position) and re-optimises the remaining variables, minimising the DSM2 delta-v alone.

### 4. Post-processing

`get_DSM_info` reconstructs the DSM2 delta-v vector (via Lambert) from a GA solution, used both to seed Stage 2 and to compute the scaling factor.

### 5. Outputs

`optimization_moon2comet` returns three structures: `result_ga`, `result_fmincon` (Stage 2), and `result_fmincon_modtof` (Stage 3). The Stage-3 `bestfeasible.x` is the one consumed downstream by `cr3bp_search` as the converged M2C solution.

---

## Design Variables

The full 9-variable vector (Stages 1–2; Stage 1 uses the first 6) is defined in `optimization_moon2comet` and unpacked identically in every objective/constraint function of this folder.

| Var | Symbol | Physical meaning | Optimization role | Defined / bounded in |
|---|---|---|---|---|
| `x(1)` | `Vinf` | Lunar hyperbolic excess-velocity magnitude | Sets the flyby energy regime | `optimization_moon2comet` (`lb/ub` from `vinf_data`) |
| `x(2)` | `θ_moon` | Lunar phase angle (position along the Moon's orbit); zero points toward Sun–Earth L2, increasing counter-clockwise | Locates the encounter along the lunar orbit | `optimization_moon2comet` |
| `x(3)` | `fpa` | Flight-path angle of the post-flyby `v_inf` | In-plane steering of the swingby | `optimization_moon2comet` |
| `x(4)` | `oop` | Out-of-plane angle of the post-flyby `v_inf` | Out-of-plane steering of the swingby | `optimization_moon2comet` |
| `x(5)` | `tof` | M2C time of flight (Moon → comet) | Controls transfer duration | `optimization_moon2comet` |
| `x(6)` | `β` | Split fraction locating DSM2 along the heliocentric leg (`β·tof` after the flyby) | Positions the deep-space manoeuvre | `optimization_moon2comet` |
| `x(7:9)` | `ΔV_DSM2` | Three components of the DSM2 impulse (scaled) | Final phasing correction to hit the comet | Added in Stage 2; seeded by `get_DSM_info` |

**Stage 3** reduces to 7 variables `[Vinf, fpa, oop, β, ΔV_DSM2(3)]`: `θ_moon` and `tof` are **fixed** (the TOF is corrected from the true Moon phasing), so they are removed from the optimised set and passed in as constants (`theta_fixed`, `new_tof`).

The DSM2 components are handled with scaling factors (`scaling1` in Stage 2, `scaling2` in Stage 3) so that all variables are of comparable magnitude — essential for well-conditioned gradient steps.

---

## Constraints

All constraints are implemented as MATLAB nonlinear-constraint functions returning `[c, ceq]` (inequalities `c ≤ 0`, equalities `ceq = 0`).

### `nonlinear_constraints_moon2comet` (Stage 2)

| Constraint | Type | Physical interpretation | Implementation |
|---|---|---|---|
| **Comet rendezvous** | `ceq` | The propagated trajectory must reach the comet's position at the encounter epoch | `ceq = arrival_state(1:3) - target_pos_synodic` |
| **Flyby→DSM2 spacing** | `c` | DSM2 must not be executed immediately after the flyby | `c = D_flyby_dsm2_adim - tof·β ≤ 0` |
| **DSM2→comet spacing** | `c` | DSM2 must leave enough time before the comet | `c = D_dsm2_comet_adim - tof·(1-β) ≤ 0` |

Mathematically, the equality constraint enforces the **two-point boundary-value condition** of the transfer; the two inequalities keep the manoeuvre timing physically meaningful (the spacing bounds come from `min_days_between`).

### `nonlinear_constraints_moon2comet_tof` (Stage 3)

With the TOF fixed, the position match is retained (`ceq = (arrival - target)·1e2`, scaled for conditioning). The maneuver spacing is instead enforced directly through the `β` bounds computed in `optimization_moon2comet` (`beta_lb`, `beta_ub`), so no explicit spacing inequality is needed. An optional delta-v cap is present but commented out.

---

## Cost Function

The M2C objective realises the **ΔV–time-of-flight trade-off** through the scalar weight `q ∈ [0,1]`, which is swept externally (in `cr3bp_search`) to generate different points of the Pareto front.

| Stage | Objective | Meaning |
|---|---|---|
| **1 (GA)** | `objective_function_moon2comet_ga` = `‖ΔV_DSM2‖ + q·tof` | DSM2 delta-v (from the Lambert velocity mismatch) plus a time penalty. |
| **2 (fmincon)** | `objective_function_moon2comet_fmincon` = `‖ΔV_DSM2‖·Vstar + q·tof` | DSM2 delta-v (dimensional) plus a time penalty. |
| **3 (fixed-TOF)** | `objective_function_moon2comet_tof` = `‖ΔV_DSM2‖·Vstar` | DSM2 delta-v only (TOF is now a constant). |

In all cases the **propellant term** is the magnitude of the DSM2 impulse — the only free manoeuvre on the M2C leg (the flyby is ballistic).

> **Discrepancy with `main.m` and the report.** Both `main.m` and the accompanying report describe the M2C cost as the *convex combination* `(1-q)·ΔV + q·ToF`. That is **not** what Stages 1 and 2 actually implement above: the code uses `ΔV + q·ToF`, with **no `(1-q)` factor on the delta-v term**.
>
> This is not a cosmetic difference — the two formulas trace the Pareto front differently as `q` is swept:
> - **As implemented (`ΔV + q·ToF`):** the delta-v term is **never discounted**; increasing `q` only adds more weight to the time-of-flight term on top of the full ΔV. Decreasing `q` toward 0 recovers a pure minimum-ΔV objective.
> - **As described in `main.m`/the report (`(1-q)·ΔV + q·ToF`):** at `q → 1` the delta-v term is driven to **zero weight**, so the optimiser would pursue a pure minimum-ToF solution regardless of propellant cost.
>
> With the formula actually implemented, empirical testing showed that **no feasible solutions were found for `q` above roughly 0.8** — because the ΔV term is never discounted away, past that point the added ToF penalty is not enough to pull the solver toward a feasible faster (and necessarily more expensive) trajectory within the mission's ΔV/spacing constraints.

---

## Folder Structure

| File | Role |
|---|---|
| `optimization_moon2comet.m` | **Entry point.** Runs the 3-stage GA → fmincon → fixed-TOF optimization; returns the three result structs. |
| `moon_state.m` | Analytic Moon state (position + velocity) in the synodic frame as a function of the phase angle `θ_moon`. |
| `objective_function_moon2comet_ga.m` | Stage-1 cost: DSM2 delta-v (via Lambert) + `q·tof`. |
| `objective_function_moon2comet_fmincon.m` | Stage-2 cost: DSM2 delta-v (dimensional) + `q·tof`. |
| `objective_function_moon2comet_tof.m` | Stage-3 cost: DSM2 delta-v only (fixed TOF). |
| `nonlinear_constraints_moon2comet.m` | Stage-2 constraints: comet position match + DSM2 spacing. |
| `nonlinear_constraints_moon2comet_tof.m` | Stage-3 constraints: comet position match (spacing enforced via `β` bounds). |
| `get_DSM_info.m` | Reconstructs the DSM2 delta-v vector from a GA solution via a Lambert arc; used to seed and scale Stage 2. |

## Related folders

- [`../Halo2Moon/`](../Halo2Moon/) — the H2M phase that consumes this module's lunar-encounter output.
- [`../Opt Manager/`](../Opt%20Manager/) — `cr3bp_search`, which calls this module inside the parametric sweep.
- [`../Auxiliar/`](../Auxiliar/) — `CR3BP`, `vinf_rotation`, `lambert`, `synodic2car`/`car2synodic` used throughout.
