# `CR3BP 2 Ephemeris Sun-Earth/` — CR3BP → Sun+Earth Ephemeris Refinement

## Table of Contents

- [Overview](#overview)
- [Optimization Pipeline](#optimization-pipeline)
- [Design Variables](#design-variables)
- [Constraints](#constraints)
- [Cost Function](#cost-function)
- [Folder Structure](#folder-structure)

---

## Overview

This module performs the **first fidelity upgrade** of the mission-design workflow: it takes a converged CR3BP trajectory (the output of the [`../Opt Manager/`](../Opt%20Manager/) search) and re-optimises it in a **real-ephemeris model** where the Earth and Moon follow their true SPICE positions, but only **Sun + Earth gravity** is integrated. The Moon is still treated as a zero-sphere-of-influence body, its effect reproduced by the instantaneous `v_inf` rotation exactly as in the CR3BP.

Its purpose is to demonstrate that CR3BP solutions **survive in a more realistic dynamical environment** while retaining approximately the same ΔV and time of flight. Because the whole trajectory (four arcs, three impulses, one flyby) is now integrated with time-varying real ephemerides, the solution is found with a **multiple-shooting** formulation, which is far more robust than a single forward shot for this stiff, sensitive problem.

This is the intermediate rung of the fidelity ladder: the CR3BP search feeds it, and its output feeds the full-ephemeris (Moon-gravity) refinement in [`../Full_ephemeris_conversion/`](../Full_ephemeris_conversion/).

---

## Optimization Pipeline

This module is driven by [`run_refinement`](../Opt%20Manager/run_refinement.m); the functions below implement each stage.

### 1. Initialization

`run_refinement` reconstructs the four CR3BP arcs of a selected solution and samples a set of intermediate nodes per phase. `build_initial_guess_MS` then assembles the **multiple-shooting initial guess**: it appends the node states of the four segments (`k_vec = [k1 k2 k3 k4]` nodes each) to the base design vector and returns the fixed node time-fractions `f_nodes` used to reconstruct node epochs.

### 2. Trajectory generation

Each shooting segment is propagated with `NBODY_J2000` (Sun + Earth ephemeris dynamics, in [`../Auxiliar/`](../Auxiliar/)). Node epochs are not free variables: they are reconstructed from the four event epochs and the fixed fractions `f_nodes` (`T = Ts + f·(Te−Ts)`), which keeps them automatically ordered.

### 3. Optimization

`optimization_ephe` runs `fmincon` (SQP) over the full multiple-shooting vector. It sets up:

- **linear inequality constraints** enforcing the chronological order of the four event epochs;
- **bounds** on the event epochs (±few days from the guess), on the halo-shift time `t_halo`, and on the position components of every control point (±`d_max`);
- the **nonlinear constraints** `nonlinear_constraints_ephe_multiple` (arc continuity, flyby, deflection, `v_inf` match, budget).

A custom `OutputFcn` stops the solver once the iterate is feasible and stagnant, avoiding wasted iterations chasing first-order optimality through adaptive integration noise.

> A single-shooting variant, `nonlinear_constraints_ephe`, exists for the `multiple_shooting = 0` path but the multiple-shooting formulation is the one used in practice.

### 4. Post-processing

`variables_organizer` assembles the optimised vector into the canonical **`out` trajectory struct** (all arcs, epochs, and manoeuvres in both synodic and J2000 form). `check_altitude` then verifies the minimum Earth altitude and the lunar flyby altitude of the refined solution. Several plotting helpers (`plot_from_state_vector_MS`, `plot_initial_guess_ephe`, `plot_synodic_ephe`) visualise the guess and the refined trajectory.

### 5. Outputs

The stage produces the `out_cr3bp` struct (via `variables_organizer`) and the optimised design vector `x_opt`, which become the input of the full-ephemeris refinement.

---

## Design Variables

The design vector has **17 base variables** plus **6·Σk_vec node states** (the multiple-shooting control points). The layout is defined in `optimization_ephe` and unpacked identically in `nonlinear_constraints_ephe_multiple` and `variables_organizer`.

| Var | Symbol | Physical meaning | Optimization role |
|---|---|---|---|
| `x(1:3)` | `ΔV_inj · 1e3` | Injection impulse [m/s, scaled] | Departure from the halo into the manifold |
| `x(4:6)` | `ΔV_DSM1 · 10` | First deep-space manoeuvre [km/s, scaled] | Deflection toward the lunar encounter |
| `x(7)` | `epoch_flyby` | Lunar flyby epoch, as *days before comet arrival* | Times the swingby |
| `x(8:10)` | `v_comet_arr` | Comet arrival velocity | Free end-state velocity for the backward legs |
| `x(11:13)` | `ΔV_DSM2 · 10` | Second deep-space manoeuvre [km/s, scaled] | Final phasing to the comet |
| `x(14)` | `t_halo` | Waiting/shift time along the halo orbit [days] | Adjusts the departure point/time |
| `x(15)` | `epoch_dep` | Departure epoch, *days before comet* | Times the departure |
| `x(16)` | `epoch_dsm1` | DSM1 epoch, *days before comet* | Times DSM1 |
| `x(17)` | `epoch_dsm2` | DSM2 epoch, *days before comet* | Times DSM2 |
| `x(18 …)` | node states | 6·k states per segment (synodic) | Multiple-shooting control points that break the trajectory into short, well-conditioned arcs |

**Design decisions worth noting.** The impulses carry scaling factors (`·1e3`, `·10`) so all variables are of order 1. Epochs are expressed as *days before comet arrival*, giving a monotonically decreasing chain (`dep > dsm1 > flyby > dsm2 > 0`) that the linear constraints exploit. Node **epochs** are *not* variables — only the four physical event epochs are — while node **states** are, which is the essence of the multiple-shooting parametrisation.

---

## Constraints

### Nonlinear constraints — `nonlinear_constraints_ephe_multiple`

| Constraint | Type | Physical interpretation | Implementation |
|---|---|---|---|
| **Node continuity** | `ceq` | Each shooting sub-arc must join the next node in position and velocity | 6D gap `S_prop(end) − S_node`, scaled by `Lscale`/`Vscale` |
| **Injection cap** | `c` | The injection impulse cannot exceed `max_dv_inj` | `‖ΔV_inj‖·1e3/max_dv_inj − 1 ≤ 0` |
| **Flyby position** | `ceq` | At the flyby epoch the spacecraft must coincide with the Moon | `S_flyby(1:3) − moon(1:3)`, scaled |
| **Maximum deflection** | `c` | The swingby cannot bend `v_inf` beyond the patched-conic limit | `cos(δ_max) − cos(θ) ≤ 0` |
| **`v_inf` magnitude match** | `ceq` | Incoming and outgoing `v_inf` must have equal magnitude (elastic flyby) | `‖v_inf,out‖ − ‖v_inf,in‖ = 0` |
| **Total ΔV budget** | `c` | The summed impulses must respect the mission budget | `Σ‖ΔV‖ − max_dv·1e-3 ≤ 0` |

The **node-continuity equalities are the multiple-shooting closure conditions**: they stitch the independently propagated sub-arcs into a single continuous trajectory. The scaling by `Lscale`/`Vscale` aligns the numerical tolerance (`ConstraintTolerance ≈ 1e-3`) with a physical tolerance of ~10 km / ~0.1 m/s.

### Linear constraints & bounds — `optimization_ephe`

- **Event ordering** (`A·x ≤ b`): enforces `dep > dsm1 > flyby > dsm2` with the minimum spacing `min_days_between`.
- **Bounds**: event epochs constrained to ±few days of the guess; `t_halo` to ±5 days; each control-point position component to ±`d_max` (velocities free).

---

## Cost Function

Implemented in `objective_function_ephe`:

```matlab
dv = ‖ΔV_inj‖/1e3 + ‖ΔV_DSM1‖/10 + ‖ΔV_DSM2‖/10;   % (multiple-shooting indices)
```

The objective is the **weighted sum of the three manoeuvre magnitudes** — the total propellant expenditure of the transfer, with each term divided by its scaling factor to recover physical km/s. Time of flight is *not* penalised here: it is essentially fixed by the CR3BP solution being refined, so the refinement's job is to close the trajectory in the higher-fidelity model at minimum additional ΔV. (An internal `tof` term is computed but set to zero, i.e. inactive.)

---

## Folder Structure

| File | Role |
|---|---|
| `optimization_ephe.m` | **Entry point.** Sets up the multiple-shooting problem (linear constraints, bounds, nonlinear constraints, custom stop) and runs `fmincon`. |
| `build_initial_guess_MS.m` | Builds the multiple-shooting initial guess: appends node states per segment and returns the fixed node time-fractions. |
| `nonlinear_constraints_ephe_multiple.m` | Nonlinear constraints for the multiple-shooting formulation (continuity, flyby, deflection, `v_inf` match, budget). |
| `nonlinear_constraints_ephe.m` | Single-shooting variant of the constraints (used only when `multiple_shooting = 0`). |
| `objective_function_ephe.m` | Cost: weighted total manoeuvre delta-v. |
| `variables_organizer.m` | Assembles the optimised vector into the canonical `out` trajectory struct (synodic + J2000, all events and manoeuvres). |
| `check_altitude.m` | Verifies minimum Earth altitude and lunar flyby altitude of a refined solution; optional synodic plot / Blender export. |
| `plot_from_state_vector_MS.m` | Visualises the multiple-shooting solution/guess (nodes + arcs) in the synodic frame. |
| `plot_initial_guess_ephe.m` | Plots the initial-guess trajectory in the Sun-centred ECLIPJ2000 frame. |
| `plot_synodic_ephe.m` | Plots the ephemeris-propagated trajectory in the synodic frame, optionally overlaying the optimised solution. |

## Related folders

- [`../Opt Manager/`](../Opt%20Manager/) — `run_refinement`, which drives this stage.
- [`../Full_ephemeris_conversion/`](../Full_ephemeris_conversion/) — the next fidelity step (adds Moon gravity + TCMs).
- [`../Auxiliar/`](../Auxiliar/) — `NBODY_J2000`, `CR3BP`, frame conversions, `vinf_rotation`.
