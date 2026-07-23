# `GODOT Refinement/` — Full-Ephemeris Refinement (Python / ESA GODOT)

## Purpose

This folder contains a **Python** implementation of a full-ephemeris trajectory refinement built on top of **ESA's [GODOT](https://godot.io.esa.int/) library** (the Agency's operational-grade orbit/trajectory propagation and optimisation toolkit). Its role is to take a trajectory produced by the MATLAB optimisation pipeline in this repository and **re-optimise it in a high-fidelity, real-ephemeris dynamical model**, beyond the fidelity reached by the MATLAB-side ephemeris refinement (`Full_ephemeris_conversion/`).

Conceptually, this is the **final rung of the fidelity ladder**: the MATLAB framework performs the trajectory *search* and the CR3BP/intermediate-ephemeris optimisation; GODOT then re-propagates and refines the selected solution with full, operational-quality dynamics (real Sun–Earth–Moon ephemerides, multiple shooting, exact-gradient nonlinear optimisation).

> **Status — work in progress.** This code is **not currently functional**: in its present state it raises errors and does not complete a refinement. It is included deliberately because the **overall structure is correct** — the problem set-up, the import of the MATLAB interface file, and the intended GODOT workflow are all in place and follow a sound design. It was **left unfinished for lack of time**, not because the approach is wrong, and is meant to be **picked up and completed** by a future developer. See [Known limitations](#known-limitations--next-steps) below.

---

## Interface with the MATLAB pipeline

The two environments communicate through a single, human-readable **text file** exported on the MATLAB side by [`write_python_inputs`](../Auxiliar/write_python_inputs.m) (and its SEROT variant `write_python_inputs_serot`). That file carries everything needed to reconstruct the trajectory as an initial guess:

- the **Halo parameters** (`m`, `T_halo`);
- the key **mission states** in the synodic non-dimensional CR3BP frame (pre-injection, post-injection, post-DSM1, flyby periapsis, post-DSM2);
- the **durations** of each arc (injection→DSM1, DSM1→flyby, flyby→DSM2, DSM2→comet);
- the **comet encounter date**.

```
MATLAB pipeline ──(write_python_inputs → python_inputs.txt)──▶ GODOT Refinement (Python)
```

See [§8 of the project README](../README.md#8-full-ephemeris-refinement-bridge-godot) for the description of the bridge from the MATLAB side.

---

## The refinement process (end to end)

The driver is [`main.py`](main.py). It runs a **two-stage** GODOT optimisation, mirroring the physical structure of the mission:

1. **Read inputs** — `main.py` parses `python_inputs.txt` into a dictionary (states, durations, dates). All mission epochs are then derived by working **backwards** from the comet-encounter date.

2. **Build the physical universe** — GODOT loads `Universe/universe.yml`, which defines the Sun, Earth and Moon (with the DE432 ephemeris), the reference frames (ICRF and the rotating **SEROT** Sun–Earth frame) and the spacecraft. This is the real-ephemeris model that replaces the idealised CR3BP.

3. **Stage A — Halo orbit optimisation.** The L2 Halo parking orbit, initially known only as an approximate CR3BP solution, is discretised into control points over one orbit and **re-optimised with real ephemerides** using a **multiple-shooting** scheme: each arc is propagated forward and continuity (*match*) constraints force the arcs to meet, while station-keeping ΔV is minimised. The optimised Halo state at the injection epoch is then extracted.

4. **Stage B — Transfer optimisation.** The four transfer arcs (injection → DSM1 → Moon flyby → DSM2 → comet) are assembled into a single multiple-shooting timeline, seeded from the CR3BP states/durations. The dynamics are then switched from **Sun+Earth (`ES_gravity`)** to the **full Sun–Earth–Moon model (`EMS_gravity`)**, so the lunar flyby becomes a genuine gravity assist. A lunar-altitude inequality constraint keeps the trajectory clear of the Moon's surface. The problem minimises the total transfer ΔV (injection + deep-space manoeuvres).

5. **Diagnostics & plots.** Before the second optimisation, `main.py` prints the minimum lunar distance along the guessed arc and produces several SEROT-frame plots (top view, Earth–Moon zoom, flyby geometry) so the initial guess can be inspected.

6. **Output.** The pre-optimisation and optimised trajectories are written to YAML: `traj_R2_preopt.yml` and `traj_R2.yml`.

The optimisation engine throughout is **IPOPT** (via **PyGMO**), driven with **exact gradients** obtained from GODOT's automatic differentiation.

---

## Files in this folder

### Python source

| File | Role |
|------|------|
| [`main.py`](main.py) | **Main driver.** Runs the full two-stage refinement (Sections 1–12): input parsing, universe set-up, Halo optimisation, hand-off, transfer set-up, `ES_gravity`→`EMS_gravity` switch, diagnostics/plots, and YAML export. |
| [`transfer_ephe.py`](transfer_ephe.py) | **Problem builders.** Defines `config_halo()` and `config_trajectory()`, which turn CR3BP states/durations into GODOT timelines (control points, manoeuvres, matches) and the corresponding optimisation problems. Also holds the shared physical constants (`D_SUN_EARTH`, `GM_SUN`, `GM_EARTH`, `mu_se`, `T_SUN_EARTH`) and the CR3BP↔physical conversions. |
| [`aux_fun.py`](aux_fun.py) | **Timeline-element helpers.** Small builders `ctr()`, `man()`, `match()` that assemble the verbose GODOT dictionaries for control points, impulsive manoeuvres and match constraints. |
| [`totalDV.py`](totalDV.py) | **Objectives & constraints.** `ScalarTimeEvaluable` classes: `totalDV_traj` / `totalDV_traj_geo` (transfer ΔV), `totalDV_halo` (Halo station-keeping ΔV), and `flybyMoonAltitude` (lunar-altitude inequality at the flyby node). |
| [`halo_ephe.py`](halo_ephe.py) | **Standalone study script** (not imported by `main.py`). Optimises **only** the Halo orbit over several revolutions and plots the result. Useful in isolation to understand and debug Stage A. |

### Data & configuration

| File / folder | Role |
|---------------|------|
| [`python_inputs.txt`](python_inputs.txt) | The **interface file** from MATLAB (synodic non-dimensional CR3BP frame). This is the input `main.py` reads. |
| [`python_inputs_serot.txt`](python_inputs_serot.txt) | Alternative interface file in the **SEROT** frame (produced by `write_python_inputs_serot`), for the SEROT variant of the workflow. |
| `Universe/` | GODOT physical model: `universe.yml` (bodies, frames, gravitational constants, spacecraft), `de432s.bsp` (JPL DE432 planetary ephemeris) and `gm_de431.tpc` (gravitational constants kernel). |
| `traj_R2_preopt.yml` | **Output** — transfer trajectory *before* the second optimisation (the initial guess), written by `main.py`. |
| `traj_R2.yml` | **Output** — transfer trajectory *after* the full-ephemeris optimisation. |
| `pre_opt.mat`, `traiettoria_cr3bp.mat` | Supporting MATLAB data (CR3BP trajectory / pre-optimisation state) kept alongside the Python code for reference. |
| `__pycache__/` | Python's compiled-bytecode cache (auto-generated; not documented here). |

---

## Installing GODOT

The GODOT library is **not** shipped with this repository — you must download and install it yourself. It is distributed by ESA, together with its installation instructions and documentation, at:

- **https://godot.io.esa.int/docs/**

Follow ESA's instructions there to set up GODOT (typically inside a dedicated conda environment) together with **PyGMO/IPOPT**. The `main.py` driver assumes the environment is activated with `conda activate godot`.

## A short home-made guide to GODOT

Alongside this code I have also included a **short guide I wrote myself** to help understand and use GODOT. It is meant as a practical, hands-on companion to the official ESA documentation — a quicker way to get oriented with the concepts and workflow used by the scripts in this folder (universes, timelines, control points, manoeuvres, matches, and the optimisation set-up). New users are encouraged to read it before diving into `main.py`.

## Running it

Once GODOT is installed and `python_inputs.txt` is present in this folder, the intended entry point is:

```bash
conda activate godot
python main.py
```

---

## Known limitations & next steps

The refinement **does not currently run to completion** — it raises errors partway through. The structure, however, is sound, so a future developer can pick it up. Reasonable starting points:

- get the import/propagation set-up running against a **known-good exported trajectory** (e.g. the provided `traj_R2_preopt.yml`), before re-enabling the optimisation;
- verify that the **gravitational constants** used here match those of the MATLAB pipeline (see the open items noted at the top of `main.py`);
- move the **integration centre** (`coi`) closer to the dominant body along each phase (Earth up to the flyby, Moon at the flyby) to improve numerical conditioning;
- tune the **node spacing** per arc (`nodes_days_vec` in `main.py`) so the delicate lunar-approach/departure phases are resolved finely and the interplanetary cruise coarsely.

---

## Related folders

- [`../Auxiliar/`](../Auxiliar/) — `write_python_inputs` / `write_python_inputs_serot`, which produce the interface file consumed here.
- [`../Full_ephemeris_conversion/`](../Full_ephemeris_conversion/) — the MATLAB-side full-ephemeris refinement (with Moon gravity and TCMs) that this Python code is intended to supersede in fidelity.
- Root [`../README.md`](../README.md) — project overview; the GODOT bridge is described in §8.
