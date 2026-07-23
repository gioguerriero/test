# First-Time Tutorial — Running and Testing the Code

This guide walks a new user through a **first, lightweight run** of the pipeline, end to end: from installing the dependencies to producing a refined trajectory. It is meant to be followed **section by section** in [`main.m`](main.m), so you can see what each stage does before moving on.

The screenshots below come from a real run with the default settings shown here; your own plots should look very similar.

---

## 0. Before you start — heavy files, dependencies and MICE

1. **Download the heavy files.** The SPICE kernels (`kernels/`) and the cached manifold/KD-tree databases (`halo_manifold.mat`, `halo_tree.mat`) are too large for GitHub and are **not included in this repository**. Download them from:

   **<https://drive.google.com/drive/folders/1UOs_jA5jqdTVACn53cZaNN0Eka7E5gfF?usp=share_link>**

   Place `kernels/` and the two `.mat` files in the **project root**, next to `main.m`.

2. **Install the required packages.** Install everything listed in the **[§14 "Getting Started" → Requirements](README.md#14-getting-started)** section of the main README (MATLAB toolboxes and the MICE/SPICE toolkit).

3. **Point the code to your MICE install.** Open [`startup.m`](startup.m) and edit the `addpath` lines so they point to the **absolute path where you saved MICE** on your machine. The code loads SPICE from that address at startup, so if the path is wrong every `cspice_*` call will fail.

Once the heavy files are in place, the packages are installed, and `startup.m` points at your MICE, you can go straight to the code.

---

## 1. How to run — section by section

Rather than running the whole `main.m` at once, **execute it one section at a time** (each `%%` block is a MATLAB section). This way you understand exactly what the code does at every step and can inspect the intermediate results.

### 1.1 Keep the first run light

In the **`USER-CONFIGURABLE PARAMETERS`** block, start with a **single value** for each sweep vector, so the global search stays fast and cheap ([`main.m` lines 52–54](main.m)):

```matlab
vinf_ub_vec    = [0.7];                 % vinf upper bound(s) [km/s]
q_vec          = [0.442857142857143];   % q parameter
maximum_dv_vec = [0.85];                % total DV budget(s) [km/s]
```

(The commented-out `linspace(...)` lines just below show how you would sweep a full grid later, once you are comfortable — that is what produces a well-populated Pareto front.)

### 1.2 Leave the rest of the parameters alone (for now)

Do **not** touch the **`Moon manifold`**, **`Matching`**, **`Match filtering`** and **`Ephemeris refinement`** blocks on your first run. You can experiment with them later, once you are familiar with the code — although these are the exact values used for **all** the reference simulations.

The **only exception** is `days_per_node` (see below), which you will tune during the refinement step.

### 1.3 About `days_per_node`

`days_per_node` ([`main.m` line 86](main.m)) sets the multiple-shooting node density used by the **MATLAB refinement**. You cannot pick a good value up front: you first need a converged trajectory to look at, then adjust it. So leave it at the default for now — we come back to it in **Step 7 (Stage 4a)**.

---

## 2. Choose a target comet

The section after the parameters lists **several comets** used in the reference studies. Pick one by assigning it to `selected_comet` ([`main.m` line 161](main.m)). You can start with **`C2019U6`**, which is already pre-selected:

```matlab
selected_comet = C2019U6;   % comet processed by the pipeline
```

Available comet variables include `C2023X1`, `C2013US10`, `C2001Q4`, `C2008A1`, `C2013K1`, `C2015P3`, `C2016VZ18`, `C2024Y1`, `C2019U6`, `Synthetic1`–`Synthetic3`, `C2025N1` (3I/ATLAS) and `C2025R2`.

---

## 3. Results folder, halo data and the target comet

Running the next sections will:

- create the **results folder** for this run (`Results/<CometName>_runX/`, see the [`Results/` README](Results/README.md));
- load the precomputed **halo matrices** (`S_halo.mat`, `unstable_dir.mat`);
- plot the **position of the target comet in the synodic plane**.

You should get a plot like this:

![Target comet in the synodic plane](Images_tutorial/comet_synodic.jpg)
*Position of the target comet in the non-dimensional synodic (Sun–Earth rotating) frame.*

---

## 4. Global search

Leave the **global-search parameters** (`search_params`) at their defaults for now and launch the **`CR3BP GLOBAL SEARCH`** section ([`main.m` line 241](main.m)).

The search returns `global_results` — a set of **optimized Halo → Moon → Comet transfers** for the initial `vinf_ub_vec` / `q_vec` / `maximum_dv_vec` values you provided. With a single value in each, expect only a handful of solutions.

---

## 5. Pareto front

The next section (**`PARETO FRONT`**) filters the solutions and extracts the non-dominated front, then plots it:

![Example Pareto front](Images_tutorial/pareto_front_example.jpg)
*Pareto front (ΔV vs ToF).*

Note this is **not** a nicely populated front — that is expected, because we are running with only a few sweep values. Running the full `linspace(...)` grid from Step 1.1 is what fills it out.

---

## 6. Plot a chosen Pareto-front solution (optional)

In the **`Plot a chosen Pareto-front solution (CR3BP)`** section you can pick which front solution to visualize via `pareto_rank` ([`main.m` line 294](main.m)); `1` is the lowest-ToF solution. You can also **skip this section** entirely.

With the default `pareto_rank = 1` you should get:

![CR3BP trajectory of Pareto solution 1](Images_tutorial/plot_solution1_pareto_front_example.jpg)
*CR3BP trajectory of the first (lowest-ToF) Pareto-front solution.*

---

## 7. Ephemeris refinement

Now we reach the **`EPHEMERIS REFINEMENT`** section, which refines the selected CR3BP solution(s) in a real-ephemeris model.

### 7.1 TCM1 timing

Optionally, choose how many **days before the flyby** to place TCM1 via `refine_params.t_tcm1` ([`main.m` line 321](main.m)). Placing it **earlier** generally makes it **cheaper** — but **do not place it before DSM1**.

### 7.2 Which solutions to refine

By default the loop refines **only the first** solution ([`main.m` line 346](main.m)):

```matlab
% for ii = 1:N_front      % <-- uncomment to iterate over the whole front
for ii = [1]              % refine only the first solution
```

To refine **all** Pareto-front solutions, comment the second line and **uncomment the line above it** (`for ii = 1:N_front`).

### 7.3 Stage 4a — tune the nodes (multiple shooting)

`run_refinement` now starts. For a good result you should first **check the multiple-shooting node subdivision** of **Stage 4a** (see the [main README, §7 Step 4a](README.md#7-detailed-description-of-each-stage)) and decide whether `days_per_node` needs tuning.

The node-layout plot is produced at **[line 252 of `Opt Manager/run_refinement.m`](Opt%20Manager/run_refinement.m)**. Put a **breakpoint** there so execution pauses on the figure *before* the (potentially long) optimization runs; inspect the discretization, adjust `days_per_node`, and re-run the refinement section.

With the initial `days_per_node = [5 5 20 60]` you get **too many nodes**:

![Too many nodes](Images_tutorial/ms_wrong_nodes.jpg)
*`days_per_node = [5 5 20 60]` — the trajectory is over-discretized (too many nodes).*

Change it to `days_per_node = [20 20 40 60]` and re-run the refinement section to get a well-balanced subdivision:

![Correct nodes](Images_tutorial/ms_right_nodes.jpg)
*`days_per_node = [20 20 40 60]` — a sensible node layout.*

> **Tip:** also check the nodes on the arcs **near the Earth at departure** — zoom into the departure region to make sure that part is properly resolved too.

### 7.4 Stage 4a converges

With the nodes set, `run_refinement` proceeds. Depending on the node count it can take a while (up to ~1 hour). The converged Stage 4a trajectory looks like:

![Stage 4a converged trajectory](Images_tutorial/refined_4a.jpg)
*Converged Stage 4a (Sun + Earth ephemeris, multiple shooting) trajectory.*

At this point the trajectory data is also **exported to text files** ([`run_refinement.m` lines 289–290](Opt%20Manager/run_refinement.m)): `python_inputs.txt` and `python_inputs_serot.txt`. These feed the external **GODOT** refiner — an alternative to Stage 4b that is **not currently functional**; see the [`GODOT Refinement/` folder and its README](GODOT%20Refinement/README.md).

### 7.5 Stage 4b — full ephemeris with TCMs (the delicate part)

The code then moves to **Stage 4b** (see the [main README, §7 Step 4b](README.md#7-detailed-description-of-each-stage)): the **full-ephemeris** refinement that adds the **two TCMs**. This method is **not very robust** — it depends on the trajectory geometry and may sometimes fail to converge.

> **TCM2 timing.** Inside `run_refinement`, at **[line 329](Opt%20Manager/run_refinement.m)**, you can change how many days **after the Moon flyby** TCM2 is placed:
> ```matlab
> t_tcm2 = 15;   % [days] post flyby
> ```
> A **larger** value converges more easily, but the correction maneuver is **more expensive**: the post-flyby error accumulates over time, so the later you correct it, the worse it gets.

---

## 8. CR3BP vs refined Pareto front

For this comet, Stage 4b converges. The result is a **refined front** (a single solution, if you refined only one). In the **`PLOT — CR3BP vs Refined Pareto front`** section you can overlay it with the CR3BP front from the global search, to compare:

![CR3BP vs refined front](Images_tutorial/pareto_front_refined.jpg)
*CR3BP Pareto front vs the ephemeris-refined solution.*

As mentioned, Stage 4b is not robust: here the ΔV cost rises by about **140 m/s** — a lot, arguably too much — while the ToF stays essentially the same.

---

## 9. Plot a refined solution

Finally, the last section lets you plot a chosen **refined** solution. For our case (the first one):

![Refined solution](Images_tutorial/refined_solution_example.jpg)
*Refined full-ephemeris trajectory of the first solution.*

Zoom into the **initial part** of the trajectory: if you zoom in far enough on the flyby, you will see that the **lunar flyby is no longer instantaneous** but **natural**, following the real dynamics.

---

## Where to go next

- Run the **full sweep** (the `linspace(...)` values in Step 1.1) to get a properly populated Pareto front.
- Refine **all** front solutions (Step 7.2) instead of just the first.
- Experiment with the manifold / matching / filtering parameters once you are comfortable.
- Read the [main README](README.md) for the full scientific and architectural description, and each folder's own README for the details of that stage.
