# =============================================================================
# HALO_EPHE.PY — Standalone Halo-orbit optimization in real ephemerides
# =============================================================================
#
# Standalone / study script (not imported by main.py). It refines ONLY the
# L2 Halo parking orbit — the same job as Section 7 of main.py — but over
# several orbits and with diagnostic plots.
#
# What it does:
#   1. Loads the physical Solar-System model (Sun, Earth, Moon).
#   2. Starts from an approximate Halo solution (simplified CR3BP model) and
#      "corrects" it using the real planetary positions.
#   3. Optimizes the orbit minimizing propellant (total station-keeping ΔV).
#   4. Plots the result (three panels) and saves it as an image.
#
# Glossary:
#   ΔV (delta-v)      = velocity change produced by the thrusters (propellant proxy)
#   CR3BP             = Circular Restricted Three-Body Problem (simplified Sun-Earth model)
#   Halo orbit        = 3D periodic orbit about the Sun-Earth L2 point
#   SEROT             = Sun-Earth rotating frame (XY plane rotates with the Earth)
#   Multiple shooting = split the trajectory into short arcs and impose that they
#                       meet at the junction (match) points
#   IPOPT             = nonlinear optimizer (Interior Point OPTimizer)
# =============================================================================


# =============================================================================
# SECTION 1 — IMPORTS
# =============================================================================

import numpy as np

from godot.core import num, tempo, astro, events, ipfwrap
# GODOT (ESA trajectory framework). Submodules used: tempo (epochs/time),
# astro (astrodynamic constants), events (event detection).

from godot.core import autodif as ad
# Automatic differentiation (gradients) — used by the optimizer.

from godot.core.autodif import bridge as br
import godot.core.util as c_util
from godot.model import interface, frames, common, prop
from godot import cosmos
from godot.cosmos import util

import matplotlib.pyplot as plt
from ruamel import yaml
import time, os, copy
import pygmo as pg
# PyGMO: provides the IPOPT interface and the "population" structure.

import pygmo_plugins_nonfree as ppnf
import midas
# CR3BP library (used in transfer_ephe.py to integrate the non-dimensional orbit).

from transfer_ephe import *
# Provides config_halo(), config_trajectory() and the physical constants
# D_SUN_EARTH, GM_SUN, GM_EARTH, mu_se, T_SUN_EARTH.

from totalDV import *
# Provides totalDV_halo / totalDV_traj (objective functions).


# =============================================================================
# SECTION 2 — UNIVERSE INITIALIZATION
# =============================================================================

c_util.suppressLogger()
# Silence GODOT's internal logging (keeps the console readable).

config_uni = util.load_yaml('Universe/universe.yml')
# universe.yml describes the bodies (Sun, Earth, Moon), ephemerides, reference
# frames, gravitational constants and the spacecraft.

universe = cosmos.Universe(config_uni)
# GODOT Universe object: the Solar-System model queried for body states and forces.


# =============================================================================
# SECTION 3 — REFERENCE HALO-ORBIT PARAMETERS (CR3BP)
# =============================================================================
# Initial point on the Halo orbit in the non-dimensional Sun-Earth CR3BP
# rotating frame: [x, y, z, vx, vy, vz]. Comes from an external CR3BP solver.

x0_cr3bp = [1.00737013, 0., -0.00292083, 0., 0.01298799, 0.]
# x0_cr3bp = [1.01098261377409, 0., 0.0042711451, 0., -0.0110683024, 0.] # own data
#   x  = 1.00737 : just beyond L2 (L2 is at ~1.01 AU from the Sun in adim units)
#   y  = 0       : in the ecliptic plane, y = 0 at t = 0
#   z  = -0.00292: small out-of-plane component (the orbit is 3D, "halo")
#   vx = 0       : zero x-velocity (symmetry condition)
#   vy = 0.01299 : in-plane velocity, needed to close the orbit
#   vz = 0       : zero z-velocity at t = 0

T_cr3bp = 3.0860603629916086
# T_cr3bp = 3.08302339296 # own data
# Halo orbit period in CR3BP non-dimensional time.
# Convert to seconds with: T_cr3bp * T_SUN_EARTH / (2*pi).


# =============================================================================
# SECTION 4 — MISSION TIME WINDOW
# =============================================================================

date_start = tempo.Epoch('2029-01-01T00:00:00.0 TDB')
# Start of the window in which the Halo orbit is sought. TDB = Barycentric
# Dynamical Time (the time scale used in astrodynamics).

date_end = tempo.Epoch('2029-12-01T00:00:00.0 TDB')
# End of the window (defines the mission context; not used directly here).


# =============================================================================
# SECTION 5 — MULTIPLE-SHOOTING DISCRETIZATION PARAMETERS
# =============================================================================

n_pt = 4
# Control points per orbit. The orbit is split into n_pt segments; each is
# forward-propagated and continuity (match) is imposed at each segment midpoint.

n_orb = 3
# Number of complete Halo orbits. Total segments = n_pt * n_orb = 12.
# More orbits = better correction but more optimization variables.


# =============================================================================
# SECTION 6 — OBJECTIVE FUNCTION: TOTAL HALO ΔV
# =============================================================================

dv_halo = totalDV_halo(n_pt * n_orb, universe)
# Sums the ΔV of all station-keeping maneuvers (man0 … man11 for 12 segments).

universe.evaluables.add('TotalDV_halo', dv_halo)
# Register the objective in the universe under the name 'TotalDV_halo' so the
# optimizer can evaluate and minimize it.


# =============================================================================
# SECTION 7 — SCALES AND PERTURBATIONS FOR THE OPTIMIZATION
# =============================================================================

dx = 100000
# Position bound half-width: 100 000 km — how far each control point may move
# from the CR3BP guess in x, y, z.

dv = 0.2
# Velocity bound half-width: 0.2 km/s.

scales_vect = [2000000, 870000, 600000, dv, dv, dv]
# Normalization scales for the optimizer variables. Positions x/y/z differ
# because the Halo orbit is asymmetric (more extended in X than in Z).

delta_vect = np.array([dx, dx, dx, dv, dv, dv])
# Bound half-widths around the CR3BP guess: [guess - delta, guess + delta].


# =============================================================================
# SECTION 8 — BUILD TRAJECTORY AND OPTIMIZATION PROBLEM
# =============================================================================

config_traj, config_prob = config_halo(
    x0_cr3bp,    # CR3BP initial state
    T_cr3bp,     # CR3BP period
    date_start,  # start epoch
    n_pt,        # control points per orbit
    n_orb,       # number of orbits
    scales_vect, # normalization scales
    delta_vect   # bound half-widths
)
# config_halo() (in transfer_ephe.py): integrates the CR3BP orbit at n_pt*n_orb
# epochs, converts each to physical km/km·s via the real ephemerides, and
# assembles the GODOT timeline (control points, maneuvers, matches) and the
# problem (min TotalDV_halo, free variables, scales, bounds).

traj = cosmos.Trajectory(universe, config_traj)
# GODOT Trajectory object (structure only; not yet propagated).

prob = cosmos.Problem(universe, [traj], config_prob, useGradient=True)
# GODOT optimization problem. useGradient=True enables automatic differentiation.


# =============================================================================
# SECTION 9 — FIGURE SETUP (INITIAL GUESS)
# =============================================================================

fig, ax = plt.subplots(1, 3)
# One figure, three side-by-side axes (ax[0], ax[1], ax[2]).

# --- Panel 0: XY plane (top view on the ecliptic) ---
ax[0].grid()
ax[0].set_aspect('equal')
ax[0].set_xlabel('$X_{rot}$ (km)')
ax[0].set_ylabel('$Y_{rot}$ (km)')

# --- Panel 1: XZ plane (side view) ---
ax[1].grid()
ax[1].set_aspect('equal')
ax[1].set_xlabel('$X_{rot}$ (km)')
ax[1].set_ylabel('$z_{rot}$ (km)')    # Z = out-of-plane component

# --- Panel 2: accumulated ΔV vs time ---
ax[2].grid()
ax[2].set_xlabel('$Y_{rot}$ (km)')    # NB: X axis here is actually time (MJD)
ax[2].set_ylabel('$z_{rot}$ (km)')    # this panel shows ΔV vs time


# =============================================================================
# SECTION 10 — INITIAL COMPUTE (BEFORE OPTIMIZATION)
# =============================================================================

traj.compute(partials=False)
# Propagate the trajectory at the current (CR3BP guess) values. partials=False:
# skip the partial derivatives (only needed for the optimization).


# =============================================================================
# SECTION 11 — OPTIMIZER SETUP (PyGMO + IPOPT)
# =============================================================================

problem = pg.problem(prob)
# Wrap the GODOT problem in a PyGMO problem (GODOT evaluates objective/constraints).

tol_con = 1e-5
# Constraint tolerance: a constraint is satisfied when its value is below this.

problem.c_tol = [tol_con] * problem.get_nc()
# Apply the tolerance to all constraints (get_nc() = number of constraints).


# =============================================================================
# SECTION 12 — INITIAL POPULATION
# =============================================================================

x0 = prob.get_x()
# Initial free-variable vector (positions, velocities, ΔV components, node timings).

pop = pg.population(problem, 0)
# Empty PyGMO population (IPOPT, gradient-based, needs a single individual).

pop.push_back(x0)
# Seed the population with the initial guess.


# =============================================================================
# SECTION 13 — RUN THE OPTIMIZATION
# =============================================================================

ip = pg.ipopt()
# IPOPT instance (nonlinear constrained optimizer).

ip.set_numeric_option("tol", 1e-3)
# IPOPT convergence tolerance (smaller = more accurate but slower).

algo = pg.algorithm(ip)
algo.set_verbosity(1)
# Print one line per IPOPT iteration (objective, constraint violation, ...).

pop = algo.evolve(pop)
# *** OPTIMIZATION HAPPENS HERE *** IPOPT evaluates objective/constraints,
# computes gradients via autodiff, and iterates until the tolerance is met.


# =============================================================================
# SECTION 14 — EXTRACT AND PRINT THE SOLUTION
# =============================================================================

traj.compute(partials=False)
# Recompute the trajectory with the optimized parameters.

traj_up = traj.applyParameterChanges()
# Collect the optimizer's parameter changes as an update dict.

conf_update = util.deep_update(config_traj, traj_up)
# Recursively update config_traj with the optimized values.

# --- Total accumulated ΔV (from the propagator) ---
dv_tot = universe.evaluables.get('SC_dv').eval(traj.point('end_point'))
# 'SC_dv' is the cumulative ΔV GODOT tracks along the trajectory; evaluated at
# the end point => total ΔV.

print(f'Total dv: {dv_tot * 1e3} m/s')   # GODOT works in km/s internally

t_tot = traj.point('end_point') - traj.point('ctr0')
# Total trajectory duration [s] (difference of the two epochs).

print(f'Total time: {t_tot / 86400} days')


# =============================================================================
# SECTION 15 — MANUAL ΔV RECOMPUTATION (independent check)
# =============================================================================
# Recompute the total ΔV by explicitly summing each maneuver vector, as an
# independent cross-check against 'SC_dv'.

sol = traj.getTimelineSolution()
# Timeline solution: all events (control points, maneuvers, matches) with their
# optimized epochs and states, as a list of lists.

dv_tot = 0

for li in sol:
    for t in li:
        if "man" in t.name and (not t.name.endswith("end") and not t.name.endswith("start")):
            # Keep only maneuver events (name contains "man"), excluding the
            # special "end"/"start" points. Kept: "man0", "man3", ...;
            # excluded: "man_injection_end", "match0_start".
            dvx = universe.evaluables.get(t.name + '_dv_x').eval(t.epoch)
            dvy = universe.evaluables.get(t.name + '_dv_y').eval(t.epoch)
            dvz = universe.evaluables.get(t.name + '_dv_z').eval(t.epoch)
            dv_tot += np.sqrt(dvx**2 + dvy**2 + dvz**2)

print(dv_tot)
# Should match Section 14 (small differences are numerical).


# =============================================================================
# SECTION 16 — SAMPLE THE TRAJECTORY FOR PLOTTING
# =============================================================================

E_grid = tempo.EpochRange(traj.point('ctr0'), traj.point('end_point')).createGrid(86400)
# Epoch grid from ctr0 to end_point, one point per day (86400 s).

XROT = []    # [epoch(MJD), x, y, z, vx, vy, vz] per sample
dv_g = []    # accumulated ΔV per sample
epog = []    # epochs in MJD (Modified Julian Date)

for E_ in E_grid:
    x_rot = universe.frames.vector6('Earth', 'SC_center', 'SEROT', E_)
    # SC state w.r.t. Earth in the SEROT frame at epoch E_ ([x y z vx vy vz]).
    XROT.append(np.concatenate(([E_.mjd()], x_rot)))
    epog.append(E_.mjd())
    dv_g.append(universe.evaluables.get('SC_dv').eval(E_))
    # Accumulated ΔV history over time.


# =============================================================================
# SECTION 17 — PLOT THE RESULTS
# =============================================================================

ax[0].plot(
    [x_[1] for x_ in XROT],   # X component (index 1, after MJD)
    [x_[2] for x_ in XROT],   # Y component (index 2)
    color='blue',
    linewidth=0.3
)
# Panel 0: XY projection of the Halo orbit in the SEROT frame.
# Indices: 0=MJD, 1=X, 2=Y, 3=Z, 4=VX, 5=VY, 6=VZ.

ax[1].plot(
    [x_[1] for x_ in XROT],   # X component
    [x_[3] for x_ in XROT],   # Z component (out-of-plane)
    color='blue',
    linewidth=0.3
)
# Panel 1: XZ projection (shows the Halo's vertical amplitude).

# ax[2].plot([x_[2] for x_ in XROT],[x_[3] for x_ in XROT] , color='blue', linewidth=0.3)
# (disabled) would be the YZ projection.

ax[2].plot(epog, dv_g)
# Panel 2: accumulated ΔV vs time (MJD). Expect a staircase: flat between
# maneuvers, with a jump at each station-keeping burn.


# =============================================================================
# SECTION 18 — FINAL ΔV CHECK (double check)
# =============================================================================

f = totalDV_halo(n_pt * n_orb, universe)
# Fresh totalDV_halo instance, to compare with the value registered in the universe.

print(f.eval(traj.point('end_point')))
print(universe.evaluables.get('TotalDV_halo').eval(traj.point('end_point')))
# If the two numbers match, the objective is consistent.


# =============================================================================
# SECTION 19 — SAVE AND SHOW THE PLOTS
# =============================================================================

plt.show()
# Open the interactive figure window (blocks until closed).

plt.savefig('halo.png', dpi=150, bbox_inches='tight')
# Save the panels as a high-resolution PNG in the current folder.

print("Plot saved as halo.png")


# =============================================================================
# LOGICAL-FLOW SUMMARY
# =============================================================================
#
#   [universe.yml] ──► cosmos.Universe
#                            │
#                   [x0_cr3bp, T_cr3bp]
#                            │
#                      config_halo()  ◄── [transfer_ephe.py]
#                            │
#                    ┌───────┴────────┐
#                 config_traj    config_prob
#                    │               │
#              cosmos.Trajectory  cosmos.Problem
#                    │               │
#                    └───────┬───────┘
#                       pg.problem()
#                            │
#                       pg.ipopt()
#                            │
#                      algo.evolve()   ← OPTIMIZATION
#                            │
#                    traj.compute()
#                            │
#               ┌────────────┴────────────┐
#          print ΔV               plots XY, XZ, ΔV(t)
#          print duration              │
#                                  plt.savefig()
# =============================================================================
