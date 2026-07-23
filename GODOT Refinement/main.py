# =============================================================================
# MAIN.PY — Full-ephemeris refinement of the CR3BP transfer (GODOT driver)
# =============================================================================
#
# End-to-end driver for the two-stage GODOT refinement:
#   Stage A (Sections 3-8)  : optimize the L2 Halo parking orbit.
#   Stage B (Sections 9-12) : optimize the Halo -> Moon flyby -> comet transfer,
#                             switching from Sun+Earth (ES_gravity) to the full
#                             Sun-Earth-Moon model (EMS_gravity).
#
# The CR3BP guess (states + arc durations) is read from python_inputs.txt, the
# .txt bridge file produced by the MATLAB pipeline. The refined trajectory is
# written to traj_R2_preopt.yml (guess) and traj_R2.yml (optimized).
#
# Environment notes:
#   - Run inside the conda env:  conda activate godot
#   - Cell separators for interactive execution:  # %%
#   - Clear console: Ctrl+L   |   reset workspace: reset   |   fullscreen: F11
#
# Possible improvements (open items):
#   - verify the gravitational constants match those used in MATLAB
#   - move the integration centre to the Moon when propagating near the Moon
#   - likewise use the Earth as centre when propagating up to the flyby
#   - number of nodes scaled by each arc length (partially done, see Section 5)
# =============================================================================

# =============================================================================
# SECTION 1 — IMPORTS
# =============================================================================

import numpy as np

from godot.core import num, tempo, astro, events, ipfwrap
from godot.core import autodif as ad
from godot.core.autodif import bridge as br
import godot.core.util as c_util
from godot.model import interface, frames, common, prop
from godot import cosmos
from godot.cosmos import util

import matplotlib.pyplot as plt
from ruamel import yaml
import time, os, copy
import pygmo as pg               # optimizer (IPOPT interface)
import midas                     # CR3BP library

from transfer_ephe import *
# Provides config_halo(), config_trajectory(), and physical constants
# (D_SUN_EARTH, GM_SUN, GM_EARTH, mu_se, T_SUN_EARTH).

from totalDV import *
# Provides totalDV_halo and totalDV_traj (objective functions).


# =============================================================================
# SECTION 1B — INPUT FILE PARSER
# =============================================================================

def _parse_inputs(filepath='python_inputs.txt'):
    """Read python_inputs.txt and return a dict of named parameters.

    Supports floats, lists of floats (bracket notation), and bare strings
    (used for the ISO date). Comment lines (#) and blank lines are skipped.
    """
    import ast
    params = {}
    with open(filepath, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' not in line:
                continue
            key, _, val = line.partition('=')
            key = key.strip()
            val = val.strip()
            try:
                params[key] = ast.literal_eval(val)
            except (ValueError, SyntaxError):
                params[key] = val   # date string falls through here
    return params

inp = _parse_inputs('python_inputs.txt')
# inp = _parse_inputs('python_inputs_serot.txt')

print("Inputs read from python_inputs.txt:", list(inp.keys()))


# =============================================================================
# SECTION 2 — UNIVERSE INITIALIZATION
# =============================================================================

c_util.suppressLogger()

config_uni = util.load_yaml('Universe/universe.yml')
universe = cosmos.Universe(config_uni)
# Loads Sun/Earth/Moon physical model with DE432 ephemeris and SC reference frames.
# Universe is re-initialized in Section 10 before the second optimization.




# =============================================================================
# SECTION 3 — HALO ORBIT PARAMETERS (from python_inputs.txt)
# =============================================================================

x0_halo_cr3bp = np.array(inp['state_pre_injection'])
# Halo state at m*T_halo before manifold injection (synodic CR3BP non-dim).

T_halo_cr3bp = float(inp['T_halo'])
# Halo orbit period in CR3BP non-dimensional time.

m = float(inp['m'])
# Fraction of Halo period the SC travels on-orbit before injecting onto the manifold.

# =============================================================================
# SECTION 4 — TRANSFER TRAJECTORY STATES & DURATIONS (from python_inputs.txt)
# =============================================================================
# Arc structure: Halo ─[inj]─► manifold ─[t_inj_dsm1]─► DSM1
#                ─[t_dsm1_flyby]─► Moon flyby ─[t_flyby_dsm2]─► DSM2
#                ─[t_dsm2_comet]─► comet

x_inj   = np.array(inp['state_after_injection'])
# CR3BP state just after manifold injection (start of arc injection→DSM1).

x_dsm1  = np.array(inp['state_after_dsm1'])
# CR3BP state just after DSM1 (start of arc DSM1→flyby).

x_flyby = np.array(inp['state_flyby'])
# CR3BP state at Moon flyby pericentre.

x_dsm2  = np.array(inp['state_after_dsm2'])
# CR3BP state just after DSM2 (start of arc DSM2→comet).

t_inj_dsm1   = float(inp['dur_inj_to_dsm1'])
# CR3BP non-dim duration: injection → DSM1.

t_dsm1_flyby = float(inp['dur_dsm1_to_flyby'])
# CR3BP non-dim duration: DSM1 → Moon flyby pericentre.

t_flyby_dsm2 = float(inp['dur_flyby_to_dsm2'])
# CR3BP non-dim duration: Moon flyby pericentre → DSM2.

t_dsm2_comet = float(inp['dur_dsm2_to_comet'])
# CR3BP non-dim duration: DSM2 → comet intercept.

dv = 0.2   # ΔV optimizer bound [km/s] (hardcoded)



# =============================================================================
# SECTION 5 — MISSION TIMELINE (derived from python_inputs.txt)
# =============================================================================
# All epochs are computed backwards from date_comet (end of trajectory).

date_comet = tempo.Epoch(inp['date_flyby'] + ' UTC')
# Target comet flyby epoch — end of the transfer trajectory.
print(f"Comet flyby date (end of traj): {date_comet}")

t_total_cr3bp        = t_inj_dsm1 + t_dsm1_flyby + t_flyby_dsm2 + t_dsm2_comet
tof_total_s          = t_total_cr3bp * T_SUN_EARTH / (2 * np.pi)
date_start_from_halo = date_comet - tof_total_s
# Injection epoch: SC leaves Halo onto the unstable manifold.
print(f"Injection epoch (start from Halo): {date_start_from_halo}")

period_s   = T_halo_cr3bp * T_SUN_EARTH / (2 * np.pi)
date_start = date_start_from_halo - m * period_s
# Start epoch: SC is on Halo orbit, m*T_halo before injection.

date_moon_flyby = date_start_from_halo + (t_inj_dsm1 + t_dsm1_flyby) * T_SUN_EARTH / (2 * np.pi)
# Moon flyby pericentre epoch (computed from timing chain).
print(f"Moon flyby epoch: {date_moon_flyby}")

# Node spacing [days between control points] per arc, in order:
#   [inj→DSM1, DSM1→flyby, flyby→DSM2, DSM2→comet].
# The number of control points on each arc is derived from its duration, so you
# refine finely in the delicate flyby phase and coarsely on the chill
# interplanetary cruise. Smaller value => more nodes => finer integration grid.
nodes_days_vec = [
    10.0,    # inj  → DSM1
    1.5,    # DSM1 → flyby   (fine: delicate lunar approach)
    0.5,    # flyby → DSM2   (fine: delicate lunar departure)
    30.0,   # DSM2 → comet   (coarse: interplanetary cruise)
]

r_moon  = universe.frames.vector6('Sun', 'Moon',  'ICRF', date_moon_flyby)
r_earth = universe.frames.vector6('Sun', 'Earth', 'ICRF', date_moon_flyby)
print(f"Moon  (ICRF, Sun-centered) at flyby: {r_moon[:3]}  km")
print(f"Earth (ICRF, Sun-centered) at flyby: {r_earth[:3]}  km")


# =============================================================================
# SECTION 6 — HALO OPTIMIZATION PARAMETERS
# =============================================================================

n_pt = 4   # control points per orbit
n_orb = 1  # 1 Halo orbit before manifold injection

dv_halo = totalDV_halo(n_pt * n_orb, universe)
universe.evaluables.add('TotalDV_halo', dv_halo)
# Registers the Halo station-keeping ΔV as objective function.

dx = 100000   # position perturbation bound: 100 000 km
dv = 0.2      # velocity perturbation bound: 0.2 km/s

scales_vect = [2000000, 870000, 600000, dv, dv, dv]
# Normalization scales for the optimizer variables.
# X/Y/Z differ because the Halo orbit is asymmetric in those directions.

delta_vect = np.array([dx, dx, dx, dv, dv, dv])
# Half-width of bounds around CR3BP guess: [guess - delta, guess + delta].

# =============================================================================
# SECTION 7 — FIRST OPTIMIZATION: HALO ORBIT
# =============================================================================

config_traj, config_prob = config_halo(
    x0_halo_cr3bp,  # CR3BP initial state
    T_halo_cr3bp,   # CR3BP period
    date_start,     # start epoch (computed backwards from flyby)
    n_pt,
    n_orb,
    scales_vect,
    delta_vect
)

traj = cosmos.Trajectory(universe, config_traj)

prob = cosmos.Problem(universe, [traj], config_prob, useGradient=True)

problem = pg.problem(prob)

tol_con = 1e-6
problem.c_tol = [tol_con] * problem.get_nc()

x0 = prob.get_x()
pop = pg.population(problem, 0)
pop.push_back(x0)

ip = pg.ipopt()
ip.set_numeric_option("tol", 1e-3)
# Loose tolerance: result is only used as initial guess for the transfer optimization.

algo = pg.algorithm(ip)
algo.set_verbosity(1)

pop = algo.evolve(pop)
# After this: traj holds the optimized Halo orbit.


# =============================================================================
# SECTION 8 — HANDOFF: EXTRACT HALO STATE AT INJECTION EPOCH
# =============================================================================
# Bridge between first and second optimization.
# The Halo orbit starting point is always x0_halo_cr3bp; m controls
# how far along the orbit (fraction of period) the SC travels before injection.

traj.compute(partials=False)

traj_up = traj.applyParameterChanges()
conf_update = util.deep_update(config_traj, traj_up)

dv_tot = universe.evaluables.get('TotalDV_halo').eval(traj.point('end_point'))
print(f'Total dv: {dv_tot * 1e3:.2e} m/s')

# Extract SC state at injection epoch in SEROT frame (Earth-centered, km, km/s).
x = universe.frames.vector6(
    'Earth',
    'SC_center',
    'SEROT',
    date_start_from_halo
)
# This optimized state feeds into config_trajectory() as the transfer arc start.
print(x)

# x = x0_halo_cr3bp

# =============================================================================
# SECTION 9 — BUILD TRANSFER TRAJECTORY
# =============================================================================

# New 4-arc structure: the arc-junction states and durations come straight from
# python_inputs.txt (parsed in Sections 3-4). Intermediate nodes are seeded by
# ES_gravity propagation inside config_trajectory().
conf_traj2, conf_prob, conf_prob2 = config_trajectory(
    x,                    # SC state at Halo injection (from Section 8)
    date_start_from_halo,
    x_inj,                # state after injection  (start of arc inj→DSM1)
    x_dsm1,               # state after DSM1        (start of arc DSM1→flyby)
    x_flyby,              # state at Moon flyby     (start of arc flyby→DSM2)
    x_dsm2,               # state after DSM2        (start of arc DSM2→comet)
    t_inj_dsm1,           # duration injection → DSM1   [adim]
    t_dsm1_flyby,         # duration DSM1 → flyby        [adim]
    t_flyby_dsm2,         # duration flyby → DSM2        [adim]
    t_dsm2_comet,         # duration DSM2 → comet        [adim]
    dv,                   # ΔV scale/bound
    nodes_days_vec,       # node spacing [days] per arc (finer near the flyby)
    False                 # with_man=False: no intermediate maneuvers on arcs
)


# =============================================================================
# SECTION 10 — UNIVERSE RE-INIT + SWITCH TO FULL EPHEMERIS (Sun-Earth-Moon)
# =============================================================================
# Fresh universe instance needed before the second optimization.
# The trajectory guess was built with the ES_gravity (Sun+Earth) dynamics to
# match how the input states were generated. For the refinement we switch every
# control point to the full Sun-Earth-Moon model (EMS_gravity), so the lunar
# flyby is now modelled by the Moon's gravity.

config_uni = util.load_yaml('Universe/universe.yml')
universe = cosmos.Universe(config_uni)

dv_traj = totalDV_traj(universe)
universe.evaluables.add('TotalDV_traj', dv_traj)

# Lunar-flyby altitude evaluable, used as an inequality constraint (>= 750 km)
# on the flyby node to keep the trajectory clear of the Moon's surface.
flyby_alt = flybyMoonAltitude(universe)
universe.evaluables.add('FlybyAlt', flyby_alt)

for _item in conf_traj2['timeline']:
    if _item.get('type') == 'control':
        _item['state'][0]['dynamics'] = 'EMS_gravity'


# =============================================================================
# SECTION 11 — SECOND OPTIMIZATION SETUP: TRANSFER TRAJECTORY
# =============================================================================

traj2 = cosmos.Trajectory(universe, conf_traj2)

prob2 = cosmos.Problem(universe, [traj2], conf_prob2, useGradient=True)

traj2.compute(partials=False)



# --- diagnostic: MINIMUM distance from the Moon along the propagated arc ---
# Not at the node, but by densely sampling the trajectory around the flyby.
_ep_fly = next(tempo.Epoch(str(_it['epoch']))
               for _it in conf_traj2['timeline']
               if _it.get('type') == 'control' and _it.get('name') == 'ctr_flybya')

# window of +/- 1 day around the flyby, step 2 minutes
_R_MOON = 1737.4
_dmin = 1e30
_t_dmin = None
for _dt_min in range(-1440, 1441, 2):          # from -1 day to +1 day, step 2 min
    _ep = _ep_fly + _dt_min * 60.0              # seconds
    try:
        _r = universe.frames.vector3('Moon', 'SC_center', 'ICRF', _ep)
        _d = float(np.linalg.norm(_r))
        if _d < _dmin:
            _dmin = _d
            _t_dmin = _dt_min
    except Exception:
        pass                                     # epoch outside the trajectory range

print('=== minimum distance from the Moon along the arc (guess) ===')
print('minimum distance [km] :', _dmin)
print('minimum altitude [km] :', _dmin - _R_MOON, '(<0 => INSIDE the Moon!)')
print('time (min from flyby node) :', _t_dmin)





# =============================================================================
# SECTION 11C — SAVE PRE-OPTIMIZATION TRAJECTORY (YAML)
# =============================================================================

traj_up_pre = traj2.applyParameterChanges()
conf_update_pre = util.deep_update(conf_traj2, traj_up_pre)
util.save_yaml(conf_update_pre, 'traj_R2_preopt.yml')
print("Salvato: traj_R2_preopt.yml")


# =============================================================================
# SECTION 11C-bis — PLOT PRE-OPTIMIZATION TRAJECTORY (Earth-centered SEROT)
# =============================================================================
# Samples SC_center and Moon in SEROT at N_plot epochs between first and last
# node, then overlays the multiple-shooting nodes. Frame: Earth-centered SEROT.
# Node positions are read directly from conf_update_pre (no vector6 for SC_center
# needed — avoids GODOT's path-finding limitation with non-SEROT frames).

def _val_to_float(v):
    return float(v.split()[0]) if isinstance(v, str) else float(v)

# --- collect node epochs and positions from config dict ---
_node_epochs = []
_node_pos    = []   # Earth-centered SEROT km
_node_names  = []
_man_nodes   = {}   # maneuver name -> ref node name

for _it in conf_update_pre['timeline']:
    if _it.get('type') == 'control':
        _ep_str = str(_it['epoch'])
        _vv     = _it['state'][0]['value']
        _node_epochs.append(tempo.Epoch(_ep_str))
        _node_pos.append([_val_to_float(_vv['pos_x']),
                          _val_to_float(_vv['pos_y']),
                          _val_to_float(_vv['pos_z'])])
        _node_names.append(_it['name'])
    elif _it.get('type') == 'manoeuvre':
        _man_nodes[_it['name']] = _it['config']['point']['reference']

_node_pos = np.array(_node_pos)

# --- sample trajectory and Moon densely between first and last node ---
# applyParameterChanges() (called above for the YAML save) invalidates the
# computed state of traj2. Recompute here so SC_center is connected again.
traj2.compute(partials=False)

_N      = 100000
_ep0    = _node_epochs[0] + 60.0
_ep1    = _node_epochs[-1] - 60.0
_tof_s  = float(_ep1 - _ep0)
_traj   = []
_moon_s = []

for _k in range(_N):
    _ep_k = _ep0 + _k * _tof_s / (_N - 1)
    _sv   = universe.frames.vector6('Earth', 'SC_center', 'SEROT', _ep_k)
    _sm   = universe.frames.vector6('Earth', 'Moon',      'SEROT', _ep_k)
    _traj.append(_sv[:3])
    _moon_s.append(_sm[:3])

_traj   = np.array(_traj)
_moon_s = np.array(_moon_s)

# Moon position at the flyby node (ctr_flybya)
_flyby_idx = _node_names.index('ctr_flybya') if 'ctr_flybya' in _node_names else None
if _flyby_idx is not None:
    _moon_flyby = universe.frames.vector6('Earth', 'Moon', 'SEROT', _node_epochs[_flyby_idx])[:3]

# --- identify special nodes for color coding ---
_junction_names = {'ctr_dsm1a', 'ctr_flybya', 'ctr_dsm2a'}
_man_ref_names  = set(_man_nodes.values())   # nodes that carry a maneuver

# --- plot (top view: X-Y plane in SEROT = plane of Sun-Earth orbit) ---
fig, ax = plt.subplots(figsize=(11, 9))

ax.plot(_traj[:,0], _traj[:,1], 'b-', lw=1.2, label='Trajectory (guess)', zorder=2)
ax.plot(_moon_s[:,0], _moon_s[:,1], color='gray', lw=0.7, ls='--',
        label='Moon track', zorder=1)

# Earth at origin
ax.scatter(0, 0, s=120, c='deepskyblue', zorder=6, label='Earth')

# Moon at flyby epoch
if _flyby_idx is not None:
    ax.scatter(_moon_flyby[0], _moon_flyby[1], s=80, c='silver',
               edgecolors='k', lw=0.8, zorder=6, label='Moon @ flyby')

# Multiple-shooting nodes
for _ni, (_nm, _np) in enumerate(zip(_node_names, _node_pos)):
    if _nm == 'ctr_flybya':
        _col, _sz, _mk = 'red', 90, 'D'
    elif _nm in _junction_names:
        _col, _sz, _mk = 'orange', 70, 's'
    elif _nm in _man_ref_names:
        _col, _sz, _mk = 'limegreen', 60, '^'
    else:
        _col, _sz, _mk = 'white', 40, 'o'
    ax.scatter(_np[0], _np[1], s=_sz, c=_col, marker=_mk,
               edgecolors='k', lw=0.6, zorder=5)

# Legend proxies for nodes
from matplotlib.lines import Line2D
_legend_extra = [
    Line2D([0],[0], marker='D', color='w', markerfacecolor='red',    markersize=9,  markeredgecolor='k', label='flyby node'),
    Line2D([0],[0], marker='s', color='w', markerfacecolor='orange', markersize=8,  markeredgecolor='k', label='DSM nodes'),
    Line2D([0],[0], marker='^', color='w', markerfacecolor='limegreen', markersize=8, markeredgecolor='k', label='maneuver nodes'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='white',  markersize=7,  markeredgecolor='k', label='internal nodes'),
]
_handles, _labels = ax.get_legend_handles_labels()
ax.legend(handles=_handles + _legend_extra, fontsize=8, loc='best')

ax.set_xlabel('X SEROT [km]')
ax.set_ylabel('Y SEROT [km]')
ax.set_title('Pre-opt guess — Earth-centered SEROT (top view)')
ax.set_aspect('equal')
ax.grid(True, lw=0.4, alpha=0.5)
plt.tight_layout()
plt.savefig('traj_preopt.png', dpi=150)
plt.show()
print('Plot saved:traj_preopt.png')

# --- ZOOM PLOT: Earth-Moon region ---
# Configurable center (km, SEROT) and half-widths (km).
# Default: Earth at origin, window large enough to contain the Moon track.
_zoom_cx  = 0.25*10**6   # centre X [km]
_zoom_cy  = 0.0   # centre Y [km]
_zoom_dx  = 1.5*10**6   # half-width  X [km]  (500 000 km ≈ 1.3 × Moon distance)
_zoom_dy  = 1.25*10**6   # half-width  Y [km]

fig2, ax2 = plt.subplots(figsize=(9, 9))

ax2.plot(_traj[:,0], _traj[:,1], 'b-', lw=1.2, label='Trajectory (guess)', zorder=2)
ax2.plot(_moon_s[:,0], _moon_s[:,1], color='gray', lw=0.7, ls='--',
         label='Moon track', zorder=1)

ax2.scatter(0, 0, s=180, c='deepskyblue', zorder=6, label='Earth')

if _flyby_idx is not None:
    ax2.scatter(_moon_flyby[0], _moon_flyby[1], s=120, c='silver',
                edgecolors='k', lw=0.8, zorder=6, label='Moon @ flyby')

for _ni, (_nm, _np) in enumerate(zip(_node_names, _node_pos)):
    if _nm == 'ctr_flybya':
        _col2, _sz2, _mk2 = 'red', 130, 'D'
    elif _nm in _junction_names:
        _col2, _sz2, _mk2 = 'orange', 90, 's'
    elif _nm in _man_ref_names:
        _col2, _sz2, _mk2 = 'limegreen', 80, '^'
    else:
        _col2, _sz2, _mk2 = 'white', 50, 'o'
    ax2.scatter(_np[0], _np[1], s=_sz2, c=_col2, marker=_mk2,
                edgecolors='k', lw=0.7, zorder=5)

_handles2, _labels2 = ax2.get_legend_handles_labels()
ax2.legend(handles=_handles2 + _legend_extra, fontsize=8, loc='best')

ax2.set_xlim(_zoom_cx - _zoom_dx, _zoom_cx + _zoom_dx)
ax2.set_ylim(_zoom_cy - _zoom_dy, _zoom_cy + _zoom_dy)
ax2.set_xlabel('X SEROT [km]')
ax2.set_ylabel('Y SEROT [km]')
ax2.set_title('Pre-opt guess — ZOOM Earth-Moon region (SEROT)')
ax2.set_aspect('equal')
ax2.grid(True, lw=0.4, alpha=0.5)
plt.tight_layout()
plt.savefig('traj_preopt_zoom.png', dpi=150)
plt.show()
print('Plot saved:traj_preopt_zoom.png')

# --- FLYBY ZOOM: Moon-centred, tight view of flyby geometry ---
# Half-widths in km: 20 000 km default gives ~10× Moon radius on each side.
_fly_dx = 20000.0   # half-width X [km]
_fly_dy = 20000.0   # half-width Y [km]
_R_moon = 1737.4    # Moon radius [km] — drawn as circle

if _flyby_idx is not None:
    _cx3 = _moon_flyby[0]
    _cy3 = _moon_flyby[1]

    fig3, ax3 = plt.subplots(figsize=(9, 9))

    ax3.plot(_traj[:,0], _traj[:,1], 'b-', lw=1.4, label='Trajectory (guess)', zorder=2)

    # Moon body (filled circle) and position marker
    _moon_circle = plt.Circle((_cx3, _cy3), _R_moon,
                               color='silver', ec='k', lw=0.8, zorder=4,
                               label='Moon (to scale)')
    ax3.add_patch(_moon_circle)
    ax3.scatter(_cx3, _cy3, s=60, c='silver', edgecolors='k', lw=0.8, zorder=5)

    # All multiple-shooting nodes (only those inside the window are visible)
    for _nm3, _np3 in zip(_node_names, _node_pos):
        if _nm3 == 'ctr_flybya':
            _c3, _s3, _m3 = 'red', 150, 'D'
        elif _nm3 in _junction_names:
            _c3, _s3, _m3 = 'orange', 100, 's'
        elif _nm3 in _man_ref_names:
            _c3, _s3, _m3 = 'limegreen', 80, '^'
        else:
            _c3, _s3, _m3 = 'white', 50, 'o'
        ax3.scatter(_np3[0], _np3[1], s=_s3, c=_c3, marker=_m3,
                    edgecolors='k', lw=0.7, zorder=6)

    # Minimum altitude circle (500 km above surface) as reference
    _h_min = 500.0
    _alt_circle = plt.Circle((_cx3, _cy3), _R_moon + _h_min,
                              color='none', ec='red', lw=0.8, ls='--', zorder=3,
                              label='h = %.0f km' % _h_min)
    ax3.add_patch(_alt_circle)

    _handles3, _ = ax3.get_legend_handles_labels()
    ax3.legend(handles=_handles3 + _legend_extra, fontsize=8, loc='best')

    ax3.set_xlim(_cx3 - _fly_dx, _cx3 + _fly_dx)
    ax3.set_ylim(_cy3 - _fly_dy, _cy3 + _fly_dy)
    ax3.set_xlabel('X SEROT [km]')
    ax3.set_ylabel('Y SEROT [km]')
    ax3.set_title('Pre-opt guess — ZOOM flyby geometry (Moon-centred SEROT)')
    ax3.set_aspect('equal')
    ax3.grid(True, lw=0.4, alpha=0.5)
    plt.tight_layout()
    plt.savefig('traj_preopt_flyby.png', dpi=150)
    plt.show()
    print('Plot saved:traj_preopt_flyby.png')
else:
    print('Node ctr_flybya not found — flyby plot skipped.')


# =============================================================================
# SECTION 11D — LAUNCH SECOND OPTIMIZATION (full-ephemeris refinement)
# =============================================================================

problem2 = pg.problem(prob2)

tol_con = 1e-3
problem2.c_tol = [tol_con] * problem2.get_nc()

x0 = prob2.get_x()
pop = pg.population(problem2, 0)
pop.push_back(x0)

# --- diagnostic: check the flyby-altitude constraint at the initial guess ---
# Evaluate fitness once (forces traj2 computation so SC_center is connected),
# then read the raw lunar altitude at the flyby node directly. The guess is
# feasible w.r.t. this constraint if the altitude is >= 750 km.
_fit0 = problem2.fitness(x0)
print('--- flyby constraint diagnostic ---')
print('total number of constraints :', problem2.get_nc())
_ep_flyby = next(tempo.Epoch(str(_it['epoch']))
                 for _it in conf_traj2['timeline']
                 if _it.get('type') == 'control' and _it.get('name') == 'ctr_flybya')
_r_fly = universe.frames.vector3('Moon', 'SC_center', 'ICRF', _ep_flyby)
_alt_fly = float(np.linalg.norm(_r_fly)) - 1737.4
print('flyby altitude at the guess [km] :', _alt_fly, '(>= 750 => feasible)')

ip = pg.ipopt()
ip.set_numeric_option("tol", 1e-3)
algo = pg.algorithm(ip)
algo.set_verbosity(1)
pop = algo.evolve(pop)


# algo = pg.scipy_optimize(method='SLSQP')
# algo.set_verbosity(1)
# pop = algo.evolve(pop)

# wh = ppnf.worhp(screen_output=True)
# algo = pg.algorithm(wh)
# pop = algo.evolve(pop)




# =============================================================================
# SECTION 12 — SAVE OPTIMIZED TRAJECTORY (YAML)
# =============================================================================

traj2.compute(partials=False)

traj_up = traj2.applyParameterChanges()
conf_update = util.deep_update(conf_traj2, traj_up)
util.save_yaml(conf_update, 'traj_R2.yml')
print("Salvato: traj_R2.yml")
