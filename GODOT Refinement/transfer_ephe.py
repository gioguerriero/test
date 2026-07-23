# =============================================================================
# TRANSFER_EPHE.PY — Trajectory configuration for optimization
# =============================================================================
#
# Contains two functions:
#   config_halo()        → builds the Halo orbit problem (ephemeris correction of CR3BP solution)
#   config_trajectory()  → builds the transfer trajectory problem (Halo → Moon flyby → comet)
#
# MULTIPLE SHOOTING:
#   Each trajectory arc is split into N short segments. Each segment has:
#     - a control point (ctr): state [pos, vel] at the node
#     - an optional maneuver (man): impulsive ΔV at that node
#     - a match constraint: requires that forward propagation from node i and
#       backward propagation from node i+1 meet at the midpoint
#
#   ctr0 ──man0──►──────match0──────◄── ctr1 ──man1──►──────match1──────◄── ctr2
#              forward ──►                    ◄── backward
#                             ^ meeting point (midpoint)
#
#   Match constraints = equations the optimizer must satisfy.
#   Free variables (positions, velocities, ΔVs) = parameters it can change.
#
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

import numpy as np               # (duplicate, harmless)
import matplotlib.pyplot as plt
from ruamel import yaml
import time, os, copy
import pygmo as pg
import pygmo_plugins_nonfree as ppnf
import midas                     # CR3BP library

from aux_fun import ctr, man, match
# Building blocks from aux_fun.py:
#   ctr()   → creates a control point (state node)
#   man()   → creates an impulsive maneuver (ΔV)
#   match() → creates a continuity constraint between two arcs


# =============================================================================
# SECTION 2 — SUN-EARTH PHYSICAL CONSTANTS
# =============================================================================

# D_SUN_EARTH = 149597870  (old value)
D_SUN_EARTH = 149600000   # Mean Sun-Earth distance in km (CR3BP length unit)

GM_SUN = 132758540018     # Sun gravitational parameter [km³/s²]

# GM_EARTH = 398600.4418  (old value used by Erwan)
GM_EARTH = 398702.4418    # Earth gravitational parameter [km³/s²]

mu_se = GM_EARTH / (GM_EARTH + GM_SUN)
# CR3BP mass parameter: ratio of Earth mass to total system mass (~3e-6).
# Fully defines the non-dimensional CR3BP dynamics.

T_SUN_EARTH = 2 * np.pi * np.sqrt(D_SUN_EARTH**3 / (GM_EARTH + GM_SUN))
# Earth orbital period [s] (Kepler's third law).
# Converts CR3BP dimensionless time to seconds: t[s] = t_cr3bp * T_SUN_EARTH / (2π)


# =============================================================================
# FUNCTION 1: config_halo()
# =============================================================================
# Takes a CR3BP Halo orbit solution and builds the GODOT optimization problem
# to correct it into a true ephemeris orbit (minimizing station-keeping ΔV).
# The CR3BP solution drifts when propagated with real physics — this fixes it.
#
# Args:
#   x0_cr3bp   : CR3BP initial state [x,y,z,vx,vy,vz] (non-dimensional)
#   T_CR3BP    : CR3BP Halo orbit period (non-dimensional)
#   date_start : start epoch (GODOT tempo.Epoch)
#   n_pt       : control points per orbit
#   n_orb      : number of complete Halo orbits
#   scales_vect: normalization scales for optimizer variables
#   delta_vect : bound half-widths around CR3BP guess
#
# Returns: (config_traj, config_prob) — GODOT-ready dicts for Trajectory and Problem
# =============================================================================

def config_halo(x0_cr3bp, T_CR3BP, date_start, n_pt, n_orb, scales_vect, delta_vect):

    # -------------------------------------------------------------------------
    # STEP 1: Initialization
    # -------------------------------------------------------------------------

    CRTBP = midas.astro.crtbp.AdimensionalCRTBP(mu_se)
    IntegratorCRTBP = midas.astro.crtbp.AdimensionalIntegratorCRTBP(mu_se)
    # CR3BP system and integrator for propagating the reference Halo orbit.

    config_uni = util.load_yaml('Universe/universe.yml')
    universe = cosmos.Universe(config_uni)
    # Real ephemeris model: needed to compute the actual Sun-Earth distance (lu)
    # at each epoch for the CR3BP → physical coordinates conversion.

    # -------------------------------------------------------------------------
    # STEP 2: Build time grid and propagate reference states
    # -------------------------------------------------------------------------

    period = T_CR3BP * T_SUN_EARTH / (2 * np.pi * 86400)  # Halo period in days
    dt = period / n_pt                                      # time step between control points [days]

    mjd_vect_ref = date_start.mjd() + np.arange(0, n_orb * period + dt, dt)
    # Epochs of control points in MJD. n_pt*n_orb+1 nodes total
    # (+1 so the closing node is included).

    x0_vect_ref = []

    for m_ in mjd_vect_ref:
        # Propagate CR3BP from x0_cr3bp to the elapsed time (mod period, since Halo is periodic).
        dt_cr3bp = (m_ - date_start.mjd()) % period
        dt_cr3bp = dt_cr3bp * 86400 * 2 * np.pi / T_SUN_EARTH  # days → CR3BP non-dim

        [t_cr3bp, x_cr3bp] = IntegratorCRTBP.integrate(x0_cr3bp, dt_cr3bp)

        x0_rot = np.copy(x_cr3bp[-1])  # final state of integration

        # Convert CR3BP non-dimensional → physical coordinates (km, km/s):
        lu = universe.frames.distance('Sun', 'Earth', tempo.Epoch(str(m_) + ' TDB'))
        # lu = actual Sun-Earth distance at this epoch (varies ~±2.5% over the year).
        tu = np.sqrt(lu**3 / (GM_SUN + GM_EARTH))

        x0_rot[0] -= 1 - mu_se   # shift barycenter-centered → Earth-centered (x=0 at Earth)
        x0_rot[:3] *= lu          # position: non-dim → km
        x0_rot[3:] *= lu / tu     # velocity: non-dim → km/s

        x0_vect_ref.append(x0_rot)

    # -------------------------------------------------------------------------
    # STEP 3: Build GODOT timeline (config_traj)
    # -------------------------------------------------------------------------

    config_traj = {
        'settings': {'relTol': 1e-12, 'steps': 10000000},
        'setup': [
            {'name': 'SC',
             'type': 'group',
             'spacecraft': 'SC',
             'input': [
                 {'name': 'center', 'type': 'point'},
                 {'name': 'dv', 'type': 'scalar', 'unit': 'm/s'}
             ]}]
    }

    TT = []
    i_ = 0
    x_point_ref = []  # physical states of control points, used later for bounds

    for m_, x_ in zip(mjd_vect_ref, x0_vect_ref):
        ct = ctr('ctr' + str(i_), tempo.Epoch(str(m_) + ' TDB'), x_)
        TT.append(ct)

        if i_ < n_orb * n_pt:
            # Add maneuver + match for all nodes except the closing one.
            mn = man('man' + str(i_), 'ctr' + str(i_), 0,
                     np.array([0.0001, 0.0001, 0.0001]))
            # Initial ΔV guess 0.1 m/s — optimizer will change this.
            TT.append(mn)

            mtch = match('match' + str(i_), 'ctr' + str(i_), dt * 86400 / 2,
                         'ctr' + str(i_ + 1), dt * 86400 / 2)
            # Meeting point at midspan (dt/2 from each side).
            TT.append(mtch)

        i_ += 1
        x_point_ref.append(x_)

    end_point = {
        'type': 'point',
        'input': 'SC',
        'name': 'end_point',
        'point': {'reference': 'ctr' + str(i_ - 1), 'dt': '1.0 s'}
    }
    # 1 s after the last node — anchor for evaluating TotalDV_halo.
    TT.append(end_point)

    config_traj['timeline'] = TT

    # -------------------------------------------------------------------------
    # STEP 4: Build optimization problem (config_prob)
    # -------------------------------------------------------------------------
    # Free variables are identified by GODOT name strings, e.g.:
    #   'ctr2_SC_center_pos_x'  → position X of control point 2
    #   'man3_dv_y'             → Y component of maneuver 3 ΔV
    #   'match1_right_dt'       → right-side timing of match constraint 1
    #   'ctr1_dt'               → timing offset of control point 1

    config_prob = {
        'objective': {
            'type': 'minimise',
            'point': 'end_point',
            'value': 'TotalDV_halo',
            'scale': 1e6    # TotalDV_halo ~ 1e-6 km/s → ×1e6 → O(1) for optimizer
        },
        'parameters': {'free': []}
    }

    FF = []   # free variable names (GODOT strings)
    SS = {}   # scales: name → value with units
    BB = {}   # bounds: name → [min, max]

    cart_vect  = ['pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z']
    units_vect = [' km',   ' km',   ' km',   ' km/s', ' km/s', ' km/s']

    for i_ in range(n_pt * n_orb):
        if not (i_ == 0):
            # ctr0 has fixed epoch; all others get a free timing offset.
            FF.append('ctr' + str(i_) + '_dt')
            FF.append('ctr' + str(i_) + '_SC_dv')
            SS['ctr' + str(i_) + '_dt'] = '1 day'
            SS['ctr' + str(i_) + '_SC_dv'] = '1 mm/s'
            BB['ctr' + str(i_) + '_dt'] = ['-1 day', '1 day']
            BB['ctr' + str(i_) + '_SC_dv'] = ['0 m/s', '1000 m/s']

        for c_, s_, u_, b_, d_ in zip(cart_vect, scales_vect, units_vect,
                                       x_point_ref[i_], delta_vect):
            FF.append('ctr' + str(i_) + '_SC_center_' + c_)
            SS['ctr' + str(i_) + '_SC_center_' + c_] = str(s_) + u_
            BB['ctr' + str(i_) + '_SC_center_' + c_] = [
                str(b_ - d_) + u_,   # lower: CR3BP reference − delta
                str(b_ + d_) + u_    # upper: CR3BP reference + delta
            ]

        # Station-keeping maneuver ΔV components.
        FF.append('man' + str(i_) + '_dv_x')
        FF.append('man' + str(i_) + '_dv_y')
        FF.append('man' + str(i_) + '_dv_z')
        SS['man' + str(i_) + '_dv_x'] = '1 mm/s'
        SS['man' + str(i_) + '_dv_y'] = '1 mm/s'
        SS['man' + str(i_) + '_dv_z'] = '1 mm/s'
        BB['man' + str(i_) + '_dv_x'] = ['-1 m/s', '1 m/s']
        BB['man' + str(i_) + '_dv_y'] = ['-1 m/s', '1 m/s']
        BB['man' + str(i_) + '_dv_z'] = ['-1 m/s', '1 m/s']

        FF.append('match' + str(i_) + '_right_dt')
        SS['match' + str(i_) + '_right_dt'] = '1 day'
        BB['match' + str(i_) + '_right_dt'] = ['1 s', str(period / n_pt + 1) + ' day']
        # Lower bound 1 s avoids dt=0 (numerical degeneracy).

    # Closing node: bounds centered on x_point_ref[0] (first node) to enforce periodicity.
    i_ = n_pt * n_orb
    FF.append('ctr' + str(i_) + '_dt')
    FF.append('ctr' + str(i_) + '_SC_dv')
    SS['ctr' + str(i_) + '_dt'] = '1 day'
    SS['ctr' + str(i_) + '_SC_dv'] = '1 mm/s'
    BB['ctr' + str(i_) + '_dt'] = ['-1 day', '1 day']
    BB['ctr' + str(i_) + '_SC_dv'] = ['0 m/s', '100 m/s']

    for c_, s_, u_, b_, d_ in zip(cart_vect, scales_vect, units_vect,
                                   x_point_ref[0], delta_vect):
        FF.append('ctr' + str(i_) + '_SC_center_' + c_)
        SS['ctr' + str(i_) + '_SC_center_' + c_] = str(s_) + u_
        BB['ctr' + str(i_) + '_SC_center_' + c_] = [str(b_ - d_) + u_,
                                                      str(b_ + d_) + u_]

    config_prob['parameters']['free'] = FF
    config_prob['parameters']['scales'] = SS
    config_prob['parameters']['bounds'] = BB

    return config_traj, config_prob


# =============================================================================
# FUNCTION 2: config_trajectory()
# =============================================================================
# Builds the transfer trajectory from the Halo orbit to the comet.
#
# Arc structure (4 arcs; per-arc node count derived from nodes_vec):
#   ctr0a ─injection─► [arc A] ─DSM1─► [arc B] ─(flyby)─► [arc C] ─DSM2─► [arc D] ─► end_pointa
#
#   arc A : inj → DSM1
#   arc B : DSM1 → Moon flyby (explicit flyby node ctr_flybya, no ΔV)
#   arc C : flyby → DSM2
#   arc D : DSM2 → comet
#
# Arcs A/B/C carry n internal propagated nodes + an explicit junction node
# (ctr_dsm1a / ctr_flybya / ctr_dsm2a); arc D ends at the comet node. Node
# spacing per arc is set by nodes_vec [days/node]. Integration center (coi)
# switches Earth → Moon (flyby) → Sun (from the arc-C midpoint onward).
#
# Args:
#   x0_halo_ephe : SC state at Halo injection (physical: km, km/s — from Section 8)
#   date_start   : injection epoch
#   x1..x3_cr3bp : CR3BP reference states at arc junctions (non-dimensional)
#   t1..t3_cr3bp : arc durations (CR3BP non-dimensional)
#   dv_cr3bp     : ΔV scale (not used directly here)
#   nodes_vec    : [days_between_nodes] per arc, length 4:
#                  [inj→DSM1, DSM1→flyby, flyby→DSM2, DSM2→comet]. The number of
#                  control points on each arc is derived from its duration so the
#                  node spacing matches the requested value — finer near the
#                  flyby, coarser on the interplanetary cruise.
#   with_man     : if True, adds intermediate maneuvers inside arcs
#
# Returns: (config_traj, config_prob, config_prob2)
#   config_prob  → feasibility problem (match constraints only, no objective)
#   config_prob2 → optimization problem (minimizes TotalDV_traj)
# =============================================================================

def config_trajectory(x0_halo_ephe, date_start,
                      x_inj, x_dsm1, x_flyby, x_dsm2,
                      t_inj_dsm1, t_dsm1_flyby, t_flyby_dsm2, t_dsm2_comet,
                      dv_cr3bp, nodes_vec, with_man=False):
    # NEW STRUCTURE (4 arcs, fewer nodes than the previous 6-arc CR3BP version):
    #   ctr0a ─injection─► [arc A: inj→DSM1] ─DSM1─► [arc B: DSM1→flyby]
    #         ─(flyby, no man)─► [arc C: flyby→DSM2] ─DSM2─► [arc D: DSM2→comet] ─► end_pointa
    #
    # The arc-junction states (x_inj, x_dsm1, x_flyby, x_dsm2) come from the input
    # .txt and are expressed in synodic adimensional coordinates, BUT they were
    # generated with a Sun-Earth full-ephemeris dynamics (NOT CR3BP). Therefore the
    # intermediate nodes inside each arc are seeded by propagating the arc-start
    # state with the ES_gravity (Sun + Earth) dynamics, not with the CR3BP
    # integrator: CR3BP would land the comet at a different point.

    # -------------------------------------------------------------------------
    # STEP 1: Base trajectory config
    # -------------------------------------------------------------------------

    config_traj = {
        'settings': {'relTol': 1e-12, 'steps': 10000000},
        'setup': [{'name': 'SC', 'type': 'group', 'spacecraft': 'SC',
                   'input': [{'name': 'center', 'type': 'point'},
                              {'name': 'dv', 'type': 'scalar', 'unit': 'm/s'}]}]
    }
    # Same structure as config_halo(). Control points use suffix 'a' (ctr0a, ctr1a, ...)
    # to distinguish them from the Halo orbit nodes (ctr0, ctr1, ...).

    # -------------------------------------------------------------------------
    # STEP 2: Define two optimization problems
    # -------------------------------------------------------------------------
    # Two-phase strategy:
    #   Phase A (config_prob):  feasibility only — assemble arcs into a continuous trajectory
    #   Phase B (config_prob2): minimize TotalDV_traj starting from Phase A solution
    # With epochs fixed, FF (Phase A) and FF2 (Phase B) carry the same free variables.

    config_prob = {'parameters': {'free': []}}

    config_prob2 = {
        'objective': {
            'type': 'minimise',
            'point': 'end_pointa',
            'value': 'TotalDV_traj',   # injection ΔV + DSM1 + DSM2
            'scale': 1e3               # TotalDV_traj ~ 0.6 km/s → ×1e3 → ~600
        },
        # Inequality constraint: SC altitude above the lunar surface at the
        # flyby node (lunar pericentre) must stay >= 750 km. Prevents the
        # optimiser from sending the propagated arc through the Moon
        # ("BodyCache evalAcc, distance ... is zero").
        # Schema (opt/problem.json, detailed form): expression/type/value/point/scale.
        
        # 'constraints': ['flyby_alt; FlybyAlt > 750 @ ctr_flybya | 0.001'],
        
        #'constraints': [{
        #     'expression': 'FlybyAlt',     # evaluable: SC altitude above Moon [km]
        #    'type':       'greater',      # FlybyAlt >= value
        #    'value':      '750 km',
        #    'point':      'ctr_flybya',   # evaluated at the flyby (pericentre) node
        #    'scale':      1e-3,           # km-scale residual -> O(1) for the optimiser
        #}],
        'parameters': {'free': []}
    }
    


    ref_state       = {}   # node name -> physical SEROT state (used for bounds)
    full_free_nodes = []   # nodes with pos+vel+dv free (internal + DSM/flyby boundaries)
    man_names       = []   # names of optional intermediate maneuvers (with_man)
    FF = []         # free variables for Phase A
    FF2 = []        # free variables for Phase B

    # -------------------------------------------------------------------------
    # STEP 3: Initialization + helpers
    # -------------------------------------------------------------------------

    uni_config = cosmos.util.load_yaml('Universe/universe.yml')
    uni = cosmos.Universe(uni_config)
    # Physical universe used ONLY for the adim→physical conversion (distance
    # queries). The ES_gravity propagation below uses its own fresh universe per
    # arc: attaching a cosmos.Trajectory registers the 'SC_center' point on the
    # universe, so a single universe cannot host several trajectories at once
    # (this is why main.py re-initialises the universe between traj and traj2).

    # Conversion factor: adimensional CR3BP time → seconds.
    T_CONV = np.sqrt(D_SUN_EARTH**3 / (GM_SUN + GM_EARTH))   # = T_SUN_EARTH / (2π)

    def _adim_to_phys(x_adim, epoch):
        # Synodic adimensional (barycenter-centered) → physical km, km/s in the
        # Earth-centered SEROT frame, evaluated at the given epoch.
        lu = uni.frames.distance('Sun', 'Earth', epoch)
        tu = np.sqrt(lu**3 / (GM_SUN + GM_EARTH))
        x = np.array(x_adim, dtype=float)
        x[0] -= 1 - mu_se          # barycenter-centered → Earth-centered
        x[:3] *= lu                # position: non-dim → km
        x[3:] *= lu / tu           # velocity: non-dim → km/s
        return x
        # return x_adim

    def _propagate_es(x_phys, epoch0, sample_epochs):
        # Propagate an Earth-centered SEROT physical state with the ES_gravity
        # (Sun + Earth) dynamics and return the SEROT states at sample_epochs.
        # Uses a throw-away cosmos.Trajectory: one seed control point + an end
        # point that defines the integration span.
        # Integrate slightly past the last sample so every query epoch (incl. the
        # last one) sits STRICTLY inside the propagated span. Without the buffer
        # the last sample lands exactly on the integration end and the float
        # round-trip through str(span_s) can push the query a fraction of a second
        # past it -> "Could not connect Point Earth to Point SC_center".
        span_s = max(float(ep - epoch0) for ep in sample_epochs) + 60.0
        cfg = {
            'settings': {'relTol': 1e-12, 'steps': 10000000},
            'setup': [{'name': 'SC', 'type': 'group', 'spacecraft': 'SC',
                       'input': [{'name': 'center', 'type': 'point'},
                                 {'name': 'dv', 'type': 'scalar', 'unit': 'm/s'}]}],
            'timeline': [
                ctr('seed', epoch0, x_phys, dynamics='ES_gravity'), # ex ES
                {'type': 'point', 'input': 'SC', 'name': 'seed_end',
                 'point': {'reference': 'seed', 'dt': str(span_s) + ' s'}}
            ]
        }
        prop_uni = cosmos.Universe(cosmos.util.load_yaml('Universe/universe.yml'))
        tmp = cosmos.Trajectory(prop_uni, cfg)
        tmp.compute(partials=False)
        return [prop_uni.frames.vector6('Earth', 'SC_center', 'SEROT', ep)
                for ep in sample_epochs]

    TT = []

    # -------------------------------------------------------------------------
    # STEP 4: Departure point (on Halo orbit) + injection maneuver
    # -------------------------------------------------------------------------

    ctr_departure = ctr('ctr0a', date_start, x0_halo_ephe, dynamics='ES_gravity')
    TT.append(ctr_departure)

    x_inj_phys = _adim_to_phys(x_inj, date_start)
    dv1_vect_kms = x_inj_phys[3:] - x0_halo_ephe[3:]
    # Injection ΔV = post-injection velocity (x_inj, converted to km/s) minus
    # current Halo velocity. This kicks the SC onto the transfer arc.

    man_injection = man('man_injection', 'ctr0a', 0, dv1_vect_kms)
    TT.append(man_injection)

    # -------------------------------------------------------------------------
    # STEP 5-12: Build the 4 transfer arcs with ES_gravity propagation
    # -------------------------------------------------------------------------
    # Node layout (departure ctr0a + per-arc nodes derived from nodes_vec):
    #   ctr0a ─inj─► [arc A: nA internal] ─► ctr_dsm1a ─man_dsm1─►
    #         [arc B: nB internal] ─► ctr_flybya ─(no man)─►
    #         [arc C: nC internal] ─► ctr_dsm2a ─man_dsm2─►
    #         [arc D: nD nodes, last = comet]
    #
    # Per-arc node count: nodes_vec[k] is the requested spacing [days between
    #   nodes] for arc k. Arc duration / spacing gives the number of segments;
    #   nA/nB/nC are (segments-1) internal nodes (the boundary node completes the
    #   arc), nD is the number of arc-D nodes. This lets you refine finely near
    #   the flyby and coarsely on the interplanetary cruise.
    #
    # Arcs A, B, C: n INTERNAL propagated nodes + 1 EXPLICIT boundary node at the
    #   arc end (dt = t_arc / (n+1), so the boundary sits exactly at the junction
    #   epoch). The boundary node is the user-supplied input state (x_dsm1 /
    #   x_flyby / x_dsm2) — this is what lets us pin the lunar-flyby geometry
    #   instead of leaving it as a free propagation endpoint.
    # Arc D: n nodes, the last one being the comet (no extra boundary node; the
    #   comet target is not a user input).
    #
    # Internal + arc-D nodes keep sequential numeric names ctr1a … (comet = last);
    # the 3 junction nodes get descriptive names (ctr_dsm1a / ctr_flybya / ctr_dsm2a).
    #
    # TOF handling at a junction: the boundary node carries TWO matches with
    # DIFFERENT segment lengths — the closing match of the arc that ENDS there
    # (dt of that arc) and the leading match of the arc that STARTS there (dt of
    # the next arc). Each match below uses the dt of its own arc, so the two
    # sides of every junction are independently consistent.

    # Center-of-integration (coi) policy — improves integration conditioning by
    # using the locally dominant body as the propagation center:
    #   - everything before the lunar flyby (arcs A, B, internal nodes)  -> Earth
    #   - the flyby node ctr_flybya (lunar pericentre)                   -> Moon
    #   - arc C up to its MIDPOINT (still Earth-dominated)               -> Earth
    #   - arc C from its midpoint onward + DSM2 boundary + arc D         -> Sun
    # The Sun-centric leg starts at the midpoint between the flyby and DSM2
    # (waiting until DSM2 would switch too late, deep in the heliocentric regime).
    # ('point' stays Earth everywhere — see ctr() in aux_fun: only the
    # integration center changes, not the state representation.)

    # Helper: internal-node count for an A/B/C arc from its duration + spacing.
    def _n_internal(t_arc_adim, days_per_node):
        arc_days = t_arc_adim * T_CONV / 86400.0
        return max(1, int(round(arc_days / days_per_node)) - 1)

    # Boundary spec for arcs A, B, C:
    #   (arc-start adim state, arc duration [adim], boundary node name,
    #    boundary-input adim state, maneuver-after name or None, tight-pos flag,
    #    coi for internal nodes ('split' = Earth→Sun at the arc midpoint),
    #    coi for the boundary node, node spacing [days])
    abc_arcs = [
        # (x_inj,   t_inj_dsm1,   'ctr_dsm1a',  x_dsm1,  'man_dsm1', False, 'Earth', 'Earth', nodes_vec[0]),
        # (x_dsm1,  t_dsm1_flyby, 'ctr_flybya', x_flyby, None,       True,  'Earth', 'Moon',  nodes_vec[1]),   # flyby: no ΔV, tight, Moon-centric
        # (x_flyby, t_flyby_dsm2, 'ctr_dsm2a',  x_dsm2,  'man_dsm2', False, 'split', 'Sun',   nodes_vec[2]),   # Sun from arc-C midpoint
        
        (x_inj,   t_inj_dsm1,   'ctr_dsm1a',  x_dsm1,  'man_dsm1', False, 'Earth', 'Earth', nodes_vec[0]),
        (x_dsm1,  t_dsm1_flyby, 'ctr_flybya', x_flyby, None,       True,  'Earth', 'Earth',  nodes_vec[1]),   # flyby: no ΔV, tight, Moon-centric
        (x_flyby, t_flyby_dsm2, 'ctr_dsm2a',  x_dsm2,  'man_dsm2', False, 'Earth', 'Earth',   nodes_vec[2]),   # Sun from arc-C midpoint
        
    ]

    gidx            = 0          # running index of numeric-named (internal/arc-D) nodes
    prev_node       = 'ctr0a'    # node to chain the next leading match from
    epoch_arc_start = date_start

    for x_arc_adim, t_arc_adim, bnd_name, x_bnd_adim, man_after, tight, coi_int, coi_bnd, days_node in abc_arcs:

        n_pt_arc = _n_internal(t_arc_adim, days_node)   # internal nodes on this arc
        t_arc_s  = t_arc_adim * T_CONV
        dt_s     = t_arc_s / (n_pt_arc + 1)   # n internal nodes + boundary => n+1 segments

        # Seed: propagate the arc start at the n internal epochs PLUS the
        # boundary epoch (the last sample = propagated arrival at the junction).
        x_arc_phys    = _adim_to_phys(x_arc_adim, epoch_arc_start)
        sample_epochs = [epoch_arc_start + i * dt_s for i in range(1, n_pt_arc + 2)]
        arc_states    = _propagate_es(x_arc_phys, epoch_arc_start, sample_epochs)
        arrival_state = arc_states[n_pt_arc]      # propagated state at the boundary epoch
        epoch_bnd     = sample_epochs[n_pt_arc]

        # Leading match: previous node -> first internal node (THIS arc's dt).
        # If prev_node carries a maneuver (DSM), forward integration applies it.
        TT.append(match('match_' + prev_node,
                        prev_node, dt_s / 2,
                        'ctr' + str(gidx + 1) + 'a', dt_s / 2))

        # Internal nodes ctr_{gidx+1}a … ctr_{gidx+n_pt_arc}a.
        for i in range(1, n_pt_arc + 1):
            gidx += 1
            name  = 'ctr' + str(gidx) + 'a'
            x_loc = arc_states[i - 1]

            # 'split' arcs (arc C): Earth before the arc midpoint, Sun after.
            if coi_int == 'split':
                node_coi = 'Sun' if (i * dt_s >= t_arc_s / 2) else 'Earth'
            else:
                node_coi = coi_int

            ref_state[name] = x_loc
            full_free_nodes.append(name)
            TT.append(ctr(name, sample_epochs[i - 1], x_loc,
                          dynamics='ES_gravity', coi=node_coi))

            if with_man:
                man_names.append('man' + str(gidx) + 'a')
                TT.append(man('man' + str(gidx) + 'a', name, 0,
                              np.array([0.00001, 0.00001, 0.00001])))

            # Match to next node: next internal node, or this arc's boundary.
            nxt = ('ctr' + str(gidx + 1) + 'a') if i < n_pt_arc else bnd_name
            TT.append(match('match' + str(gidx) + 'a', name, dt_s / 2, nxt, dt_s / 2))

        # Explicit boundary node (the user-supplied junction state).
        x_bnd_phys = _adim_to_phys(x_bnd_adim, epoch_bnd)
        if tight:
            # Flyby: no maneuver, so arrival == departure. Store x_flyby fully —
            # the closing match of this arc absorbs the propagation residual.
            x_bnd_node = x_bnd_phys
        else:
            # DSM: the stored value is the PRE-maneuver (arrival) state, required
            # for the match mechanism (backward integration into this arc uses
            # the arrival velocity; the maneuver below is applied forward into the
            # next arc). Position is taken from the input (continuous across the
            # DSM); velocity is the propagated arrival.
            x_bnd_node = np.concatenate([x_bnd_phys[:3], arrival_state[3:]])

        ref_state[bnd_name] = x_bnd_node
        full_free_nodes.append(bnd_name)
        TT.append(ctr(bnd_name, epoch_bnd, x_bnd_node,
                      dynamics='ES_gravity', coi=coi_bnd))

        # Impulsive maneuver at the boundary (DSM1 / DSM2): brings the arrival
        # velocity up to the next arc's departure velocity (x_bnd = next start).
        if man_after is not None:
            dsm_vect = x_bnd_phys[3:] - arrival_state[3:]
            TT.append(man(man_after, bnd_name, 0, dsm_vect))

        prev_node       = bnd_name
        epoch_arc_start = epoch_bnd

    # ---- Arc D (DSM2 -> comet): nD nodes, last = comet, dt = t / nD ----------
    # Heliocentric cruise to the comet: Sun-centric integration (coi=Sun).
    # coi_arcD = 'Sun'
    coi_arcD = 'Earth'
    arc_days = t_dsm2_comet * T_CONV / 86400.0
    n_pt_arc = max(1, int(round(arc_days / nodes_vec[3])))   # arc-D nodes (incl. comet)
    t_arc_s  = t_dsm2_comet * T_CONV
    dt_s     = t_arc_s / n_pt_arc
    x_arc_phys  = _adim_to_phys(x_dsm2, epoch_arc_start)
    node_epochs = [epoch_arc_start + i * dt_s for i in range(1, n_pt_arc + 1)]
    arc_states  = _propagate_es(x_arc_phys, epoch_arc_start, node_epochs)

    # Leading match: ctr_dsm2a -> first arc-D node (arc D's dt; through man_dsm2).
    TT.append(match('match_' + prev_node,
                    prev_node, dt_s / 2,
                    'ctr' + str(gidx + 1) + 'a', dt_s / 2))

    for i in range(1, n_pt_arc + 1):
        gidx += 1
        name  = 'ctr' + str(gidx) + 'a'
        x_loc = arc_states[i - 1]

        ref_state[name] = x_loc
        TT.append(ctr(name, node_epochs[i - 1], x_loc,
                      dynamics='ES_gravity', coi=coi_arcD))

        if i < n_pt_arc:
            # Intermediate node (comet excluded): pos+vel free, optional maneuver.
            full_free_nodes.append(name)
            if with_man:
                man_names.append('man' + str(gidx) + 'a')
                TT.append(man('man' + str(gidx) + 'a', name, 0,
                              np.array([0.00001, 0.00001, 0.00001])))
            TT.append(match('match' + str(gidx) + 'a', name, dt_s / 2,
                            'ctr' + str(gidx + 1) + 'a', dt_s / 2))

    comet_name = 'ctr' + str(gidx) + 'a'   # last arc-D node (comet intercept)

    # -------------------------------------------------------------------------
    # STEP 13: End point (comet intercept)
    # -------------------------------------------------------------------------

    end_point = {
        'type': 'point',
        'input': 'SC',
        'name': 'end_pointa',
        'point': {'reference': comet_name, 'dt': '0 s'}
    }
    # dt='0 s': end_pointa coincides exactly with the last arc-D node (comet).
    TT.append(end_point)

    config_traj['timeline'] = TT
    # Full timeline: ctr0a + man_injection +
    #   [arc A internal + match] + ctr_dsm1a + man_dsm1 +
    #   [arc B internal + match] + ctr_flybya (no man) +
    #   [arc C internal + match] + ctr_dsm2a + man_dsm2 +
    #   [arc D nodes + match] (last = comet) + end_pointa

    # -------------------------------------------------------------------------
    # STEP 14: Define free variables for both problems
    # -------------------------------------------------------------------------
    # Epochs are kept FIXED in both phases (no '_dt' variables): this shrinks the
    # search space and keeps the mission timeline locked to the current dates.
    # Free variables:
    #   - position + velocity of every intermediate node (internal nodes + the
    #     explicit DSM1 / flyby / DSM2 junction nodes), collected in full_free_nodes
    #   - velocity only of the comet arrival node (position = fixed target)
    #   - injection / DSM1 / DSM2 impulsive ΔV
    # FF (Phase A, feasibility) and FF2 (Phase B, min TotalDV) are now identical.

    # Injection ΔV free in both (timing fixed).
    FF.append('man_injection_dv_x')
    FF.append('man_injection_dv_y')
    FF.append('man_injection_dv_z')
    FF2.append('man_injection_dv_x')
    FF2.append('man_injection_dv_y')
    FF2.append('man_injection_dv_z')

    # Fully-free nodes (internal + DSM1/flyby/DSM2 boundaries): position + velocity
    # (+ SC_dv). Epoch fixed → no 'ctr{i}a_dt'.
    for name in full_free_nodes:
        for lista in [FF, FF2]:
            lista.append(name + '_SC_dv')
            lista.append(name + '_SC_center_pos_x')
            lista.append(name + '_SC_center_pos_y')
            lista.append(name + '_SC_center_pos_z')
            lista.append(name + '_SC_center_vel_x')
            lista.append(name + '_SC_center_vel_y')
            lista.append(name + '_SC_center_vel_z')

    # DSM1 and DSM2 ΔV components free in both problems (timing fixed).
    FF.append('man_dsm1_dv_x')
    FF.append('man_dsm1_dv_y')
    FF.append('man_dsm1_dv_z')
    FF2.append('man_dsm1_dv_x')
    FF2.append('man_dsm1_dv_y')
    FF2.append('man_dsm1_dv_z')

    FF.append('man_dsm2_dv_x')
    FF.append('man_dsm2_dv_y')
    FF.append('man_dsm2_dv_z')
    FF2.append('man_dsm2_dv_x')
    FF2.append('man_dsm2_dv_y')
    FF2.append('man_dsm2_dv_z')

    if with_man:
        for mname in man_names:
            FF.append(mname + '_dv_x')
            FF.append(mname + '_dv_y')
            FF.append(mname + '_dv_z')
            FF2.append(mname + '_dv_x')
            FF2.append(mname + '_dv_y')
            FF2.append(mname + '_dv_z')

    # Last node (comet): only velocity is free — position is fixed (comet target).
    FF.append(comet_name + '_SC_dv')
    FF.append(comet_name + '_SC_center_vel_x')
    FF.append(comet_name + '_SC_center_vel_y')
    FF.append(comet_name + '_SC_center_vel_z')
    FF2.append(comet_name + '_SC_dv')
    FF2.append(comet_name + '_SC_center_vel_x')
    FF2.append(comet_name + '_SC_center_vel_y')
    FF2.append(comet_name + '_SC_center_vel_z')


    # -------------------------------------------------------------------------
    # STEP 15: Scales and bounds
    # -------------------------------------------------------------------------

    SS = {}
    BB = {}

    delta_pos       = 200     # refinement position bound half-width [km]
    delta_pos_flyby = 100     # tighter half-width for the lunar-flyby node [km]
    delta_vel       = 0.0001      # node velocity bound half-width [km/s]
    dv_man          = 1.5      # impulsive maneuver bound (±500 m/s) [km/s]

    # Injection ΔV: ±500 m/s (symmetric).
    for _c in ['x', 'y', 'z']:
        SS['man_injection_dv_' + _c] = '0.1 km/s'
        BB['man_injection_dv_' + _c] = [f'{-dv_man:.3f} km/s', f'{dv_man:.3f} km/s']

    # (Epochs fixed: no 'ctr0a_dt' / 'ctr{i}a_dt' variables.)

    # Fully-free nodes (internal + DSM1/flyby/DSM2 boundaries): position + velocity.
    for _pfx in full_free_nodes:
        SS[_pfx + '_SC_dv'] = '1 km/s'
        BB[_pfx + '_SC_dv'] = ['0 km/s', '1 km/s']

        _ref = ref_state[_pfx]

        # Tighter position bounds for the explicit lunar-flyby control point.
        if _pfx == 'ctr_flybya':
            _dp = delta_pos_flyby
        else:
            _dp = delta_pos

        for _j, _c in enumerate(['pos_x', 'pos_y', 'pos_z']):
            _pos_scale = max(abs(float(_ref[_j])), 1e5)   # min scale 100 000 km
            SS[_pfx + '_SC_center_' + _c] = f'{_pos_scale:.3e} km'
            BB[_pfx + '_SC_center_' + _c] = [
                str(float(_ref[_j]) - _dp) + ' km',
                str(float(_ref[_j]) + _dp) + ' km'
            ]

        # Velocity: bounds = guess ± delta_vel (relaxed).
        for _j, _c in enumerate(['vel_x', 'vel_y', 'vel_z']):
            _vel_val = float(_ref[_j + 3])
            _vel_scale = max(abs(_vel_val), 0.001)
            SS[_pfx + '_SC_center_' + _c] = f'{_vel_scale:.6f} km/s'
            BB[_pfx + '_SC_center_' + _c] = [
                f'{_vel_val - delta_vel:.6f} km/s',
                f'{_vel_val + delta_vel:.6f} km/s'
            ]

    # DSM1 and DSM2: ±500 m/s.
    for _dsm in ['man_dsm1', 'man_dsm2']:
        for _comp in ['x', 'y', 'z']:
            SS[_dsm + '_dv_' + _comp] = '0.1 km/s'
            BB[_dsm + '_dv_' + _comp] = [f'{-dv_man:.3f} km/s', f'{dv_man:.3f} km/s']

    # Optional intermediate maneuvers: ±150 m/s.
    if with_man:
        for _mname in man_names:
            for _comp in ['x', 'y', 'z']:
                SS[_mname + '_dv_' + _comp] = '0.1 km/s'
                BB[_mname + '_dv_' + _comp] = ['-0.15 km/s', '0.15 km/s']

    # Last node (comet): only velocity free, ±10 km/s (position is the fixed target).
    _pfx_last = comet_name
    SS[_pfx_last + '_SC_dv'] = '1 km/s'
    BB[_pfx_last + '_SC_dv'] = ['0 km/s', '1 km/s']
    _ref_last = ref_state[comet_name]
    for _j, _c in enumerate(['vel_x', 'vel_y', 'vel_z']):
        _vel_val = float(_ref_last[_j + 3])
        _vel_scale = max(abs(_vel_val), 0.001)
        SS[_pfx_last + '_SC_center_' + _c] = f'{_vel_scale:.6f} km/s'
        BB[_pfx_last + '_SC_center_' + _c] = [
            f'{_vel_val - delta_vel:.6f} km/s',
            f'{_vel_val + delta_vel:.6f} km/s'
        ]


    # -------------------------------------------------------------------------
    # STEP 16: Assemble and return
    # -------------------------------------------------------------------------

    config_prob['parameters']['free'] = FF
    config_prob2['parameters']['free'] = FF2
    config_prob['parameters']['scales'] = SS
    config_prob['parameters']['bounds'] = BB
    config_prob2['parameters']['scales'] = SS
    config_prob2['parameters']['bounds'] = BB
    # SS and BB are shared: GODOT only applies entries that appear in each
    # problem's 'free' list, so extra entries are silently ignored.

    return config_traj, config_prob, config_prob2


# =============================================================================
# COMPLETE FLOW SUMMARY
# =============================================================================
#
#  Physical constants (D, GM_SUN, GM_EARTH, mu_se, T_SUN_EARTH)
#       used everywhere to convert CR3BP non-dim ↔ physical km/km/s
#                           │
#              ┌────────────┴────────────┐
#              ▼                         ▼
#       config_halo()             config_trajectory()
#              │                         │
#   [CR3BP → ephemeris coords]   [x0_halo + CR3BP guesses → ephemeris]
#              │                         │
#   [timeline: ctr+man+match]   [timeline: inj + 4 arcs + 2 DSMs]
#              │                         │
#   [problem: min TotalDV_halo]  [prob:  feasibility only   ]
#              │                  [prob2: min TotalDV_traj   ]
#              ▼                         ▼
#        (config_traj,           (config_traj,
#         config_prob)            config_prob, config_prob2)
#              │                         │
#           main.py Section 7      main.py Section 9
# =============================================================================


