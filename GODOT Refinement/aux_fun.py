# =============================================================================
# AUX_FUN.PY — GODOT timeline-element builders
# =============================================================================
#
# Helper functions that assemble the dictionaries GODOT expects for the three
# building blocks of a multiple-shooting timeline:
#   ctr()   -> control point (a state to be propagated)
#   man()   -> impulsive manoeuvre (a ΔV applied at a reference event)
#   match() -> match/continuity constraint (two arcs must meet)
#
# They exist only to avoid repeating verbose nested dictionaries in the driver
# scripts (main.py, halo_ephe.py, transfer_ephe.py).
# =============================================================================


def ctr(name , date , x0_loc , dynamics='EMS_gravity', coi='Earth') :
    """Build a GODOT 'control' point.

    A control point pins the spacecraft state at a given epoch; the optimizer is
    allowed to move it (within bounds set elsewhere) and GODOT propagates from it.

    Inputs:
        name     : unique event name (e.g. 'ctr0', 'ctr_flyby').
        date     : epoch of the control point (tempo.Epoch or string).
        x0_loc   : Earth-centered SEROT state [x, y, z, vx, vy, vz] (km, km/s).
        dynamics : force model used to propagate from this point
                   ('EMS_gravity' = Earth+Moon+Sun, 'ES_gravity' = Earth+Sun).
        coi      : centre of integration ('Earth', 'Moon' or 'Sun').

    Output:
        dict : GODOT 'control' timeline element.
    """
    # 'point' is always Earth: the 'value' state (and hence the bounds and the
    # match comparison, which happens in the Earth frame) is Earth-centered SEROT.
    # 'coi' (centre of integration) is parametrized: switching it along the
    # trajectory (Earth -> Moon at the flyby -> Sun in cruise) improves numerical
    # conditioning without changing how the state is represented.
    ctr = {
        'type': 'control',
        'epoch': str(date),
        'name': name,
        'state': [
            {'name': 'SC_center',
            'point': 'Earth',
            'coi': coi,
            'axes': 'SEROT',
            'project': False,
            'dynamics': dynamics,
            'value': {
                'pos_x': str(x0_loc[0]) + ' km',
                'pos_y': str(x0_loc[1]) + ' km',
                'pos_z': str(x0_loc[2]) + ' km',
                'vel_x': str(x0_loc[3]) + ' km/s',
                'vel_y': str(x0_loc[4]) + ' km/s',
                'vel_z': str(x0_loc[5]) + ' km/s'
            }},
            {'name': 'SC_dv',
             'value': '0 m/s'}
        ]}
    return ctr

def man(name , reference , dt_s , dv_vect_kms) :
    """Build a GODOT impulsive 'manoeuvre'.

    Applies an instantaneous ΔV to the spacecraft at a fixed time offset from a
    reference event, expressed in the SEROT axes.

    Inputs:
        name        : unique manoeuvre name (e.g. 'man_dsm1', 'man0').
        reference   : name of the event the burn is anchored to.
        dt_s        : time offset from the reference event [s].
        dv_vect_kms : ΔV vector [dv_x, dv_y, dv_z] in SEROT axes [km/s].

    Output:
        dict : GODOT 'manoeuvre' timeline element.
    """
    man = {'type': 'manoeuvre',
            'name': name ,
            'model': 'impulsive',
            'input': 'SC',
            'thruster' : 'main',
            'config': { 'point': {'reference': reference, 'dt' : str(dt_s) + ' s'} ,
                'direction': {
                           'axes': 'SEROT',
                           'dv_x': str(dv_vect_kms[0]) +' km/s' ,
                           'dv_y': str(dv_vect_kms[1]) +' km/s' ,
                           'dv_z': str(dv_vect_kms[2]) +' km/s'
                       }}}
    return man

def match(name , ref_left , dt_left , ref_right , dt_right) :
    """Build a GODOT 'match' (continuity) constraint.

    Requires the two arcs propagated from a left and a right reference event to
    coincide (Cartesian position and velocity) at their meeting point — the core
    constraint of a multiple-shooting scheme.

    Inputs:
        name      : unique constraint name (e.g. 'match0').
        ref_left  : name of the left reference event.
        dt_left   : time offset from the left reference to the match point [s].
        ref_right : name of the right reference event.
        dt_right  : time offset from the right reference to the match point [s].

    Output:
        dict : GODOT 'match' timeline element (Cartesian, in the Earth frame).
    """
    match = {'type': 'match' ,
             'name': name ,
             'input': 'SC' ,
             'left': {
                 'reference': ref_left,
                 'dt': str(dt_left) + ' s'},
             'right': {
                 'reference': ref_right,
                 'dt': str(dt_right) + ' s'},
             'body': 'Earth',
             'vars': 'cart'}
    return match
    