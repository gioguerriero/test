# =============================================================================
# TOTALDV.PY — Objective and constraint functions for the GODOT optimization
# =============================================================================
#
# Defines the ScalarTimeEvaluable objects that GODOT/PyGMO minimize or constrain.
# A ScalarTimeEvaluable exposes an eval(epoch) method returning a single scalar,
# built with autodif (ad) operations so the optimizer gets exact gradients.
#
# Classes:
#   totalDV_traj        - total transfer ΔV (injection + DSM1 + DSM2)
#   totalDV_halo        - total station-keeping ΔV of an n-manoeuvre Halo orbit
#   flybyMoonAltitude   - lunar-altitude inequality constraint at the flyby node
#   totalDV_traj_geo    - single-DSM transfer ΔV (geocentric variant)
# =============================================================================

import godot.model.common as common
from godot.cosmos import util
import godot.core.util as c_util
from godot.core import autodif as ad
import numpy as np

c_util.suppressLogger()   # keep the console clean of GODOT's internal logging

class totalDV_traj( common.ScalarTimeEvaluable ) :
    """Total transfer ΔV: injection + DSM1 + DSM2 manoeuvres [km/s]."""

    def __init__(self, universe):
        super().__init__()
        self.universe = universe

    def eval(self, e ):
        # Sum the magnitudes of the three transfer manoeuvres. NB: dv_tot is
        # reset to 0.0 just before the return (see below), so this objective
        # currently evaluates to zero — kept as-is to preserve existing logic.

        dv_tot = 0.0

        dvx_inj = self.universe.evaluables.get(f'man_injection_dv_x').eval(e)
        dvy_inj = self.universe.evaluables.get(f'man_injection_dv_y').eval(e)
        dvz_inj = self.universe.evaluables.get(f'man_injection_dv_z').eval(e)
        dv_tot += ad.sqrt(dvx_inj*dvx_inj + dvy_inj*dvy_inj + dvz_inj*dvz_inj)

        dvx_dsm = self.universe.evaluables.get(f'man_dsm1_dv_x').eval(e)
        dvy_dsm = self.universe.evaluables.get(f'man_dsm1_dv_y').eval(e)
        dvz_dsm = self.universe.evaluables.get(f'man_dsm1_dv_z').eval(e)
        dv_tot += ad.sqrt(dvx_dsm*dvx_dsm + dvy_dsm*dvy_dsm + dvz_dsm*dvz_dsm)

        dvx_dsm = self.universe.evaluables.get(f'man_dsm2_dv_x').eval(e)
        dvy_dsm = self.universe.evaluables.get(f'man_dsm2_dv_y').eval(e)
        dvz_dsm = self.universe.evaluables.get(f'man_dsm2_dv_z').eval(e)
        dv_tot += ad.sqrt(dvx_dsm*dvx_dsm + dvy_dsm*dvy_dsm + dvz_dsm*dvz_dsm)

        dv_tot = 0.0   # objective overridden to zero (kept unchanged)

        return dv_tot


class totalDV_halo( common.ScalarTimeEvaluable ) :
    """Total station-keeping ΔV of a Halo orbit with n_man manoeuvres [km/s].

    Inputs:
        n_man    : number of impulsive manoeuvres (named man0 … man{n_man-1}).
        universe : GODOT Universe holding the per-component ΔV evaluables.
    """

    def __init__(self, n_man , universe):
        super().__init__()
        self.n_man = n_man
        self.universe = universe

    def eval(self, e ):
        # Sum |ΔV| over all n_man manoeuvres of the Halo timeline.

        dv_tot = 0.0
        for i__ in range(self.n_man) :
            dvx = self.universe.evaluables.get(f'man{i__}_dv_x').eval(e)
            dvy = self.universe.evaluables.get(f'man{i__}_dv_y').eval(e)
            dvz = self.universe.evaluables.get(f'man{i__}_dv_z').eval(e)
            dv_tot += ad.sqrt(dvx*dvx + dvy*dvy + dvz*dvz)

        return dv_tot



class flybyMoonAltitude( common.ScalarTimeEvaluable ) :
    """Altitude of the spacecraft above the lunar surface at epoch e [km].

    Used as an inequality constraint at the flyby node (lunar pericentre) to
    stop the optimizer from driving the trajectory into the Moon (which would
    trigger the BodyCache evalAcc "distance is zero" failure).

    Inputs:
        universe  : GODOT Universe (queried for the Moon->SC vector).
        r_moon_km : lunar mean radius [km] (default 1737.4).
    """
    def __init__(self, universe, r_moon_km=1737.4):
        super().__init__()
        self.universe = universe
        self.r_moon   = r_moon_km

    def eval(self, e):
        # Moon->SC vector in the inertial ICRF frame; its magnitude minus the
        # lunar radius is the surface altitude (>0 outside the Moon).
        r = self.universe.frames.vector3('Moon', 'SC_center', 'ICRF', e)
        dist = ad.sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2])
        return dist - self.r_moon


class totalDV_traj_geo( common.ScalarTimeEvaluable ) :
    """Transfer ΔV for the single-DSM geocentric variant (|man_dsm|) [km/s]."""

    def __init__(self, universe):
        super().__init__()
        self.universe = universe
# SC_dv
    def eval(self, e ):
        # Magnitude of the single deep-space manoeuvre 'man_dsm'.
        dv_tot = 0.0

        dvx_dsm = self.universe.evaluables.get(f'man_dsm_dv_x').eval(e)
        dvy_dsm = self.universe.evaluables.get(f'man_dsm_dv_y').eval(e)
        dvz_dsm = self.universe.evaluables.get(f'man_dsm_dv_z').eval(e)
        dv_tot = ad.sqrt(dvx_dsm*dvx_dsm + dvy_dsm*dvy_dsm + dvz_dsm*dvz_dsm)
        return dv_tot
