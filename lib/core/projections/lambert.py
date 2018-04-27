import numpy as np

from core.common import HALF_PI, rotate_vectors_deg
from core.projections.projection import Projection
from core.projections.rotors import Rotor, RotorZ, RotorY


class LambertConformalProjection(Projection):
    """
    Links the geographical coordinate system with a Cartesian coordinate system
    defined on a Lambert conformal conic projection plane. The intersection of
    the 45th northern parallel and the Greenwich Meridian is always in the
    center of the projection. The X-axis is aligned from the center to the
    eastern hemisphere along the equator and the Y-axis is aligned from the
    same origin to the northern hemisphere along the Greenwich Meridian (both
    of the directions are true only in the vicinity of the center point).
    """

    _CENTER_LAT = 45.0
    short_name = 'lambert'
    long_name = 'Lambert conformal conic'
    standard_name = 'lambert_conformal_conic'

    def __init__(self, true_lat1, true_lat2, earth_radius):
        """
        Constructor of the class.
        :param true_lat1: Southern latitude of true scale (in degrees, must be
        less or equal to 45).
        :param true_lat2: Northern latitude of true scale (in degrees, must be
        greater or equal to 45).
        :param earth_radius: Earth radius (in meters).
        """
        self.earth_radius = earth_radius
        self.true_scale_lats = [true_lat1, true_lat2]

        phi1_rad = np.radians(true_lat1)
        phi2_rad = np.radians(true_lat2)

        self._cos_phi1 = np.cos(phi1_rad)
        self._tan_quart_pi_half_phi1 = np.tan((HALF_PI + phi1_rad) / 2.0)

        if np.fabs(true_lat1 - true_lat2) < np.finfo(float).eps:
            self._n = np.sin(phi1_rad)
        else:
            self._n = (np.log(self._cos_phi1 / np.cos(phi2_rad)) /
                       np.log(
                           np.tan((HALF_PI + phi2_rad) / 2.0) /
                           self._tan_quart_pi_half_phi1))

        self._f = (self._cos_phi1 *
                   np.power(self._tan_quart_pi_half_phi1, self._n) / self._n)

        tan_quart_pi_half_center_lat = np.tan(
            (HALF_PI + np.radians(LambertConformalProjection._CENTER_LAT)) /
            2.0)

        self._rho0 = (self._f *
                      np.power(tan_quart_pi_half_center_lat, -self._n))

    @classmethod
    def unified_init(cls, earth_radius, true_lats):
        if len(true_lats) == 0:
            true_lats = [LambertConformalProjection._CENTER_LAT,
                         LambertConformalProjection._CENTER_LAT]
        elif len(true_lats) != 2:
            raise Exception('The list of true scales for Lambert conformal '
                            'projection must either contain exactly two '
                            'values or be empty.')

        return LambertConformalProjection(true_lats[0], true_lats[1],
                                          earth_radius)

    @classmethod
    def build_rotor(cls, orig_lat, orig_lon, add_angle_deg):
        return Rotor.chain(RotorZ(180.0 - orig_lon),
                           RotorY(90.0 - orig_lat),
                           RotorZ(add_angle_deg),
                           RotorY(45.0))

    def convert_points(self, lats, lons):
        xx_unitless, yy_unitless = self._convert_to_unitless(lats, lons)
        return xx_unitless * self.earth_radius, yy_unitless * self.earth_radius

    def restore_points(self, xx, yy):
        xx, yy = np.asanyarray(xx), np.asanyarray(yy)
        xx_unitless = xx / self.earth_radius
        yy_unitless = yy / self.earth_radius
        return self._restore_from_unitless(xx_unitless, yy_unitless)

    def convert_vectors(self, uu, vv, lats, lons, return_points=False):
        rot_uu, rot_vv = rotate_vectors_deg(uu, vv, lons * self._n)
        if return_points:
            xx, yy = self.convert_points(lats, lons)
            return rot_uu, rot_vv, xx, yy
        else:
            return rot_uu, rot_vv

    def restore_vectors(self, uu, vv, xx, yy, return_points=False):
        lats, lons = self.restore_points(xx, yy)
        rot_uu, rot_vv = rotate_vectors_deg(uu, vv, -lons * self._n)
        if return_points:
            return rot_uu, rot_vv, lats, lons
        else:
            return rot_uu, rot_vv

    def get_scale_factors(self, lats, lons):
        lats_rad = np.radians(lats)

        k = (self._cos_phi1 / np.cos(lats_rad) *
             np.power(self._tan_quart_pi_half_phi1 /
                      np.tan((HALF_PI + lats_rad) / 2.0),
                      self._n))

        return k, k

    def _convert_to_unitless(self, lats, lons):
        rho = (self._f *
               np.power(np.tan((HALF_PI + np.radians(lats)) / 2.0), -self._n))

        nlo = self._n * np.radians(lons)
        x = rho * np.sin(nlo)
        y = self._rho0 - rho * np.cos(nlo)
        return x, y

    def _restore_from_unitless(self, xx, yy):
        rho = np.sign(self._n) * np.sqrt(np.power(xx, 2.0) +
                                         np.power(self._rho0 - yy, 2.0))

        theta = np.arctan(xx / (self._rho0 - yy))

        lats = (2.0 * np.arctan(np.power(self._f / rho, 1.0 / self._n)) -
                HALF_PI)

        lons = theta / self._n
        return np.degrees(lats), np.degrees(lons)
