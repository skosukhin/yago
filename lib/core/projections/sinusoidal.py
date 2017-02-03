import numpy as np

from core.projections.projection import Projection
from core.projections.rotors import Rotor, RotorZ, RotorY


class SinusoidalProjection(Projection):
    """
    Links the geographical coordinate system with a Cartesian coordinate system
    defined on a Sinusoidal projection plane. The intersection of the equator
    and the Greenwich Meridian is always in the center of the projection. The
    X-axis is aligned from the center to the eastern hemisphere along the
    equator and the Y-axis is aligned from the same origin to the northern
    hemisphere along the Greenwich Meridian.
    """

    short_name = 'sinus'
    long_name = 'Sinusoidal'
    standard_name = 'sinusoidal'

    _TRUE_SCALE_LAT = 0.0

    def __init__(self, earth_radius):
        """
        Constructor of the class.
        :param earth_radius: Earth radius (in meters).
        """
        self.earth_radius = earth_radius
        self.true_scale_lats = [SinusoidalProjection._TRUE_SCALE_LAT]

    @classmethod
    def unified_init(cls, earth_radius, true_lats):
        if len(true_lats) == 1:
            eps = np.finfo(SinusoidalProjection._TRUE_SCALE_LAT).eps
            diff = np.fabs(true_lats[0] - SinusoidalProjection._TRUE_SCALE_LAT)
            if diff >= eps:
                raise Exception('True scale for Sinusoidal projection must '
                                'be equal to ' +
                                str(SinusoidalProjection._TRUE_SCALE_LAT) +
                                '.')

        elif len(true_lats) != 0:
            raise Exception('The list of true scales for Sinusoidal '
                            'projection must either contain exactly one '
                            'value or be empty.')

        return SinusoidalProjection(earth_radius)

    @classmethod
    def build_rotor(cls, orig_lat, orig_lon, add_angle_deg):
        return Rotor.chain(RotorZ(180.0 - orig_lon),
                           RotorY(90.0 - orig_lat),
                           RotorZ(add_angle_deg),
                           RotorY(90.0))

    def convert_points(self, lats, lons):
        lats_rad = np.radians(lats)
        yy = self.earth_radius * lats_rad
        xx = self.earth_radius * np.cos(lats_rad) * np.radians(lons)
        return xx, yy

    def restore_points(self, xx, yy):
        lats_rad = yy / self.earth_radius
        lons = np.degrees(xx / (self.earth_radius * np.cos(lats_rad)))
        return np.degrees(lats_rad), lons

    def convert_vectors(self, uu, vv, lats, lons, return_points=False):
        if return_points:
            xx, yy = self.convert_points(lats, lons)
            return uu, vv, xx, yy
        else:
            return uu, vv

    def restore_vectors(self, uu, vv, xx, yy, return_points=False):
        if return_points:
            lats, lons = self.restore_points(xx, yy)
            return uu, vv, lats, lons
        else:
            return uu, vv
