import numpy as np

from core.common import HALF_PI
from core.projections.projection import Projection


class MercatorProjection(Projection):
    """
    Links the geographical coordinate system with a Cartesian coordinate system
    defined on a Mercator projection plane. The intersection of the equator and
    the Greenwich Meridian is always in the center of the projection. The
    X-axis is aligned from the center to the eastern hemisphere along the
    equator and the Y-axis is aligned from the same origin to the northern
    hemisphere along the Greenwich Meridian.
    """

    short_name = 'mercator'
    long_name = 'Mercator'
    standard_name = 'mercator'

    def __init__(self, true_lat, earth_radius):
        """
        Constructor of the class.
        :param true_lat: Latitude of true scale (in degrees).
        :param earth_radius: Earth radius (in meters).
        """
        self.true_scale_lats = [true_lat]
        self.earth_radius = earth_radius

        self._k = np.cos(np.radians(true_lat))

    @classmethod
    def unified_init(cls, earth_radius, true_lats):
        if len(true_lats) == 0:
            true_lats = [0.0]
        elif len(true_lats) != 1:
            raise Exception('The list of true scales for Mercator projection '
                            'must either contain exactly one '
                            'value or be empty.')

        return MercatorProjection(true_lats[0], earth_radius)

    @property
    def reference_point(self):
        return 0.0, 0.0

    def convert_points(self, lats, lons):
        xx = self._k * self.earth_radius * np.radians(lons)
        yy = self._k * self.earth_radius * np.log(
            np.tan((HALF_PI + np.radians(lats)) / 2.0))
        return xx, yy

    def restore_points(self, xx, yy):
        xx, yy = np.asanyarray(xx), np.asanyarray(yy)
        lons = np.degrees(xx / (self.earth_radius * self._k))
        lats = np.degrees(2.0 * np.arctan(
            np.exp(yy / (self.earth_radius * self._k))) - HALF_PI)
        return lats, lons

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

    def get_scale_factors(self, lats, lons):
        k = self._k * np.cos(np.radians(lats))
        return k, k
