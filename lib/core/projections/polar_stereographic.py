import numpy as np

from core.common import cos_sin_deg, HALF_PI, rotate_vectors
from core.projections.projection import Projection
from core.projections.rotors import Rotor, RotorZ, RotorY


class PolarStereographicProjection(Projection):
    """
    Links the geographical coordinate system with a Cartesian coordinate system
    defined on a Polar stereographic projection plane. The North Pole is always
    in the center of the projection. The X-axis is aligned from the North Pole
    to the South Pole over the 90th meridian and the Y-axis is aligned from the
    same origin to the South Pole over the 180th meridian.
    """

    short_name = 'stereo'
    long_name = 'Polar stereographic'
    standard_name = 'polar_stereographic'

    def __init__(self, true_lat, earth_radius):
        """
        Constructor of the class.
        :param true_lat: Latitude of true scale (in degrees).
        :param earth_radius: Earth radius (in meters).
        """
        self.z = np.sin(np.radians(true_lat)) + 1.0
        self.true_scale_lats = [true_lat]
        self.earth_radius = earth_radius

    @classmethod
    def unified_init(cls, earth_radius, true_lats):
        if len(true_lats) == 0:
            true_lats = [90.0]
        elif len(true_lats) != 1:
            raise Exception('The list of true scales for Polar stereographic '
                            'projection must either contain exactly one '
                            'value or be empty.')

        return PolarStereographicProjection(true_lats[0], earth_radius)

    @classmethod
    def build_rotor(cls, center_lat, center_lon, z_angle):
        return Rotor.chain(RotorZ(180.0 - center_lon),
                           RotorY(90.0 - center_lat),
                           RotorZ(z_angle))

    def convert_points(self, lats, lons):
        c_lats, s_lats = cos_sin_deg(lats)
        c_lons, s_lons = cos_sin_deg(lons)
        rr = self.z * self.earth_radius * c_lats / (1.0 + s_lats)
        xx = rr * s_lons
        yy = -rr * c_lons
        return xx, yy

    def restore_points(self, xx, yy):
        rr = np.sqrt(xx * xx + yy * yy) / self.earth_radius
        lats = np.degrees(HALF_PI - 2.0 * np.arctan(rr / self.z))
        lons = np.degrees(np.arctan2(xx, -yy))
        return lats, lons

    def convert_vectors(self, uu, vv, lats, lons, return_points=False):
        rot_uu, rot_vv = rotate_vectors(uu, vv, lons)
        if return_points:
            xx, yy = self.convert_points(lats, lons)
            return rot_uu, rot_vv, xx, yy
        else:
            return rot_uu, rot_vv

    def restore_vectors(self, uu, vv, xx, yy, return_points=False):
        lats, lons = self.restore_points(xx, yy)
        c_lons, s_lons = cos_sin_deg(-lons)
        rot_matrices = np.asanyarray([[c_lons, -s_lons], [s_lons, c_lons]])
        rot_vecs = np.einsum('ij...,j...', rot_matrices, np.stack([uu, vv]))
        rot_uu = rot_vecs[..., 0]
        rot_vv = rot_vecs[..., 1]
        if return_points:
            return rot_uu, rot_vv, lats, lons
        else:
            return rot_uu, rot_vv
