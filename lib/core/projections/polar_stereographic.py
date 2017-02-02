import numpy as np

from core.common import build_2d_rotation_z_rad, cos_sin_deg, HALF_PI
from core.projections.projection import Projection
from core.rotors import Rotor, RotorZ, RotorY


class PolarStereographicProjection(Projection):
    """
    Links spherical lat/lon coordinates and the corresponding vectors with
    Cartesian coordinates onto a Polar stereographic projection plane. The
    North Pole is always in the center of the projection. The X-axis is aligned
    from the North Pole to the South Pole over the 90th meridian and the Y-axis
    is aligned from the same origin to the South Pole over the 180th meridian.
    """

    short_name = 'stereo'
    long_name = 'Polar stereographic'
    standard_name = 'polar_stereographic'

    def __init__(self, true_lat, earth_radius):
        """
        The constructor of the class.
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

    def build_rotor(self, orig_lat, orig_lon, add_angle_deg):
        """
        The function generates an instance of class Rotor to be used in
        conjunction with Polar stereographic projection. The obtained
        instance is the result of three consecutive rotations: around Z-axis,
        around Y-axis, and again around Z-axis. The first two rotations shift
        the given point to the North Pole by rotating the coordinate system.
        The second rotation around Z-axis (the last one among the three) is
        optional to help users to adjust the orientation of the coordinate
        grid to account either for the features of the following projection
        procedure or for the plotting needs.
        :param orig_lat: Latitude (in degrees) of the origin point of the
        projection.
        :param orig_lon: Latitude (in degrees) of the origin point of the
        projection.
        :param add_angle_deg: Angle (in degrees) of the optional rotation
        around Z-axis.
        :return: Returns an instance of class Rotor that rotates the
        geographical coordinate system to move the origin point to the North
        Pole.
        """

        return Rotor.chain(RotorZ(180.0 - orig_lon), RotorY(90.0 - orig_lat),
                           RotorZ(add_angle_deg))

    def convert_point(self, la, lo):
        """
        Calculates the Cartesian coordinates of the given point.
        :param la: Latitude of the point (in degrees).
        :param lo: Longitude of the point (in degrees).
        :return: A tuple (x, y) of the Cartesian coordinates of the given point
        on the projection plane.
        """
        c_la, s_la = cos_sin_deg(la)
        c_lo, s_lo = cos_sin_deg(lo)
        r = self.z * self.earth_radius * c_la / (1.0 + s_la)
        x = r * s_lo
        y = -r * c_lo
        return x, y

    def restore_point(self, x, y):
        """
        Calculates the spherical coordinates of the given point.
        :param x: Coordinate along the X-axis.
        :param y: Coordinate along the Y-axis.
        :return: A tuple (lat, lon) of the spherical coordinates of the given
        point.
        """
        r = np.sqrt(x * x + y * y) / self.earth_radius
        la = np.degrees(HALF_PI - 2.0 * np.arctan(r / self.z))
        lo = np.degrees(np.arctan2(x, -y))
        return la, lo

    def convert_vector(self, u, v, la, lo, return_point=False):
        """
        Calculates the components along the X- and Y- axes of the given vector
        that originates from the given point.
        :param u: Zonal component of the vector.
        :param v: Meridional component of the vector.
        :param la: Latitude (in degrees) of the vector's origin.
        :param lo: Longitude (in degrees) of the vector's origin.
        :param return_point: Flag that tells the method to return the Cartesian
        coordinates of the vector's origin along with its components.
        :return: A tuple of the components of the vector along the X- and Y-
        axes respectively. If the flag return_point is True, than the result
        tuple is extended with the Cartesian (x, y) coordinates of the vector's
        origin.
        """
        rot_angle = np.radians(lo)
        rot_matrix = build_2d_rotation_z_rad(rot_angle)
        rot_vec = np.dot(rot_matrix, np.array([u, v]))
        rot_u = rot_vec[0]
        rot_v = rot_vec[1]
        if return_point:
            x, y = self.convert_point(la, lo)
            return rot_u, rot_v, x, y
        else:
            return rot_u, rot_v

    def restore_vector(self, u, v, x, y, return_point=False):
        """
        Calculates the zonal and meridional components of the given vector that
        originates from the given point.
        :param u: Component of the vector along the X-axis.
        :param v: Component of the vector along the Y-axis.
        :param x: Coordinate of the vector's origin along the X-axis.
        :param y: Coordinate of the vector's origin along the Y-axis.
        :param return_point: Flag that tells the method to return the spherical
        coordinates of the vector's origin along with its components.
        :return: A tuple of the zonal and meridional components of the vector
        respectively. If the flag return_point is True, than the result tuple
        is extended with the spherical (lat, lon) coordinates of the vector's
        origin.
        """
        la, lo = self.restore_point(x, y)
        rot_angle = -np.radians(lo)
        rot_matrix = build_2d_rotation_z_rad(rot_angle)
        rot_vec = np.dot(rot_matrix, np.array([u, v]))
        rot_u = rot_vec[0]
        rot_v = rot_vec[1]
        if return_point:
            return rot_u, rot_v, la, lo
        else:
            return rot_u, rot_v
