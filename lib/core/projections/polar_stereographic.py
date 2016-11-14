import numpy as np

from core.common import build_2d_rotation_z_rad, cos_sin_deg


class PolarStereographicProjector(object):
    """
    Links spherical lat/lon coordinates and the corresponding vectors with
    Cartesian coordinates onto a Polar stereographic projection plane. The
    North Pole is always in the center of the projection. The X-axis is aligned
    from the North Pole to the South Pole over the 90th meridian and the Y-axis
    is aligned from the same origin to the South Pole over the 180th meridian.
    """
    half_pi = np.pi / 2.0
    short_name = 'stereo'
    long_name = 'Polar stereographic'
    standard_name = 'polar_stereographic'

    def __init__(self, true_lat=90.0, earth_radius=6370997.0):
        """
        The constructor of the class.
        :param true_lat: Latitude of true scale (in degrees).
        :param earth_radius: Earth radius (in meters).
        """
        self.z = np.sin(np.radians(true_lat)) + 1.0
        self.earth_radius = earth_radius

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
        la = np.degrees(self.half_pi - 2.0 * np.arctan(r / self.z))
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
