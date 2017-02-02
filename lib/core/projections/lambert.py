import numpy as np

from core.common import build_2d_rotation_z_rad, QUARTER_PI
from core.projections.projection import Projection
from core.rotors import Rotor, RotorZ, RotorY


class LambertConformalProjection(Projection):
    """
    Links spherical lat/lon coordinates and the corresponding vectors with
    Cartesian coordinates onto a Lambert conformal conic projection plane. The
    intersection of the 45th northern parallel and the Greenwich Meridian is
    always in the center of the projection. The X-axis is aligned from the
    center to the eastern hemisphere along the equator and the Y-axis is
    aligned from the same origin to the northern hemisphere along the Greenwich
    Meridian (both of the given directions are true only in the vicinity of the
    point of origin).
    """

    _CENTER_LAT = 45.0
    short_name = 'lambert'
    long_name = 'Lambert conformal conic'
    standard_name = 'lambert_conformal_conic'

    def __init__(self, true_lat1, true_lat2, earth_radius):
        """
        The constructor of the class.
        :param true_lat1: The southern latitude (in degrees) of true scale
        (must be less or equal to 45 degrees)
        :param true_lat2: The northern latitude (in degrees) of true scale
        (must be greater or equal to 45 degrees)
        :param earth_radius: Earth radius (in meters).
        """
        phi1_rad = np.radians(true_lat1)
        phi2_rad = np.radians(true_lat2)
        self.earth_radius = earth_radius

        if np.fabs(true_lat1 - true_lat2) < np.finfo(float).eps:
            self.n = np.sin(phi1_rad)
        else:
            self.n = \
                (np.log(np.cos(phi1_rad) / np.cos(phi2_rad)) /
                 (np.log(np.tan(QUARTER_PI + phi2_rad / 2.0) /
                         np.tan(QUARTER_PI + phi1_rad / 2.0))))
        self.n_sign = np.sign(self.n)
        self.f = \
            (np.cos(phi1_rad) *
             np.power(np.tan(QUARTER_PI + phi1_rad / 2.0), self.n) / self.n)
        self.rho0 = \
            (self.f *
             np.power(1.0 /
                      np.tan(QUARTER_PI +
                             np.radians(
                                 LambertConformalProjection._CENTER_LAT) /
                             2.0),
                      self.n))

    @classmethod
    def init(cls, earth_radius, true_lats):
        return LambertConformalProjection(true_lats[0], true_lats[1],
                                          earth_radius)

    def build_rotor(self, orig_lat, orig_lon, add_angle_deg):
        """
        The function generates an instance of class Rotor to be used in
        conjunction with Lambert projection.
        :param orig_lat: Latitude (in degrees) of the origin point of the
        projection.
        :param orig_lon: Latitude (in degrees) of the origin point of the
        projection.
        :param add_angle_deg: Angle (in degrees) of the optional rotation
        around Z-axis.
        :return: Returns an instance of class Rotor that rotates the
        geographical coordinate system to move the origin point to the
        intersection of the 45th northern parallel and the Greenwich Meridian.
        """

        return Rotor.chain(RotorZ(180.0 - orig_lon), RotorY(90.0 - orig_lat),
                           RotorZ(add_angle_deg), RotorY(45.0))

    def convert_point(self, la, lo):
        """
        Calculates the Cartesian coordinates of the given point.
        :param la: Latitude of the point (in degrees).
        :param lo: Longitude of the point (in degrees).
        :return: A tuple (x, y) of the Cartesian coordinates of the given point
        on the projection plane.
        """
        x_unit, y_unit = self._convert_point_unit(la, lo)
        return x_unit * self.earth_radius, y_unit * self.earth_radius

    def restore_point(self, x, y):
        """
        Calculates the spherical coordinates of the given point.
        :param x: Coordinate along the X-axis.
        :param y: Coordinate along the Y-axis.
        :return: A tuple (lat, lon) of the spherical coordinates of the given
        point.
        """
        x_unit = x / self.earth_radius
        y_unit = y / self.earth_radius
        return self._restore_point_unit(x_unit, y_unit)

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
        rot_angle = np.radians(lo) * self.n
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
        rot_angle = -np.radians(lo) * self.n
        rot_matrix = build_2d_rotation_z_rad(rot_angle)
        rot_vec = np.dot(rot_matrix, np.array([u, v]))
        rot_u = rot_vec[0]
        rot_v = rot_vec[1]
        if return_point:
            return rot_u, rot_v, la, lo
        else:
            return rot_u, rot_v

    def _convert_point_unit(self, la, lo):
        rho = \
            self.f * \
            np.power(1.0 / np.tan(QUARTER_PI + np.radians(la) / 2.0), self.n)
        nlo = self.n * np.radians(lo)
        x = rho * np.sin(nlo)
        y = self.rho0 - rho * np.cos(nlo)
        return x, y

    def _restore_point_unit(self, x, y):
        rho = \
            self.n_sign * \
            np.sqrt(np.power(x, 2.0) + np.power(self.rho0 - y, 2.0))
        theta = np.arctan(x / (self.rho0 - y))
        la = \
            2.0 * np.arctan(np.power(self.f / rho, 1.0 / self.n)) - np.pi / 2.0
        lo = theta / self.n
        return np.degrees(la), np.degrees(lo)
