import numpy as np

from core.projections.projection import Projection
from core.rotors import Rotor, RotorZ, RotorY


class SinusoidalProjection(Projection):
    """
    Class that links spherical lat/lon coordinates and the corresponding
    vectors with the Cartesian coordinates onto the Sinusoidal projection
    plane. The intersection of the equator and the Greenwich Meridian is always
    in the center of the projection. The X-axis is aligned from the center to
    the eastern hemisphere along the equator and the Y-axis is aligned from the
    same origin to the northern hemisphere along the Greenwich Meridian.
    """

    short_name = 'sinusoidal'
    long_name = 'Sinusoidal'
    standard_name = 'sinusoidal'

    def __init__(self, earth_radius):
        """
        The constructor of the class.
        :param earth_radius: Earth radius (in meters).
        """
        self.earth_radius = earth_radius

    @classmethod
    def init(cls, earth_radius, true_lats):
        return SinusoidalProjection(earth_radius)

    def build_rotor(self, orig_lat, orig_lon, add_angle_deg):
        """
        The function generates an instance of class Rotor to be used in
        conjunction with Sinusoidal projection.
        :param orig_lat: Latitude (in degrees) of the origin point of the
        projection.
        :param orig_lon: Latitude (in degrees) of the origin point of the
        projection.
        :param add_angle_deg: Angle (in degrees) of the optional rotation
        around Z-axis.
        :return: Returns an instance of class Rotor that rotates the
        geographical coordinate system to move the origin point to the
        intersection of the equator and the Greenwich Meridian.
        """

        return Rotor.chain(RotorZ(180.0 - orig_lon), RotorY(90.0 - orig_lat),
                           RotorZ(add_angle_deg), RotorY(90.0))

    def convert_point(self, la, lo):
        """
        Calculates the Cartesian coordinates of the given point.
        :param la: Latitude of the point (in degrees).
        :param lo: Longitude of the point (in degrees).
        :return: A tuple (x, y) of the Cartesian coordinates of the given point
        on the projection plane.
        """
        la_rad = np.radians(la)
        y = self.earth_radius * la_rad
        x = self.earth_radius * np.cos(la_rad) * np.radians(lo)
        return x, y

    def restore_point(self, x, y):
        """
        Calculates the spherical coordinates of the given point.
        :param x: Coordinate along the X-axis.
        :param y: Coordinate along the Y-axis.
        :return: A tuple (lat, lon) of the spherical coordinates of the given
        point.
        """
        la_rad = y / self.earth_radius
        lo = np.degrees(x / (self.earth_radius * np.cos(la_rad)))
        return np.degrees(la_rad), lo

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
        if return_point:
            x, y = self.convert_point(la, lo)
            return u, v, x, y
        else:
            return u, v

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
        if return_point:
            la, lo = self.restore_point(x, y)
            return u, v, la, lo
        else:
            return u, v
