import numpy as np

from core.common import QUARTER_PI, HALF_PI
from core.projections.projection import Projection
from core.rotors import Rotor, RotorZ, RotorY


class MercatorProjection(Projection):
    """
    Links spherical lat/lon coordinates and the corresponding vectors with
    Cartesian coordinates onto a Mercator projection plane. The intersection
    of the equator and the Greenwich Meridian is always in the center of the
    projection. The X-axis is aligned from the center to the eastern hemisphere
    along the equator and the Y-axis is aligned from the same origin to the
    northern hemisphere along the Greenwich Meridian.
    """

    short_name = 'mercator'
    long_name = 'Mercator'
    standard_name = 'mercator'

    def __init__(self, true_lat, earth_radius):
        """
        The constructor of the class.
        :param true_lat: Latitude of true scale (in degrees).
        :param earth_radius: Earth radius (in meters).
        """
        self.k = np.cos(np.radians(true_lat))
        self.earth_radius = earth_radius

    @classmethod
    def init(cls, earth_radius, true_lats):
        return MercatorProjection(true_lats[0], earth_radius)

    def build_rotor(self, orig_lat, orig_lon, add_angle_deg):
        """
        The function generates an instance of class Rotor to be used in
        conjunction with Mercator projection.
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
        x = self.k * self.earth_radius * np.radians(lo)
        y = self.k * self.earth_radius * np.log(
            np.tan(QUARTER_PI + np.radians(la) / 2.0))
        return x, y

    def restore_point(self, x, y):
        """
        Calculates the spherical coordinates of the given point.
        :param x: Coordinate along the X-axis.
        :param y: Coordinate along the Y-axis.
        :return: A tuple (lat, lon) of the spherical coordinates of the given
        point.
        """
        lo = np.degrees(x / (self.earth_radius * self.k))
        la = np.degrees(2.0 * np.arctan(
            np.exp(y / (self.earth_radius * self.k))) - HALF_PI)
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
