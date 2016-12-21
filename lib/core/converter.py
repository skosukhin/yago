import numpy as np


class Converter(object):
    """
    Class that combines the rotation and projection transformations.
    """

    def __init__(self, rotor, projector):
        """
        Constructor of the class.
        :param rotor: An instance of the class Rotor that performs rotation.
        :param projector: An instance of the class that performs projection.
        """
        self._rotor = rotor
        self._projector = projector

    @property
    def rotor(self):
        return self._rotor

    @property
    def projector(self):
        return self._projector

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return NotImplemented

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(tuple(sorted(self.__dict__.items())))

    def convert_point(self, lat, lon):
        """
        Calculates the Cartesian coordinates of the given point.
        :param lat: Latitude of the point (in degrees).
        :param lon: Longitude of the point (in degrees).
        :return: A tuple (x, y) of the Cartesian coordinates of the given point
        on the projection plane.
        """
        return self._projector.convert_point(
            *self._rotor.convert_point(lat, lon))

    def restore_point(self, x, y):
        """
        Calculates the spherical coordinates of the given point.
        :param x: Coordinate along the X-axis.
        :param y: Coordinate along the Y-axis.
        :return: A tuple (lat, lon) of the spherical coordinates of the given
        point.
        """
        return self._rotor.restore_point(*self._projector.restore_point(x, y))

    def convert_vector(self, u, v, lat, lon, return_point=False):
        """
        Calculates the components along the X- and Y- axes of the given vector
        that originates from the given point.
        :param u: Zonal component of the vector.
        :param v: Meridional component of the vector.
        :param lat: Latitude (in degrees) of the vector's origin.
        :param lon: Longitude (in degrees) of the vector's origin.
        :param return_point: Flag that tells the method to return the Cartesian
        coordinates of the vector's origin along with its components.
        :return: A tuple of the components of the vector along the X- and Y-
        axes respectively. If the flag return_point is True, than the result
        tuple is extended with the Cartesian (x, y) coordinates of the vector's
        origin.
        """
        rot_u, rot_v, rot_lat, rot_lon = self._rotor.convert_vector(u, v, lat,
                                                                    lon, True)
        return self._projector.convert_vector(rot_u, rot_v, rot_lat, rot_lon,
                                              return_point)

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
        rot_u, rot_v, rot_lat, rot_lon = self._projector.restore_vector(u, v, x,
                                                                        y, True)
        return self._rotor.restore_vector(rot_u, rot_v, rot_lat, rot_lon,
                                          return_point)


def restore_points(xx, yy, converter):
    """
    Converts Cartesian coordinates into spherical.
    :param xx: List of lists with the X-coordinates.
    :param yy: List of lists with the Y-coordinates.
    :param converter: Converter to use for the transformations.
    :return: Two 2D arrays with latitudes and longitudes of the given points.
    """
    res_lat, res_lon = \
        np.zeros((len(xx), len(xx[0]))), np.zeros((len(xx), len(xx[0])))
    for i in range(len(xx)):
        for j in range(len(xx[i])):
            res_lat[i, j], res_lon[i, j] = \
                converter.restore_point(xx[i][j], yy[i][j])
    return res_lat, res_lon


def convert_points(la, lo, converter, progress_callback=None):
    """
    Converts spherical coordinates into Cartesian.
    :param la: List of lists with latitudes.
    :param lo: List of lists with longitudes.
    :param converter: Converter to use for the transformations.
    :param progress_callback: Pointer to a function that is called after each
    row is processed.
    :return: Two 2D arrays with latitudes and longitudes of
    the given points.
    """
    res_x, res_y = \
        np.zeros((len(la), len(la[0]))), np.zeros((len(la), len(la[0])))
    row_count = len(la)
    for i in range(row_count):
        if progress_callback is not None:
            progress_callback(i, row_count)
        for j in range(len(la[i])):
            res_x[i, j], res_y[i, j] = \
                converter.convert_point(la[i][j], lo[i][j])
    progress_callback(row_count, row_count)
    return res_x, res_y
