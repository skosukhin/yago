import numpy as np


class Converter(object):
    """
    Class that combines the rotation and projection transformations.
    """

    def __init__(self, rotor, projection):
        """
        Constructor of the class.
        :param rotor: An instance of the class Rotor
        :param projection: An instance of the class Projection.
        """
        self._rotor = rotor
        self._projection = projection

    @property
    def rotor(self):
        return self._rotor

    @property
    def projection(self):
        return self._projection

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return NotImplemented

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(tuple(sorted(self.__dict__.items())))

    def convert_points(self, lats, lons):
        """
        Calculates Cartesian coordinates of geographical points.
        :param lats: 2D array of geographical latitudes of the points
        (in degrees).
        :param lons: 2D array of geographical longitudes of the points
        (in degrees).
        :return: Tuple of 2D arrays (xx, yy) of Cartesian coordinates of the
        points.
        """
        res_x, res_y = self._rotor.convert_points(lats, lons)
        for i in xrange(lats.shape[0]):
            for j in xrange(lats.shape[1]):
                res_x[i, j], res_y[i, j] = \
                    self._projection.convert_point(res_x[i, j], res_y[i, j])
        return res_x, res_y

    def restore_points(self, xx, yy):
        """
        Calculates geographical coordinates of Cartesian point.
        :param xx: 2D array of x-coordinates of the points (in meters).
        :param yy: 2D array of y-coordinates of the points (in meters).
        :return: Tuple of 2D arrays (lat, lon) of geographical coordinates of
        the points.
        """
        res_lats, res_lons = \
            np.zeros(xx.shape), np.zeros(xx.shape)
        for i in xrange(xx.shape[0]):
            for j in xrange(xx.shape[1]):
                res_lats[i, j], res_lons[i, j] = \
                    self._projection.restore_point(xx[i, j], yy[i, j])
        res_lats, res_lons = self._rotor.restore_points(res_lats, res_lons)
        return res_lats, res_lons

    def convert_vectors(self, uu, vv, lats, lons, return_point=False):
        """
        Calculates X- and Y-components of given vectors.
        :param uu: 2D array of zonal components of the vectors.
        :param vv: 2D array of meridional components of the vectors.
        :param lats: 2D array of geographical latitudes of the vectors' origins
        (in degrees).
        :param lons: 2D array of geographical longitudes of vectors' origins
        (in degrees).
        :param return_point: Flag that tells the method to include Cartesian
        coordinates of the vectors' origins into the output tuple.
        :return: Tuple of X- and Y-components of the vectors. If the flag
        return_point is True, than the result tuple is extended with Cartesian
        (xx, yy) coordinates of the vectors' origins.
        """

        rot_uu, rot_vv, rot_lats, rot_lons = \
            self._rotor.convert_vectors(uu, vv, lats, lons, True)

        if return_point:
            for i in xrange(uu.shape[0]):
                for j in xrange(uu.shape[1]):
                    (rot_uu[i, j], rot_vv[i, j],
                     rot_lats[i, j], rot_lons[i, j]) = \
                        self._projection.convert_vector(rot_uu[i, j],
                                                        rot_vv[i, j],
                                                        rot_lats[i, j],
                                                        rot_lons[i, j],
                                                        True)
            return rot_uu, rot_vv, rot_lats, rot_lons
        else:
            for i in xrange(uu.shape[0]):
                for j in xrange(uu.shape[1]):
                    rot_uu[i, j], rot_vv[i, j] = \
                        self._projection.convert_vector(rot_uu[i, j],
                                                        rot_vv[i, j],
                                                        rot_lats[i, j],
                                                        rot_lons[i, j],
                                                        False)
            return rot_uu, rot_vv

    def restore_vectors(self, uu, vv, xx, yy, return_point=False):
        """
        Calculates zonal and meridional components of given vectors.
        :param uu: 2D array of x-components of the vectors.
        :param vv: 2D array of y-components of the vectors.
        :param xx: 2D array of x-coordinates of the vector's origins
        (in meters).
        :param yy: 2D array of y-coordinates of the vector's origins.
        :param return_point: Flag that tells the method to include geographical
        coordinates of the vectors' origins into the output tuple.
        :return: Tuple of zonal and meridional components of the vectors. If
        the flag return_point is True, than the result tuple is extended with
        the geographical (lats, lons) coordinates of the vectors' origins.
        """
        rot_uu, rot_vv = np.zeros(uu.shape), np.zeros(vv.shape)
        rot_lats, rot_lons = np.zeros(uu.shape), np.zeros(vv.shape)

        for i in xrange(uu.shape[0]):
            for j in xrange(uu.shape[1]):
                (rot_uu[i, j], rot_vv[i, j],
                 rot_lats[i, j], rot_lons[i, j]) = \
                    self._projection.restore_vector(uu[i, j],
                                                    vv[i, j],
                                                    xx[i, j],
                                                    yy[i, j],
                                                    True)

        return self._rotor.restore_vectors(rot_uu, rot_vv, rot_lats, rot_lons,
                                           return_point)
