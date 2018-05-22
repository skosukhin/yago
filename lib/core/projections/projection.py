class Projection(object):
    short_name = None
    long_name = None
    standard_name = None

    true_scale_lats = None
    earth_radius = None

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return NotImplemented

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(tuple(sorted(self.__dict__.items())))

    @classmethod
    def unified_init(cls, earth_radius, true_lats):
        """
        Creates an instance of the class given a unified set of arguments.
        :param earth_radius: Earth radius (in meters).
        :param true_lats: List of latitudes of true scale (in degrees).
        :return: Instance of the class Projection.
        """
        raise NotImplementedError()

    @property
    def reference_point(self):
        """
        Returns geographical coordinates of the reference point of the
        projection.
        :return: Tuple of latitude and longitude (in degrees) of the reference
        point.
        """
        raise NotImplementedError()

    def convert_points(self, lats, lons):
        """
        Calculates Cartesian coordinates of points.
        :param lats: Scalar or array of geographical latitudes of points
        (in degrees).
        :param lons: Scalar or array of geographical longitudes of points
        (in degrees).
        :return: Tuple of scalars or arrays of cartesian coordinates of the
        points (in meters).
        """
        raise NotImplementedError()

    def restore_points(self, xx, yy):
        """
        Calculates geographical coordinates of points.
        :param xx: Scalar or array of x-coordinates of points (in meters).
        :param yy: Scalar or array of y-coordinates of points (in meters).
        :return: Tuple of scalars or arrays of geographical coordinates of the
        points (in degrees).
        """
        raise NotImplementedError()

    def convert_vectors(self, uu, vv, lats, lons, return_points=False):
        """
        Calculates X- and Y-components of vectors.
        :param uu: Scalar or array of zonal components of vectors.
        :param vv: Scalar or array of meridional components of vectors.
        :param lats: Scalar or array of geographical latitudes of vectors'
        origins (in degrees).
        :param lons: Scalar or array of geographical longitudes of vectors'
        origins (in degrees).
        :param return_points: Boolean flag that tells the method to include
        Cartesian coordinates of the vectors' origins into the output tuple.
        :return: Tuple of scalars or arrays of x- and y-components of the
        vectors. If the flag return_points is set to True, than the result
        tuple is extended with scalars or arrays of Cartesian coordinates of
        the vectors' origins (in meters).
        """
        raise NotImplementedError()

    def restore_vectors(self, uu, vv, xx, yy, return_points=False):
        """
        Calculates zonal and meridional components of vectors.
        :param uu: Scalar or array of x-components of vectors.
        :param vv: Scalar or array of y-components of vectors.
        :param xx: Scalar or array of x-coordinates of vectors' origins
        (in meters).
        :param yy: Scalar or array of Y-coordinates of vectors' origins
        (in meters).
        :param return_points: Boolean flag that tells the method to include
        geographical coordinates of the vectors' origins into the output tuple.
        :return: Tuple of scalars or arrays of zonal and meridional components
        of the vectors. If the flag return_points is set to True, than the
        result tuple is extended with scalars or arrays of geographical
        coordinates of the vectors' origins (in degrees).
        """
        raise NotImplementedError()

    def get_scale_factors(self, lats, lons):
        """
        Calculates scale distortions induces by the projection in given points.
        :param lats: Scalar or array of geographical latitudes of points
        (in degrees).
        :param lons: Scalar or array of geographical longitudes of points
        (in degrees).
        :return: Tuple of scalars or arrays of zonal and meridional distortions
        respectively. In the case of a conformal projection, both elements of
        the tuple reference the same array.
        """
        raise NotImplementedError()
