import numpy as np

from core.common import cos_sin_deg


class Rotor(object):
    """
    Class that links regular spherical (lat/lon) coordinates and the
    corresponding vectors with the lat/lon coordinates in the rotated lat/lon
    system. Rotations are performed around axes of a 3D right handed Cartesian
    coordinate system. The origin of the Cartesian system is in the center of
    the sphere that approximates the Earth. X-axis points to the intersection
    of the equator and the Greenwich Meridian; Y-axis points to the
    intersection of the equator and the 90th eastern meridian; Z-axis points to
    the Northern Pole.
    """

    def __init__(self):
        self._rot_axes_ids = ()
        self._rot_angles_deg = ()
        self._rot_matrix_to = np.array(
            [[1.0, 0.0, 0.0],
             [0.0, 1.0, 0.0],
             [0.0, 0.0, 1.0]])
        self._rot_matrix_from = np.transpose(self._rot_matrix_to)
        self._rot_matrix_from.flags.writeable = False
        self._rot_matrix_to.flags.writeable = False

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self._rot_axes_ids == other._rot_axes_ids and
                    self._rot_angles_deg == other._rot_angles_deg and
                    np.array_equal(self._rot_matrix_to, other._rot_matrix_to))
        return NotImplemented

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((self._rot_axes_ids, self._rot_angles_deg,
                     self._rot_matrix_to.data))

    @property
    def rot_axes_ids(self):
        return self._rot_axes_ids

    @property
    def rot_angles_deg(self):
        return self._rot_angles_deg

    def convert_point(self, lat, lon):
        """
        Calculates the rotated lat/lon coordinates of the given point.
        :param lat: Latitude of the point (in degrees) in the regular lat/lon
        system.
        :param lon: Longitude of the point (in degrees) in the regular lat/lon
        system.
        :return: A tuple (rot_lat, rot_lon) of the rotated lat/lon coordinates
        of the given point.
        """
        return Rotor._rotate_point(self._rot_matrix_to, lat, lon)

    def restore_point(self, rot_lat, rot_lon):
        """
        Calculates the regular lat/lon coordinates of the given point.
        :param rot_lat: Latitude of the point (in degrees) in the rotated
        lat/lon system.
        :param rot_lon: Longitude of the point (in degrees) in the rotated
        lat/lon system.
        :return: A tuple (lat, lon) of the regular lat/lon coordinates of the
        given point.
        """
        return Rotor._rotate_point(self._rot_matrix_from, rot_lat, rot_lon)

    def convert_vector(self, u, v, lat, lon, return_point=False):
        """
        Calculates the rotated zonal and meridional components of the given
        vector that originates from the given point.
        :param u: True zonal component of the vector.
        :param v: True meridional component of the vector.
        :param lat: True latitude (in degrees) of the vector's origin.
        :param lon: True longitude (in degrees) of the vector's origin.
        :param return_point: Flag that tells the method to return the rotated
        lat/lon coordinates of the vector's origin along with its components.
        :return: A tuple of the rotated zonal and meridional components of
        the vector respectively. If the flag return_point is True, than the
        result tuple is extended with the rotated (lat, lon) coordinates of the
        vector's origin.
        """
        return Rotor._rotate_vector(
            self._rot_matrix_to, lat, lon, u, v, return_point
        )

    def restore_vector(self, rot_u, rot_v, rot_lat, rot_lon,
                       return_point=False):
        """
        Calculates the zonal and meridional components of the given vector that
        originates from the given point.
        :param rot_u: False zonal component of the vector.
        :param rot_v: False meridional component of the vector.
        :param rot_lat: False latitude (in degrees) of the vector's origin.
        :param rot_lon: False longitude (in degrees) of the vector's origin.
        :param return_point: Flag that tells the method to return the spherical
        coordinates of the vector's origin along with its components.
        :return: A tuple of the true zonal and true meridional components of
        the vector respectively. If the flag return_point is True, than the
        result tuple is extended with the regular spherical (lat, lon)
        coordinates of the vector's origin.
        """
        return Rotor._rotate_vector(self._rot_matrix_from, rot_lat, rot_lon,
                                    rot_u,
                                    rot_v, return_point)

    @staticmethod
    def chain(*rotors):
        result = Rotor()
        for rotor in rotors:
            result._rot_axes_ids = result._rot_axes_ids + rotor._rot_axes_ids
            result._rot_angles_deg = \
                result._rot_angles_deg + rotor._rot_angles_deg
            result._rot_matrix_to = np.dot(
                rotor._rot_matrix_to,
                result._rot_matrix_to)
        result._rot_matrix_from = np.transpose(result._rot_matrix_to)
        result._rot_matrix_from.flags.writeable = False
        result._rot_matrix_to.flags.writeable = False
        return result

    @staticmethod
    def _rotate_point(rot_matrix, lat, lon):
        lat, lon = Rotor._resolve_polar_point(lat, lon)
        orig_normal = Rotor._build_normal(lat, lon)
        rot_normal = np.dot(rot_matrix, orig_normal).ravel()

        return Rotor._resolve_polar_point(
            np.degrees(np.arcsin(rot_normal[2])),
            np.degrees(np.arctan2(rot_normal[1], rot_normal[0]))
        )

    @staticmethod
    def _rotate_vector(rot_matrix, lat, lon, u, v, return_point=False):
        lat, lon = Rotor._resolve_polar_point(lat, lon)
        east, north = Rotor._build_east(lon), Rotor._build_north(lat, lon)

        vector_3d = Rotor._from_uv_to_3d(u, v, east, north)
        rot_vector_3d = np.dot(rot_matrix, vector_3d)
        rot_lat, rot_lon = Rotor._rotate_point(rot_matrix, lat, lon)

        rot_east = Rotor._build_east(rot_lon)
        rot_north = Rotor._build_north(rot_lat, rot_lon)

        rot_u, rot_v = Rotor._from_3d_to_uv(rot_vector_3d, rot_east, rot_north)

        if return_point:
            return rot_u, rot_v, rot_lat, rot_lon
        else:
            return rot_u, rot_v

    @staticmethod
    def _build_normal(lat, lon):
        c_lat, s_lat = cos_sin_deg(lat)
        c_lon, s_lon = cos_sin_deg(lon)
        return np.array([c_lat * c_lon, c_lat * s_lon, s_lat])

    @staticmethod
    def _build_east(lon):
        c_lon, s_lon = cos_sin_deg(lon)
        return np.array([-s_lon, c_lon, 0.0])

    @staticmethod
    def _build_north(lat, lon):
        c_lat, s_lat = cos_sin_deg(lat)
        c_lon, s_lon = cos_sin_deg(lon)
        return np.array([-s_lat * c_lon, -s_lat * s_lon, c_lat])

    @staticmethod
    def _from_uv_to_3d(u, v, east, north):
        return np.add(east * u, north * v)

    @staticmethod
    def _from_3d_to_uv(vector, east, north):
        return np.dot(vector, east), np.dot(vector, north)

    @staticmethod
    def _resolve_polar_point(lat, lon):
        eps = np.finfo(float).eps
        if np.fabs(np.fabs(lat) - 90.0) <= eps:
            lon = 0.0
        return lat, lon


class RotorX(Rotor):
    def __init__(self, angle_deg):
        self._rot_axes_ids = ('X',)
        self._rot_angles_deg = (angle_deg,)
        c, s = cos_sin_deg(angle_deg)
        self._rot_matrix_to = np.array(
            [[1, 0.0, 0.0],
             [0.0, c, -s],
             [0.0, s, c]])
        self._rot_matrix_from = np.transpose(self._rot_matrix_to)
        self._rot_matrix_from.flags.writeable = False
        self._rot_matrix_to.flags.writeable = False


class RotorY(Rotor):
    def __init__(self, angle_deg):
        self._rot_axes_ids = ('Y',)
        self._rot_angles_deg = (angle_deg,)
        c, s = cos_sin_deg(angle_deg)
        self._rot_matrix_to = np.array(
            [[c, 0.0, s],
             [0.0, 1.0, 0.0],
             [-s, 0.0, c]])
        self._rot_matrix_from = np.transpose(self._rot_matrix_to)
        self._rot_matrix_from.flags.writeable = False
        self._rot_matrix_to.flags.writeable = False


class RotorZ(Rotor):
    def __init__(self, angle_deg):
        self._rot_axes_ids = ('Z',)
        self._rot_angles_deg = (angle_deg,)
        c, s = cos_sin_deg(angle_deg)
        self._rot_matrix_to = np.array(
            [[c, -s, 0.0],
             [s, c, 0.0],
             [0.0, 0.0, 1.0]])
        self._rot_matrix_from = np.transpose(self._rot_matrix_to)
        self._rot_matrix_from.flags.writeable = False
        self._rot_matrix_to.flags.writeable = False
