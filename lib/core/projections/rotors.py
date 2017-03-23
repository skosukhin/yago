import numpy as np

from core.common import cos_sin_deg, NORTH_POLE_TOLERANCE


class Rotor(object):
    """
    Links the geographical coordinate system with a rotated lat/lon coordinate
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

    def convert_points(self, lats, lons):
        """
        Calculates rotated lat/lon coordinates of points.
        :param lats: Scalar or array of geographical latitudes of points
        (in degrees).
        :param lons: Scalar or array of geographical longitudes of points
        (in degrees).
        :return: Tuple of scalars or arrays of rotated lat/lon coordinates of
        the points (in degrees).
        """
        lats, lons = np.asanyarray(lats), np.asanyarray(lons)
        return Rotor._rotate_points(self._rot_matrix_to, lats, lons)

    def restore_points(self, rlats, rlons):
        """
        Calculates geographical coordinates of points.
        :param rlats: Scalar or array of rotated latitudes of points
        (in degrees).
        :param rlons: Scalar or array of rotated longitudes of points
        (in degrees).
        :return: Tuple of scalar or arrays of geographical coordinates of the
        points (in degrees).
        """
        rlats, rlons = np.asanyarray(rlats), np.asanyarray(rlons)
        return Rotor._rotate_points(self._rot_matrix_from, rlats, rlons)

    def convert_vectors(self, uu, vv, lats, lons, return_points=False):
        """
        Calculates rotated zonal and meridional components of vectors.
        :param uu: Scalar or array of zonal components of vectors.
        :param vv: Scalar or array of meridional components of vectors.
        :param lats: Scalar or array of geographical latitudes of vectors'
        origins (in degrees).
        :param lons: Scalar or array of geographical longitudes of vectors'
        origins (in degrees).
        :param return_points: Boolean flag that tells the method to include
        rotated lat/lon coordinates of the vectors' origins into the output
        tuple.
        :return: Tuple of scalar or arrays of rotated zonal and meridional
        components of the vectors. If the flag return_points is set to True,
        than the result tuple is extended with scalars or arrays of rotated
        lat/lon coordinates of the vectors' origins (in degrees).
        """
        uu, vv = np.asanyarray(uu), np.asanyarray(vv)
        lats, lons = np.asanyarray(lats), np.asanyarray(lons)
        return Rotor._rotate_vectors(self._rot_matrix_to,
                                     uu, vv, lats, lons,
                                     return_points)

    def restore_vectors(self, rot_uu, rot_vv, rot_lats, rot_lons,
                        return_points=False):
        """
        Calculates zonal and meridional components of vectors.
        :param rot_uu: Scalar or array of rotated zonal component of vectors.
        :param rot_vv: Scalar or array of rotated meridional component of
        vectors.
        :param rot_lats: Scalar or array of rotated latitudes of vectors'
        origins (in degrees).
        :param rot_lons: Scalar or array of rotated longitudes of vectors'
        origins (in degrees).
        :param return_points: Boolean flag that tells the method to include
        geographical coordinates of the vectors' origins into the output tuple.
        :return: Tuple of scalars or arrays of zonal and meridional components
        of the vectors. If the flag return_points is set to True, than the
        result tuple is extended with geographical coordinates of the vectors'
        origins (in degrees).
        """
        rot_uu, rot_vv = np.asanyarray(rot_uu), np.asanyarray(rot_vv)
        rot_lats, rot_lons = np.asanyarray(rot_lats), np.asanyarray(rot_lons)
        return Rotor._rotate_vectors(self._rot_matrix_from,
                                     rot_uu, rot_vv, rot_lats, rot_lons,
                                     return_points)

    @staticmethod
    def chain(*rotors):
        """
        Combines a sequence of rotors into a single one.
        :param rotors: Rotors to be combined in a single rotation object.
        :return: Instance of the class Rotor that represents a sequence of the
        given rotations.
        """
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
    def _rotate_points(rot_matrix, lats, lons):
        np_tol = lats.dtype.type(NORTH_POLE_TOLERANCE)
        lats, lons = Rotor._resolve_polar_points(lats, lons, np_tol)
        c_lats, s_lats = cos_sin_deg(lats)
        c_lons, s_lons = cos_sin_deg(lons)
        orig_normals = Rotor._build_normals(c_lats, s_lats, c_lons, s_lons)
        rot_normals = np.einsum('ij,...j', rot_matrix, orig_normals)
        rot_lats = np.degrees(np.arcsin(rot_normals[..., 2]))
        rot_lons = np.degrees(np.arctan2(rot_normals[..., 1],
                                         rot_normals[..., 0]))
        return Rotor._resolve_polar_points(rot_lats, rot_lons, np_tol)

    @staticmethod
    def _rotate_vectors(rot_matrix, uu, vv, lats, lons, return_point=False):
        np_tol = lats.dtype.type(NORTH_POLE_TOLERANCE)
        lats, lons = Rotor._resolve_polar_points(lats, lons, np_tol)
        c_lats, s_lats = cos_sin_deg(lats)
        c_lons, s_lons = cos_sin_deg(lons)
        easts = Rotor._build_easts(c_lons, s_lons)
        norths = Rotor._build_norths(c_lats, s_lats, c_lons, s_lons)
        vectors_3d = \
            np.einsum('...i,...', easts, uu) + \
            np.einsum('...i,...', norths, vv)

        rot_vectors_3d = np.einsum('ij,...j', rot_matrix, vectors_3d)
        rot_lats, rot_lons = Rotor._rotate_points(rot_matrix, lats, lons)

        c_rot_lats, s_rot_lats = cos_sin_deg(rot_lats)
        c_rot_lons, s_rot_lons = cos_sin_deg(rot_lons)
        rot_easts = Rotor._build_easts(c_rot_lons, s_rot_lons)
        rot_norths = Rotor._build_norths(c_rot_lats, s_rot_lats,
                                         c_rot_lons, s_rot_lons)

        rot_uu = np.einsum('...i,...i', rot_vectors_3d, rot_easts)
        rot_vv = np.einsum('...i,...i', rot_vectors_3d, rot_norths)

        if return_point:
            return rot_uu, rot_vv, rot_lats, rot_lons
        else:
            return rot_uu, rot_vv

    @staticmethod
    def _build_normals(c_lats, s_lats, c_lons, s_lons):
        return np.stack([c_lats * c_lons,
                         c_lats * s_lons,
                         s_lats],
                        axis=-1)

    @staticmethod
    def _build_easts(c_lons, s_lons):
        return np.stack([-s_lons,
                         c_lons,
                         np.zeros(c_lons.shape, dtype=c_lons.dtype)],
                        axis=-1)

    @staticmethod
    def _build_norths(c_lats, s_lats, c_lons, s_lons):
        return np.stack([-s_lats * c_lons,
                         -s_lats * s_lons,
                         c_lats],
                        axis=-1)

    @staticmethod
    def _resolve_polar_points(lats, lons, eps):
        t90 = lats.dtype.type(90)
        mask = np.fabs(np.fabs(lats) - t90) <= eps
        lons = np.where(mask,
                        np.zeros(lons.shape, lons.dtype), lons)
        lats = np.where(mask,
                        np.ones(lats.shape, lats.dtype)
                        * t90
                        * np.sign(lats),
                        lats)
        return lats, lons


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
