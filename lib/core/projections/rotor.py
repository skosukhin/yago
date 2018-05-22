import itertools

import numpy as np

from core.common import cos_sin_deg, almost_equal


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

    _ROT_MATRICES = {
        'X':
            lambda c, s: np.array([[1, 0.0, 0.0],
                                   [0.0, c, -s],
                                   [0.0, s, c]]),
        'Y':
            lambda c, s: np.array([[c, 0.0, s],
                                   [0.0, 1.0, 0.0],
                                   [-s, 0.0, c]]),
        'Z':
            lambda c, s: np.array([[c, -s, 0.0],
                                   [s, c, 0.0],
                                   [0.0, 0.0, 1.0]])}

    def __init__(self, axes_names, angles_deg, simplify=False):
        """
        Constructor of the class.
        :param axes_names: iterable containing chars that indicate axes
        of rotations: 'X', 'Y', 'Z'.
        :param angles_deg: iterable containing angles
        of rotations (counter-clockwise, in degrees).
        :param simplify: flag that indicates whether or not the
        given sequence of rotations must be simplified
        """

        rot_axes_names = []
        rot_angles_deg = []
        for name, angle in itertools.izip(axes_names, angles_deg):
            if name not in Rotor._ROT_MATRICES:
                raise Exception('Unexpected axis name: ' + name +
                                '. Expected values are: ' +
                                ', '.join(Rotor._ROT_MATRICES.keys()))
            rot_axes_names.append(name)
            rot_angles_deg.append(angle)

        if simplify:
            rot_axes_names, rot_angles_deg =\
                Rotor._simplify_sequence(rot_axes_names, rot_angles_deg)

        if rot_axes_names:
            name, angle = rot_axes_names[0], rot_angles_deg[0]
            rot_matrix_to = Rotor._ROT_MATRICES[name](*cos_sin_deg(angle))
            for name, angle in itertools.izip(rot_axes_names[1:],
                                              rot_angles_deg[1:]):
                rot_matrix_to = np.dot(
                    Rotor._ROT_MATRICES[name](*cos_sin_deg(angle)),
                    rot_matrix_to)
        else:
            rot_matrix_to = np.eye(3)

        self._rot_axes_names = tuple(rot_axes_names)
        self._rot_angles_deg = tuple(rot_angles_deg)
        self._rot_matrix_to = rot_matrix_to
        self._rot_matrix_from = np.transpose(rot_matrix_to)
        self._rot_matrix_from.flags.writeable = False
        self._rot_matrix_to.flags.writeable = False

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self._rot_axes_names == other._rot_axes_names and
                    self._rot_angles_deg == other._rot_angles_deg and
                    np.array_equal(self._rot_matrix_to, other._rot_matrix_to))
        return NotImplemented

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((self._rot_axes_names, self._rot_angles_deg,
                     self._rot_matrix_to.data))

    @property
    def rot_axes_names(self):
        return self._rot_axes_names

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
    def build_rotor(geo_lat, geo_lon, rot_lat, rot_lon, norm_angle):
        """
        Creates an instance of the class Rotor that represents a sequence of
        rotations of the geographical coordinate system so that a point with
        given latitude and longitude gets given coordinates in the rotated
        lat/lon coordinate system.
        :param geo_lat: Latitude of the point in the geographical coordinate
        system (in degrees).
        :param geo_lon:  Longitude of the point in the geographical coordinate
        system (in degrees).
        :param rot_lat: Latitude of the point in the rotated lat/lon coordinate
        system (in degrees).
        :param rot_lon: Longitude of the point in the rotated lat/lon
        coordinate system (in degrees).
        :param norm_angle: Additional angle of rotation of the coordinate
        system around the normal at the point (counter-clockwise, in degrees).
        :return: Instance of the class Rotor.
        """
        return Rotor('ZYZYZ', [-geo_lon, geo_lat - 90.0,
                               norm_angle,
                               90.0 - rot_lat, rot_lon],
                     simplify=True)

    @staticmethod
    def _rotate_points(rot_matrix, lats, lons):
        c_lats, s_lats = cos_sin_deg(lats)
        c_lons, s_lons = cos_sin_deg(lons)
        orig_normals = Rotor._build_normals(c_lats, s_lats, c_lons, s_lons)
        rot_normals = np.einsum('ij,...j', rot_matrix, orig_normals)
        rot_lats = np.degrees(np.arcsin(rot_normals[..., 2]))
        rot_lons = np.degrees(np.arctan2(rot_normals[..., 1],
                                         rot_normals[..., 0]))
        return rot_lats, rot_lons

    @staticmethod
    def _rotate_vectors(rot_matrix, uu, vv, lats, lons, return_point=False):
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

        mask = np.ma.getmask(uu)
        if mask is not np.ma.nomask:
            rot_uu = np.ma.masked_where(mask, rot_uu)
            rot_vv = np.ma.masked_where(mask, rot_vv)

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
    def _normalize_angle_deg(angle_deg):
        return (angle_deg + 360) % 360

    @staticmethod
    def _simplify_sequence(axes_names, angles_deg):
        result_axes_names, result_angles_deg = [], []

        for name, angle in itertools.izip(axes_names, angles_deg):

            angle = Rotor._normalize_angle_deg(angle)
            if almost_equal(angle, 0):
                continue

            if result_axes_names and result_axes_names[-1] == name:
                new_angle = Rotor._normalize_angle_deg(
                    result_angles_deg[-1] + angle)

                if almost_equal(new_angle, 0):
                    result_axes_names.pop()
                    result_angles_deg.pop()
                else:
                    result_axes_names[-1] = name
                    result_angles_deg[-1] = new_angle
            else:
                result_axes_names.append(name)
                result_angles_deg.append(angle)

        return result_axes_names, result_angles_deg
