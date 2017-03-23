import numpy as np

HALF_PI = np.pi / 2.0
QUARTER_PI = np.pi / 4.0
NORTH_POLE_TOLERANCE = 1e-5


def cos_sin_deg(angle_deg):
    """
    Calculates cosines and sinuses of angles.
    :param angle_deg: Scalar or array of angles (in degrees).
    :return: Tuple of cosines and sinuses of the angles.
    """
    angle_rad = np.radians(angle_deg)
    return np.cos(angle_rad), np.sin(angle_rad)


def gen_rot_matrices(angle_degs):
    """
    Generates rotation matrices.
    :param angle_degs: Scalar or array of rotation angles (in degrees).
    :return: Rotation matrices.
    """
    c_lons, s_lons = cos_sin_deg(angle_degs)
    return np.asanyarray([[c_lons, -s_lons], [s_lons, c_lons]])


def apply_rot_matrices(uu, vv, rot_matrices):
    """
    Performs rotation of given vectors applying given array of rotation
    matrices.
    :param uu: Scalar or array of the first vectors' components.
    :param vv: Scalar or array of the second vectors' components.
    :param rot_matrices: Rotation matrices.
    :return: Tuple of scalars or arrays of components of the rotated vectors.
    """
    rot_vecs = np.einsum('ij...,j...', rot_matrices, np.stack([uu, vv]))
    rot_uu = rot_vecs[..., 0]
    rot_vv = rot_vecs[..., 1]
    return rot_uu, rot_vv


def rotate_vectors(uu, vv, angles):
    """
    Rotates 2D vectors by given angles.
    :param uu: Scalar or array of the first vectors' components.
    :param vv: Scalar or array of the second vectors' components.
    :param angles: Scalar or array of rotation angles (in degrees).
    :return: Tuple of scalars or arrays of components of the rotated vectors.
    """
    return apply_rot_matrices(uu, vv, gen_rot_matrices(angles))
