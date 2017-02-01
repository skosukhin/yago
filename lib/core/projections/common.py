import numpy as np

HALF_PI = np.pi / 2.0
QUARTER_PI = np.pi / 4.0


def cos_sin_deg(angle_deg):
    """
    Calculates cosines and sinuses of angles.
    :param angle_deg: Scalar or array of angles (in degrees).
    :return: Tuple of cosines and sinuses of the angles.
    """
    angle_rad = np.radians(angle_deg)
    return np.cos(angle_rad), np.sin(angle_rad)


def rotate_vectors(uu, vv, angles):
    """
    Rotates 2D vectors by given angles.
    :param uu: Scalar or array of the first vectors' components.
    :param vv: Scalar or array of the second vectors' components.
    :param angles: Scalar or array of rotation angles (in degrees).
    :return: Tuple of scalars or arrays of components of the rotated vectors.
    """
    c_lons, s_lons = cos_sin_deg(angles)
    rot_matrices = np.asanyarray([[c_lons, -s_lons], [s_lons, c_lons]])
    rot_vecs = np.einsum('ij...,j...', rot_matrices, np.stack([uu, vv]))
    rot_uu = rot_vecs[..., 0]
    rot_vv = rot_vecs[..., 1]
    return rot_uu, rot_vv
