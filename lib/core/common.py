import numpy as np

HALF_PI = np.pi / 2.0
POLE_TOLERANCE = 1e-5


def cos_sin_deg(angle_degs):
    """
    Calculates cosines and sinuses of angles.
    :param angle_degs: Scalar or array of angles (in degrees).
    :return: Tuple of cosines and sinuses of the angles.
    """
    # We want float64 precision here.
    angle_degs = np.asanyarray(angle_degs)
    angle_rads = np.empty(angle_degs.shape)
    np.radians(angle_degs, angle_rads)
    return np.cos(angle_rads), np.sin(angle_rads)


def gen_rot_matrices_deg(angle_degs, gen_reverse=False):
    """
    Generates rotation matrices.
    :param angle_degs: Scalar or array of rotation angles (in degrees).
    :param gen_reverse: Boolean flag that specifies whether an array of
    matrices for backward rotation is return along with an array of matrices
    for forward rotation.
    :return: Rotation matrices.
    """
    cosines, sines = cos_sin_deg(angle_degs)
    return gen_rot_matrices_sincos(sines, cosines, gen_reverse)


def gen_rot_matrices_rad(angle_rad, gen_reverse=False):
    return gen_rot_matrices_sincos(np.sin(angle_rad), np.cos(angle_rad),
                                   gen_reverse)


def gen_rot_matrices_sincos(sines, cosines, gen_reverse=False):
    result = np.asanyarray([[cosines, -sines], [sines, cosines]])
    if not gen_reverse:
        return result
    else:
        return result, np.swapaxes(result, 0, 1)


def apply_rot_matrices(uu, vv, rot_matrices):
    """
    Performs rotation of given vectors applying given array of rotation
    matrices.
    :param uu: Scalar or array of the first vectors' components.
    :param vv: Scalar or array of the second vectors' components.
    :param rot_matrices: Rotation matrices.
    :return: Tuple of scalars or arrays of components of the rotated vectors.
    """
    stacked = np.ma.concatenate([uu[np.newaxis, ...], vv[np.newaxis, ...]])
    rot_vecs = np.einsum('ij...,j...', rot_matrices, stacked)

    rot_uu = rot_vecs[..., 0]
    rot_vv = rot_vecs[..., 1]

    mask = np.ma.getmask(uu)
    if mask is not np.ma.nomask:
        rot_uu = np.ma.masked_where(mask, rot_uu)
        rot_vv = np.ma.masked_where(mask, rot_vv)

    return rot_uu, rot_vv


def rotate_vectors_deg(uu, vv, angle_degs):
    """
    Rotates 2D vectors by given angles.
    :param uu: Scalar or array of the first vectors' components.
    :param vv: Scalar or array of the second vectors' components.
    :param angle_degs: Scalar or array of rotation angles (in degrees).
    :return: Tuple of scalars or arrays of components of the rotated vectors.
    """
    return apply_rot_matrices(uu, vv, gen_rot_matrices_deg(angle_degs))
