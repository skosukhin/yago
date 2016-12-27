import numpy as np


def cos_sin_rad(angle_rad):
    return np.array([np.cos(angle_rad), np.sin(angle_rad)])


def cos_sin_deg(angle_deg):
    return cos_sin_rad(np.radians(angle_deg))


def build_2d_rotation_z_rad(angle_rad):
    c, s = cos_sin_rad(angle_rad)
    return np.array([[c, -s], [s, c]])


def adjust_angle_180(angle):
    while angle >= 180.0:
        angle -= 360.0
    while angle < -180.0:
        angle += 360.0
    return angle


def min_max(arr):
    mi = arr[0]
    ma = mi
    for val in arr[1:]:
        if val > ma:
            ma = val
        if val < mi:
            mi = val
    return mi, ma
