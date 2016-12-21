import os
import sys
import time

import numpy as np

from core.converter import Converter
from core.projections.lambert import LambertConformalProjector
from core.projections.mercator import MercatorProjector
from core.projections.polar_stereographic import PolarStereographicProjector
from core.rotors import Rotor, RotorZ, RotorY, RotorX


def parse_list_of_floats(string):
    return np.array([np.float64(val) for val in parse_list_of_strings(string)])


def parse_list_of_strings(string):
    return string.split(',')


def chain_rotors(*rotors):
    result = Rotor()
    for rotor in rotors:
        result.rot_axes_ids.extend(rotor.rot_axes_ids)
        result.rot_angles_deg.extend(rotor.rot_angles_deg)
        result.rot_matrix_to = np.dot(
            rotor.rot_matrix_to,
            result.rot_matrix_to)
    result.rot_matrix_from = np.transpose(result.rot_matrix_to)
    return result


def build_rotor_for_polar_stereographic(orig_lat, orig_lon,
                                        add_angle_deg):
    """
    The function generates an instance of the class Rotor to be used in
    conjunction with the Polar stereographic projection. The obtained instance
    is the result of three consecutive rotations: around Z-axis, around Y-axis,
    and again around Z-axis. The first two rotations shift the given point to
    the North Pole by rotating the coordinate system. The second rotation
    around Z-axis (the last one among the three) is optional to help users to
    adjust the orientation of the coordinate grid to account either for the
    features of the following projection procedure or for the plotting needs.
    :param orig_lat: Real latitude (in degrees) of the point that is supposed
    to be in the center of the projection.
    :param orig_lon: Real longitude (in degrees) of the point that is
    supposed to be in the center of the projection.
    :param add_angle_deg: Angle (in degrees) of the last rotation around the
    Z-axis.
    :return: Returns an instance of the class Rotor that rotates the regular
    lat/lon coordinate system to move a point that a user wants to be in the
    middle of the Polar stereographic projection to the North pole.
    """

    return chain_rotors(
        RotorZ(180.0 - orig_lon),
        RotorY(90.0 - orig_lat),
        RotorZ(add_angle_deg))


def build_rotor_for_mercator(orig_lat, orig_lon, add_angle_deg):
    rotor_adjust_center = build_rotor_for_polar_stereographic(orig_lat,
                                                              orig_lon,
                                                              add_angle_deg)
    return chain_rotors(rotor_adjust_center, RotorY(90.0))


def build_rotor_for_lambert(orig_lat, orig_lon, add_angle_deg):
    rotor_adjust_center = build_rotor_for_polar_stereographic(orig_lat,
                                                              orig_lon,
                                                              add_angle_deg)

    return chain_rotors(rotor_adjust_center, RotorY(45.0))


def generate_cartesian_grid(x_count, x_step, y_count, y_step):
    x_start = _calc_cartesian_start(x_count, x_step)
    y_start = _calc_cartesian_start(y_count, y_step)

    x_array = np.linspace(x_start, x_start + (x_count - 1) * x_step,
                          num=x_count)
    y_array = np.linspace(y_start, y_start + (y_count - 1) * y_step,
                          num=y_count)

    return np.meshgrid(x_array, y_array)


def copy_nc_attributes(src_var, dst_var):
    for attr_name in src_var.ncattrs():
        dst_var.setncattr(attr_name, src_var.getncattr(attr_name))


def set_generic_lat_attributes(lat_var):
    lat_var.units = 'degrees_north'
    lat_var.long_name = 'latitude coordinate'
    lat_var.standard_name = 'latitude'


def set_generic_lon_attributes(lon_var):
    lon_var.units = 'degrees_east'
    lon_var.long_name = 'longitude coordinate'
    lon_var.standard_name = 'longitude'


def gen_hist_string(ignored=None):
    if ignored:
        tup = tuple(ignored)
        args = [arg for arg in sys.argv[1:] if not arg.startswith(tup)]
    else:
        args = sys.argv[1:]

    return (time.ctime(time.time()) + ': ' + os.path.basename(
        sys.argv[0]) + ' ' + ' '.join(args))


def init_converter_from_proj_var(proj_var):
    r = _decode_rotor(proj_var.rot_axes, proj_var.rot_angles_deg)

    if proj_var.short_name == 'stereo':
        p = PolarStereographicProjector(proj_var.standard_parallel,
                                        proj_var.earth_radius)
    elif proj_var.short_name == 'mercator':
        p = MercatorProjector(proj_var.standard_parallel,
                              proj_var.earth_radius)
    elif proj_var.short_name == 'lambert':
        p = LambertConformalProjector(proj_var.standard_parallel[0],
                                      proj_var.standard_parallel[1],
                                      proj_var.earth_radius)
    else:
        raise Exception('Unknown projection.')

    return Converter(r, p)


def _calc_cartesian_start(count, step):
    # If the number of grid points along the axis is odd than we put the
    # central point to the origin of the coordinate system.
    if count % 2 == 1:
        start = -(count - 1) / 2.0 * step
    # If the number of grid points is even than we put the origin of the
    # coordinate system between two central points of the grid.
    else:
        start = (-count / 2.0 + 0.5) * step

    return start


def _decode_rotor(rot_axes_ids, rot_angles_deg):
    rotors = []
    for i, c in enumerate(rot_axes_ids):
        angle = rot_angles_deg[i]
        if c == 'X':
            rotors.append(RotorX(angle))
        elif c == 'Y':
            rotors.append(RotorY(angle))
        elif c == 'Z':
            rotors.append(RotorZ(angle))
        else:
            raise Exception('Unknown rotation axis ID: \'' + c + '\'.')
    return chain_rotors(*rotors)
