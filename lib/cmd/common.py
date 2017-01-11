import os
import sys
import time

import numpy as np

import cmd.name_constants as names
from core.converter import Converter
from core.projections import projections
from core.rotors import Rotor, RotorZ, RotorY, RotorX


def parse_list_of_floats(string):
    return np.array([np.float64(val) for val in parse_list_of_strings(string)])


def parse_list_of_strings(string):
    return string.split(';')


def parse_pos_intp(string):
    result = np.intp(string)
    if result <= 0:
        raise Exception('Argument must be positive.')
    return result


def parse_pos_float(string):
    result = np.float64(string)
    if result <= 0:
        raise Exception('Argument must be positive.')
    return result


def generate_cartesian_grid(x_start, x_count, x_step,
                            y_start, y_count, y_step):
    x_array = np.linspace(x_start, x_start + (x_count - 1) * x_step,
                          num=x_count)
    y_array = np.linspace(y_start, y_start + (y_count - 1) * y_step,
                          num=y_count)

    return np.meshgrid(x_array, y_array)


def copy_nc_attributes(src_var, dst_var):
    for attr_name in src_var.ncattrs():
        dst_var.setncattr(attr_name, src_var.getncattr(attr_name))


def copy_dim_var(src_ds, dst_ds, dim_var_name, return_data=False):
    src_dim = src_ds.dimensions[dim_var_name]
    src_var = src_ds.variables[dim_var_name]
    dst_ds.createDimension(dim_var_name, src_dim.size)
    dst_var = dst_ds.createVariable(dim_var_name, src_var.dtype,
                                    dimensions=(dim_var_name,))
    src_var_list = src_var[:]
    dst_var[:] = src_var_list
    copy_nc_attributes(src_var, dst_var)

    if return_data:
        return src_var_list


def set_generic_lat_attributes(lat_var):
    lat_var.units = 'degrees_north'
    lat_var.long_name = 'latitude coordinate'
    lat_var.standard_name = 'latitude'


def set_generic_lon_attributes(lon_var):
    lon_var.units = 'degrees_east'
    lon_var.long_name = 'longitude coordinate'
    lon_var.standard_name = 'longitude'


def add_or_append_history(dataset, ignored_args=None):
    add_history(dataset, get_history(dataset), ignored_args)


def add_history(dataset, history_suffix=None, ignored_args=None):
    history = gen_hist_string(ignored_args)
    if history_suffix:
        history = [history + '\n', history_suffix]
    dataset.setncattr(names.ATTR_HISTORY, history)


def get_history(dataset):
    if hasattr(dataset, names.ATTR_HISTORY):
        return dataset.getncattr(names.ATTR_HISTORY)


def gen_hist_string(ignored_args=None):
    if ignored_args:
        tup = tuple(ignored_args)
        args = [arg for arg in sys.argv[1:] if not arg.startswith(tup)]
    else:
        args = sys.argv[1:]

    return (time.ctime(time.time()) + ': ' + os.path.basename(
        sys.argv[0]) + ' ' + ' '.join(args))


def check_preprocessed(dataset):
    if (not hasattr(dataset, names.ATTR_HISTORY) or
            'arctic preproc' not in dataset.getncattr(names.ATTR_HISTORY)):
        raise Exception()


def init_converter_from_proj_var(proj_var):
    r = _decode_rotor(proj_var.rot_axes, proj_var.rot_angles_deg)
    p = projections[proj_var.short_name].init(proj_var.earth_radius,
                                              proj_var.standard_parallel)

    return Converter(r, p)


def init_converter_from_args(args):
    p = projections[args.proj_name].init(args.earth_radius,
                                         args.true_scale_lats)
    r = p.build_rotor(args.orig_lat, args.orig_lon, args.adjust_angle)

    return Converter(r, p)


def create_dir_for_file(filename):
    try:
        os.makedirs(os.path.dirname(filename))
    except:
        pass


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
    return Rotor.chain(*rotors)
