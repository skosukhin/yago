from datetime import datetime

import numpy as np

from core.projections import projections
from core.projections.converter import Converter
from core.projections.rotor import Rotor
from core.projections.translator import Translator


class ListParser(object):
    def __init__(self, val_type=None, separator=';'):
        self.separator = separator
        self.val_type = val_type

    def __call__(self, *args, **kwargs):
        str_arr = args[0].split(self.separator)

        if self.val_type:
            return [self.val_type(val) for val in str_arr if val]
        else:
            return str_arr


class PairParser(object):
    def __init__(self, key_type=None, val_type=None, separator=':'):
        self.separator = separator
        self.key_type = key_type
        self.val_type = val_type

    def __call__(self, *args, **kwargs):
        arr = args[0].split(self.separator)

        key_str = arr[0]
        val_str = self.separator.join(arr[1:])

        key = key_str if self.key_type is None else self.key_type(key_str)
        val = val_str if self.val_type is None else self.val_type(val_str)

        return key, val


class DateTimeParser(object):
    def __init__(self, fmt='%Y%m%d%H%M%S'):
        self.fmt = fmt

    def __call__(self, *args, **kwargs):
        return datetime.strptime(args[0], self.fmt)


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


def parse_slice(string):
    elem_arr_str = string.split(':')
    if len(elem_arr_str) > 3:
        raise Exception('Argument is not a valid slice.')

    try:
        slice_params = map(lambda x: int(x.strip()) if x.strip() else None,
                           elem_arr_str)
    except ValueError:
        raise Exception('Argument is not a valid slice.')

    if len(slice_params) == 1:
        if slice_params[0] is None:
            raise Exception('Argument is not a valid slice.')
        else:
            return slice(slice_params[0], slice_params[0] + 1, None)
    else:
        return slice(*slice_params)


def init_converter_from_args(args):
    proj_cls = projections[args.proj_name]

    p = proj_cls.unified_init(args.earth_radius, args.true_scale_lats)

    ref_lat, ref_lon = p.reference_point
    r = Rotor.build_rotor(args.orig_lat, args.orig_lon,
                          ref_lat, ref_lon,
                          args.adjust_angle)

    t = Translator(args.easting, args.northing)

    return Converter(r, p, t)


def split_scalar_and_vector_vars(var_names):
    scalar_vars = []
    vector_vars = []

    if var_names:
        for var_name_list in var_names:
            var_names = var_name_list.split('+')
            if len(var_names) == 1:
                scalar_vars.extend(var_names)
            elif len(var_names) == 2:
                vector_vars.append(var_names)
            elif len(var_names) > 2:
                raise Exception()

    return scalar_vars, vector_vars
