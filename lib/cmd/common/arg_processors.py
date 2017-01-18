from datetime import datetime

import numpy as np

from core.converter import Converter
from core.projections import projections


class ListParser(object):
    def __init__(self, val_type=None, separator=';'):
        self.separator = separator
        self.val_type = val_type

    def __call__(self, *args, **kwargs):
        str_arr = args[0].split(self.separator)

        if self.val_type:
            return map(self.val_type, str_arr)
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
    p = projections[args.proj_name].init(args.earth_radius,
                                         args.true_scale_lats)
    r = p.build_rotor(args.orig_lat, args.orig_lon, args.adjust_angle)

    return Converter(r, p)
