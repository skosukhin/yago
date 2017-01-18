import os
import sys
import time

import numpy as np
from datetime import datetime
from netcdftime import utime

import cmd.name_constants as names
from core.converter import Converter
from core.projections import projections
from core.rotors import Rotor, RotorZ, RotorY, RotorX


# Maximum number of dimensions to copy during one read/write operation
MAX_COPY_DIM_COUNT = 2


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


def get_time_converter(time_var):
    try:
        time_var.units
    except AttributeError:
        raise AttributeError('netcdf time variable is missing a \'units\' '
                             'attribute')
    return utime(time_var.units,
                 calendar=getattr(time_var, 'calendar', 'standard'))


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
    if not hasattr(dataset, names.ATTR_HISTORY) or 'arctic preproc' not in \
            dataset.getncattr(names.ATTR_HISTORY):
        raise Exception()


def init_converter_from_proj_var(proj_var):
    r = _decode_rotor(proj_var.rot_axes, proj_var.rot_angles_deg)

    standard_parallel_list = proj_var.standard_parallel
    try:
        standard_parallel_list[0]
    except:
        standard_parallel_list = [standard_parallel_list]

    p = projections[proj_var.short_name].init(proj_var.earth_radius,
                                              standard_parallel_list)

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


class DimIterator(object):
    def __init__(self, shape, slices=None, iter_mask=None):
        self._empty = False
        self._slices, self._iter_mask, self._index_lists = None, None, None

        if len(shape) == 0:
            self._empty = True
            return

        for s in shape:
            if s == 0:
                self._empty = True
                return

        if iter_mask is None:
            iter_mask = [True] * len(shape)
        else:
            if len(iter_mask) != len(shape):
                raise Exception()

        index_lists = [None] * len(shape)
        if slices is None:
            slices = [slice(0, stop, 1) for stop in shape]
        else:
            if len(slices) != len(shape):
                raise Exception()
            for i, s in enumerate(slices):
                if slices[i] is None:
                    slices[i] = slice(0, shape[i], 1)
                else:
                    s = slices[i]
                    if isinstance(s, slice):
                        s = slice(*s.indices(shape[i]))
                        if s.start == s.stop:
                            self._empty = True
                            return
                    else:
                        try:
                            index_list = list(s)
                            if len(index_list) == 0:
                                self._empty = True
                                return
                            index_lists[i] = index_list
                            s = slice(0, len(index_list), 1)
                        except TypeError:
                            index_lists[i] = [s]
                            s = slice(0, 1, 1)

                    slices[i] = s

        if not self._empty:
            self._slices = slices
            self._iter_mask = iter_mask
            self._index_lists = index_lists

    def __len__(self):
        if self._empty:
            return 0

        result = 1
        for idx, slc in enumerate(self._slices):
            if self._iter_mask[idx]:
                result *= len(xrange(slc.start, slc.stop, slc.step))
        return result

    def slice_tuples(self):
        if not self._empty:
            current_slice = [s.start if self._iter_mask[i] else s for i, s in
                             enumerate(self._slices)]

            while True:
                result = \
                    [s if self._index_lists[i] is None else
                     self._index_lists[i][s]
                     for i, s in enumerate(current_slice)]
                yield tuple(result)

                stop = True
                for idx in xrange(len(self._slices) - 1, -1, -1):
                    if self._iter_mask[idx]:
                        next_idx = current_slice[idx] + self._slices[idx].step
                        if (self._slices[idx].step > 0
                            and self._slices[idx].stop > next_idx) \
                                or (self._slices[idx].step < 0
                                    and self._slices[idx].stop < next_idx):
                            current_slice[idx] = next_idx
                            stop = False
                            break
                        else:
                            current_slice[idx] = self._slices[idx].start
                if stop:
                    break


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
