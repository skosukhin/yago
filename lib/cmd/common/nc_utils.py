import os
import sys
import time

import numpy as np
from netcdftime import utime

import cmd.common.name_constants as names
from core.grids.rectilinear import RectilinearGrid
from core.grids.structured import StructuredGrid
from core.projections import projections
from core.projections.converter import Converter
from core.projections.rotors import Rotor, RotorX, RotorY, RotorZ

# Maximum number of dimensions to copy during one read/write operation
from core.projections.translator import Translator

MAX_COPY_DIM_COUNT = 2


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


def copy_nc_attributes(src_var, dst_var):
    for attr_name in src_var.ncattrs():
        dst_var.setncattr(attr_name, src_var.getncattr(attr_name))


def create_nc_var_like_other(out_ds, other_var, **overrides):
    kwargs = {'dimensions': other_var.dimensions}

    other_filters = other_var.filters()

    if other_filters is not None:
        kwargs.update(other_filters)

    kwargs.update(overrides)

    result = out_ds.createVariable(other_var.name, other_var.dtype, **kwargs)

    copy_nc_attributes(other_var, result)

    return result


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


def init_converter_from_proj_var(proj_var):
    r = _decode_rotor(proj_var.rot_axes, proj_var.rot_angles_deg)

    standard_parallel_list = np.atleast_1d(proj_var.standard_parallel)
    p = projections[proj_var.short_name].unified_init(proj_var.earth_radius,
                                                      standard_parallel_list)

    t = Translator(proj_var.false_easting, proj_var.false_northing)

    return Converter(r, p, t)


def find_dim_indices(names_to_find, dimensions):
    result = [None] * len(names_to_find)
    for result_idx, name_to_find in enumerate(names_to_find):
        for dim_idx, dim_name in enumerate(dimensions):
            if name_to_find == dim_name:
                if result[result_idx] is not None:
                    raise Exception()
                result[result_idx] = dim_idx
    return result


def reorder_axes(arr, new_order):
    current_order = range(len(arr.shape))
    undo_order = [None] * len(new_order)
    for current_idx, new_idx in enumerate(new_order):
        undo_order[new_idx] = current_idx
    return np.moveaxis(arr, current_order, new_order), undo_order


def rename_dimensions(dim_name_tuple, rename_dict):
    result = []
    for dim_name in dim_name_tuple:
        if dim_name in rename_dict:
            result.append(rename_dict[dim_name])
        else:
            result.append(dim_name)
    return tuple(result)


def add_missing_dim_vars(src_ds, dst_ds, dim_names):
    for dim_name in dim_names:
        if dim_name not in src_ds.dimensions:
            continue
        if dim_name not in dst_ds.dimensions:
            src_dim = src_ds.dimensions[dim_name]
            dst_ds.createDimension(dim_name,
                                   None if src_dim.isunlimited()
                                   else src_dim.size)
            if dim_name in src_ds.variables:
                src_dim_var = src_ds.variables[dim_name]
                if 1 == len(src_dim_var.dimensions) \
                        and src_dim_var.dimensions[0] == dim_name:
                    dst_dim_var = \
                        dst_ds.createVariable(
                            dim_name,
                            src_dim_var.dtype,
                            dimensions=src_dim_var.dimensions)
                    copy_nc_attributes(src_dim_var, dst_dim_var)
                    dst_dim_var[:] = src_dim_var[:]


def init_grid_from_vars(x_var, y_var):
    if len(x_var.shape) != len(y_var.shape) or len(x_var.shape) > 2:
        raise Exception()

    if len(x_var.shape) == 2:
        if x_var.dimensions != y_var.dimensions:
            raise Exception()
        return StructuredGrid(xx=x_var[:], yy=y_var[:]), x_var.dimensions
    else:
        return RectilinearGrid(x=x_var[:], y=y_var[:]), \
               y_var.dimensions + x_var.dimensions


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
