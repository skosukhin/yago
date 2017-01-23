import sys
from itertools import izip

import numpy as np
from netCDF4 import Dataset

import cmd.common.name_constants as names
from cmd.common.misc import create_dir_for_file
from cmd.common.nc_utils import add_or_append_history, \
    find_dim_indices, rename_dimensions, add_missing_dim_vars, DimIterator, \
    reorder_axes
from cmd.common.arg_processors import ListParser

description = 'applies provided weights to perform interpolation'


def setup_parser(parser):
    mandatory_args = parser.add_argument_group('mandatory arguments')
    mandatory_args.add_argument('--input-file',
                                help='name of netcdf file that contains '
                                     'fields that need to be interpolated',
                                required=True)
    mandatory_args.add_argument('--output-file',
                                help='output filename',
                                required=True)
    mandatory_args.add_argument('--weight-file',
                                help='name of netcdf file that interpolation '
                                     'weights will be read from',
                                required=True)
    list_parser = ListParser()
    mandatory_args.add_argument('--var-names',
                                help='\'%s\'-separated list of names of '
                                     'netcdf variables that need to be '
                                     'interpolated'
                                     % list_parser.separator,
                                type=list_parser)


def cmd(args):
    weight_ds = Dataset(args.weight_file, 'r')
    expected_in_dims = tuple(weight_ds.variables[names.VAR_INPUT_DIMS][:])
    expected_in_shape = tuple(weight_ds.variables[names.VAR_INPUT_SHAPE][:])

    in_ds = Dataset(args.input_file, 'r')
    for dim_idx, dim_name in enumerate(expected_in_dims):
        if dim_name not in in_ds.dimensions or \
                        in_ds.dimensions[dim_name].size != \
                        expected_in_shape[dim_idx]:
            raise Exception('Weight file does not match the input file.')

    weight_var = weight_ds.variables[names.VAR_WEIGHTS]
    out_dims = weight_var.dimensions[:-1]
    weights = weight_var[:]
    indices = _split_and_squeeze(
        np.ma.filled(weight_ds.variables[names.VAR_INDICES][:], 0), -2)

    weight_ds.close()

    create_dir_for_file(args.output_file)
    out_ds = Dataset(args.output_file, 'w')
    dim_rename_dict = {}
    for dim_idx, dim_name in enumerate(out_dims):
        out_ds.createDimension(dim_name, weights.shape[dim_idx])
        dim_rename_dict[expected_in_dims[dim_idx]] = dim_name

    for var_name in args.var_names:
        print var_name
        in_var = in_ds.variables[var_name]
        in_var_dim_indices = find_dim_indices(expected_in_dims,
                                              in_var.dimensions)
        in_field_dim_indices = np.argsort(in_var_dim_indices)

        iter_mask = np.ones((len(in_var.shape, )), dtype=bool)
        for dim_idx in in_var_dim_indices:
            if dim_idx is None:
                raise Exception()
            iter_mask[dim_idx] = False

        out_dim_tuple = rename_dimensions(in_var.dimensions, dim_rename_dict)
        add_missing_dim_vars(in_ds, out_ds, out_dim_tuple)

        out_var = out_ds.createVariable(var_name, in_var.dtype,
                                        dimensions=out_dim_tuple)

        read_iter = DimIterator(in_var.shape, None, iter_mask)
        write_iter = DimIterator(out_var.shape, None, iter_mask)
        write_op_count = len(read_iter)
        for write_op_num, (read_slc, write_slc) in enumerate(
                izip(read_iter.slice_tuples(), write_iter.slice_tuples())):
            _progress(write_op_num, write_op_count)
            in_field = in_var[read_slc]

            in_field, undo_order = reorder_axes(in_field, in_field_dim_indices)

            in_field = in_field[indices]

            out_field = np.ma.masked_where(
                np.ma.count_masked(in_field, axis=2) > 0,
                np.ma.sum(in_field * weights, axis=-1), copy=False)

            out_field, _ = reorder_axes(out_field, undo_order)

            out_var[write_slc] = out_field
        _progress(write_op_count, write_op_count)

    in_ds.close()

    add_or_append_history(out_ds)
    out_ds.close()


def _split_and_squeeze(arr, axis):
    return [np.squeeze(sub_arr, axis=axis) for sub_arr in
            np.split(arr, arr.shape[axis], axis=axis)]


def _progress(row_num, row_count):
    print '%.2f%%' % (float(row_num) / float(row_count) * 100.0)
    sys.stdout.flush()
