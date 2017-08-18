import sys

import numpy as np
from netCDF4 import Dataset

from cmd.common.arg_processors import ListParser, split_scalar_and_vector_vars
from cmd.common.misc import create_dir_for_file
from cmd.common.nc_utils import add_history, get_history, find_dim_indices, \
    MAX_COPY_DIM_COUNT, DimIterator, add_missing_dim_vars, copy_nc_attributes
from core.common import gen_rot_matrices_rad, apply_rot_matrices

description = 'rotates vector fields by specified angle'


def setup_parser(parser):
    mandatory_args = parser.add_argument_group('mandatory arguments')
    mandatory_args.add_argument('--input-file',
                                help='name of netcdf file that contains '
                                     'vector fields that need to be rotated',
                                required=True)
    mandatory_args.add_argument('--output-file',
                                help='output filename',
                                required=True)

    string_list_parser = ListParser()
    mandatory_args.add_argument('--var-names',
                                help='\'%s\'-separated list of names of '
                                     'netcdf variables that will be rotated '
                                     'and/or saved to the output file; if a '
                                     'couple of variables contain vector '
                                     'field components than they should '
                                     'appear as a single entry in the list, '
                                     'joined with the symbol \'+\' '
                                     '(e.g. uwnd+vwnd); at least one vector '
                                     'field must be specified'
                                     % string_list_parser.separator,
                                type=string_list_parser,
                                required=True)
    mandatory_args.add_argument('--angle-name',
                                help='name of netcdf variable that contains '
                                     'angles of rotation (in radians)',
                                required=True)

    parser.add_argument('--angle-file',
                        help='name of netcdf file that contains angles of '
                             'rotation; if the parameter is not specified then'
                             ' the angles are expected in the input file')
    parser.add_argument('--dim-names',
                        help='\'%s\'-separated list of names of netcdf '
                             'dimensions that define subfields of the input '
                             'fields that need to be rotated; the dimensions '
                             'must define the same shape as the angle '
                             'variable has; if the parameter is not specified '
                             'then the names of dimensions of the angle '
                             'variable are used to identify the rotated '
                             'subfields'
                             % string_list_parser.separator,
                        type=string_list_parser)


def cmd(args):
    scalar_vars, vector_vars = split_scalar_and_vector_vars(args.var_names)

    if len(vector_vars) < 1:
        raise Exception('No input vector fields specified.')

    if args.angle_file is None:
        args.angle_file = args.input_file

    angle_ds = Dataset(args.angle_file, 'r')

    if args.angle_name not in angle_ds.variables:
        raise Exception('Angle file \'%s\' does not contain angle variable '
                        '\'%s\'' % (args.angle_file, args.angle_name))

    angle_var = angle_ds.variables[args.angle_name]
    angle_shape = angle_var.shape

    in_ds = Dataset(args.input_file, 'r')

    if args.dim_names is None:
        args.dim_names = angle_var.dimensions
    else:
        if len(args.dim_names) != len(angle_var.dimensions):
            raise Exception('Number of dimensions specified by the parameter '
                            '--dim-names is not equal to the number of '
                            'dimensions of the --angle-name variable.')
        for dim_idx, dim_name in enumerate(args.dim_names):
            if dim_name not in in_ds.dimensions:
                raise Exception('Input file does not contain dimension '
                                '\'%s\'.' % dim_name)
            if in_ds.dimensions[dim_name].size != angle_shape[dim_idx]:
                raise Exception('Size of the dimension \'%s\' of the input '
                                'file is not equal to the size of the '
                                'dimension \'%s\' of the angle file.'
                                % (dim_name, angle_var.dimensions[dim_idx]))

    rot_matrices = gen_rot_matrices_rad(angle_var[:])
    angle_ds.close()

    create_dir_for_file(args.output_file)
    out_ds = Dataset(args.output_file, 'w')

    if len(scalar_vars) > 0:
        print 'Processing scalar fields:'

    for var_name in scalar_vars:
        print var_name
        in_var = in_ds.variables[var_name]
        add_missing_dim_vars(in_ds, out_ds, in_var.dimensions)

        out_var = out_ds.createVariable(var_name,
                                        in_var.dtype,
                                        dimensions=in_var.dimensions)
        copy_nc_attributes(in_var, out_var)

        iter_mask = np.ones((len(in_var.shape, )), dtype=bool)
        iter_mask[-MAX_COPY_DIM_COUNT:] = False

        dim_iterator = DimIterator(in_var.shape, None, iter_mask)
        write_op_count = len(dim_iterator)
        for write_op_num, slc in enumerate(dim_iterator.slice_tuples()):
            _progress(write_op_num, write_op_count)
            out_var[slc] = in_var[slc]
        _progress(write_op_count, write_op_count)

    if len(vector_vars) > 0:
        print 'Processing vector fields:'

    angle_iter_mask = np.ones((len(args.dim_names),), dtype=bool)
    dims_per_operation = len(args.dim_names) \
        if len(args.dim_names) < MAX_COPY_DIM_COUNT \
        else MAX_COPY_DIM_COUNT
    angle_iter_mask[-dims_per_operation:] = False
    angle_dim_iterator = DimIterator(angle_shape, None, angle_iter_mask)

    angle_iter_count = len(angle_dim_iterator)

    for var_name_pair in vector_vars:
        print var_name_pair
        in_u_var = in_ds.variables[var_name_pair[0]]
        in_v_var = in_ds.variables[var_name_pair[1]]

        if in_u_var.dimensions != in_v_var.dimensions:
            raise Exception()

        angle_dim_indices = find_dim_indices(args.dim_names,
                                             in_u_var.dimensions)

        for idx, dim_name in enumerate(args.dim_names):
            if angle_dim_indices[idx] is None:
                raise Exception('Variable \'%s\' is not specified along '
                                'dimension \'%s\'.'
                                % (var_name_pair[0], dim_name))

        add_missing_dim_vars(in_ds, out_ds, in_u_var.dimensions)

        out_u_var = out_ds.createVariable(
            var_name_pair[0],
            in_u_var.dtype,
            dimensions=in_u_var.dimensions)

        out_v_var = out_ds.createVariable(
            var_name_pair[1],
            in_v_var.dtype,
            dimensions=in_u_var.dimensions)

        fixed_indices = angle_dim_indices[-dims_per_operation:]

        dim_var_order = np.arange(len(fixed_indices))
        dim_angle_order = np.argsort(fixed_indices)

        var_iter_mask = np.ones((len(in_u_var.shape, )), dtype=bool)
        var_iter_mask[fixed_indices] = False

        for angle_iter_num, angle_slc in \
                enumerate(angle_dim_iterator.slice_tuples()):
            var_fixed_slices = np.repeat(None, len(in_u_var.shape))
            var_fixed_slices[angle_dim_indices] = angle_slc
            var_iter = DimIterator(in_u_var.shape, var_fixed_slices,
                                   var_iter_mask)
            var_iter_count = len(var_iter)
            total_iter_count = angle_iter_count * var_iter_count
            for var_iter_num, var_slc in enumerate(var_iter.slice_tuples()):
                _progress(var_iter_num + angle_iter_num * var_iter_count,
                          total_iter_count)
                in_u_field = in_u_var[var_slc]
                in_v_field = in_v_var[var_slc]

                in_u_field = np.moveaxis(in_u_field, dim_var_order,
                                         dim_angle_order)
                in_v_field = np.moveaxis(in_v_field, dim_var_order,
                                         dim_angle_order)

                out_u_field, out_v_field = \
                    apply_rot_matrices(in_u_field, in_v_field,
                                       rot_matrices[(slice(None), slice(None))
                                                    + angle_slc])

                out_u_field = np.moveaxis(out_u_field, dim_angle_order,
                                          dim_var_order)

                out_v_field = np.moveaxis(out_v_field, dim_angle_order,
                                          dim_var_order)

                out_u_var[var_slc] = out_u_field
                out_v_var[var_slc] = out_v_field

        _progress(total_iter_count, total_iter_count)

    add_history(out_ds, get_history(in_ds))

    in_ds.close()
    out_ds.close()


def _progress(row_num, row_count):
    print '%.2f%%' % (float(row_num) / float(row_count) * 100.0)
    sys.stdout.flush()
