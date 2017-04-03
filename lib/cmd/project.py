import sys

import numpy as np
from netCDF4 import Dataset

import cmd.common.name_constants as names
from cmd.common.arg_processors import ListParser, split_scalar_and_vector_vars
from cmd.common.misc import create_dir_for_file
from cmd.common.nc_utils import copy_nc_attributes, \
    init_converter_from_proj_var, DimIterator, add_history, get_history, \
    MAX_COPY_DIM_COUNT, find_dim_indices, add_missing_dim_vars, \
    init_grid_from_vars

description = 'projects fields specified on a rectilinear grid in ' \
              'geographical coordinates to a Cartesian plane'


def setup_parser(parser):
    mandatory_args = parser.add_argument_group('mandatory arguments')
    mandatory_args.add_argument('--input-file',
                                help='name of netcdf file that contains data '
                                     'that need to be projected',
                                required=True)
    mandatory_args.add_argument('--grid-file',
                                help='name of netcdf file that contains '
                                     'projection description',
                                required=True)
    mandatory_args.add_argument('--output-file',
                                help='output filename',
                                required=True)
    mandatory_args.add_argument('--lat-name',
                                help='name of a netcdf variable that '
                                     'contains latitudes',
                                required=True)
    mandatory_args.add_argument('--lon-name',
                                help='name of a netcdf variable that '
                                     'contains longitudes',
                                required=True)

    string_list_parser = ListParser()
    parser.add_argument('--var-names',
                        help='\'%s\'-separated list of names of netcdf '
                             'variables that will be projected and/or saved '
                             'to the output file; if the list is empty than '
                             'only lat/lon coordinates will be projected; if '
                             'a couple of variables contain vector field '
                             'components aligned with meridians and parallels '
                             'than they should appear as a single entry in '
                             'the list, joined with the symbol \'+\' '
                             '(e.g. uwnd+vwnd)'
                             % string_list_parser.separator,
                        type=string_list_parser)
    parser.add_argument('--x-name',
                        help='name to be given to the netcdf variable that '
                             'contains x-coordinates of the applied '
                             'projection',
                        default=names.DIMVAR_X)
    parser.add_argument('--y-name',
                        help='name to be given to the netcdf variable that '
                             'contains y-coordinates of the applied '
                             'projection',
                        default=names.DIMVAR_Y)


def cmd(args):
    in_ds = Dataset(args.input_file, 'r')

    in_grid, in_grid_dim_names = \
        init_grid_from_vars(in_ds.variables[args.lon_name],
                            in_ds.variables[args.lat_name])

    scalar_vars, vector_vars = split_scalar_and_vector_vars(args.var_names)

    create_dir_for_file(args.output_file)
    out_ds = Dataset(args.output_file, 'w')

    add_missing_dim_vars(in_ds, out_ds, in_grid_dim_names)

    grid_ds = Dataset(args.grid_file, 'r')
    grid_proj_var = grid_ds.variables[names.VAR_PROJECTION]
    converter = init_converter_from_proj_var(grid_proj_var)

    print 'Calculating coordinates of grid points:'
    xx, yy = [], []
    for i in xrange(in_grid.shape[0]):
        _progress(i, in_grid.shape[0])
        row_xx, row_yy = converter.convert_points(in_grid[i, :, 1],
                                                  in_grid[i, :, 0])
        xx.append(row_xx)
        yy.append(row_yy)
    xx = np.concatenate(xx)
    yy = np.concatenate(yy)
    _progress(in_grid.shape[0], in_grid.shape[0])

    out_proj_var = out_ds.createVariable(names.VAR_PROJECTION,
                                         grid_proj_var.dtype)
    copy_nc_attributes(grid_proj_var, out_proj_var)

    out_x_var = out_ds.createVariable(args.x_name,
                                      xx.dtype,
                                      dimensions=in_grid_dim_names)
    copy_nc_attributes(grid_ds.variables[names.DIMVAR_X], out_x_var)
    out_x_var[:, :] = xx

    out_y_var = out_ds.createVariable(args.y_name,
                                      yy.dtype,
                                      dimensions=in_grid_dim_names)
    copy_nc_attributes(grid_ds.variables[names.DIMVAR_Y], out_y_var)
    out_y_var[:, :] = yy

    grid_ds.close()

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

    for var_name_pair in vector_vars:
        print var_name_pair
        in_u_var = in_ds.variables[var_name_pair[0]]
        in_v_var = in_ds.variables[var_name_pair[1]]

        if in_u_var.dimensions != in_v_var.dimensions:
            raise Exception()

        grid_dim_indices = find_dim_indices(
            in_grid_dim_names,
            in_u_var.dimensions)

        if any(idx is None for idx in grid_dim_indices):
            raise Exception()

        add_missing_dim_vars(in_ds, out_ds, in_u_var.dimensions)

        out_x_var = out_ds.createVariable(
            '_'.join(var_name_pair) + '_' + args.x_name,
            in_u_var.dtype,
            dimensions=in_u_var.dimensions)
        out_x_var.projection = out_proj_var.grid_mapping_name

        out_y_var = out_ds.createVariable(
            '_'.join(var_name_pair) + '_' + args.y_name,
            in_v_var.dtype,
            dimensions=in_u_var.dimensions)
        out_y_var.projection = out_proj_var.grid_mapping_name

        swap_axes = grid_dim_indices[0] > grid_dim_indices[1]

        iter_mask = np.ones((len(in_u_var.shape, )), dtype=bool)
        iter_mask[grid_dim_indices] = False

        dim_iterator = DimIterator(in_u_var.shape, None, iter_mask)
        write_op_count = len(dim_iterator)
        for write_op_num, slc in enumerate(dim_iterator.slice_tuples()):
            _progress(write_op_num, write_op_count)

            in_u_field = in_u_var[slc]
            in_v_field = in_v_var[slc]

            if swap_axes:
                in_u_field = np.swapaxes(in_u_field, grid_dim_indices[0],
                                         grid_dim_indices[1])
                in_v_field = np.swapaxes(in_v_field, grid_dim_indices[0],
                                         grid_dim_indices[1])

            out_x_field, out_y_field, = \
                converter.convert_vectors(in_u_field, in_v_field,
                                          in_grid[..., 1], in_grid[..., 0])

            if swap_axes:
                out_x_field = np.swapaxes(out_x_field, grid_dim_indices[0],
                                          grid_dim_indices[1])
                out_y_field = np.swapaxes(out_y_field, grid_dim_indices[0],
                                          grid_dim_indices[1])

            out_x_var[slc] = out_x_field
            out_y_var[slc] = out_y_field
        _progress(write_op_count, write_op_count)

    add_history(out_ds, get_history(in_ds))

    in_ds.close()
    out_ds.close()


def _progress(row_num, row_count):
    print '%.2f%%' % (float(row_num) / float(row_count) * 100.0)
    sys.stdout.flush()
