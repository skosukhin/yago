import sys
from itertools import izip

import numpy as np
from netCDF4 import Dataset

import cmd.name_constants as names
from cmd.common import init_converter_from_proj_var, copy_nc_attributes, \
    get_history, add_history, ListParser, MAX_COPY_DIM_COUNT, DimIterator
from core.converter import convert_points, convert_vectors

description = 'projects geographical coordinates and corresponding ' \
              'fields to a plane'


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
                                help='name of 1D netcdf variable that '
                                     'contains latitudes',
                                required=True)
    mandatory_args.add_argument('--lon-name',
                                help='name of 1D netcdf variable that '
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
    parser.add_argument('--add-north-pole',
                        help='flag that enables inclusion of the North Pole '
                             'grid point into input data before applying the '
                             'projection; flag should be set to \'1\' if the '
                             'input file does not contain grid point for the '
                             'North Pole to avoid \'holes\' in the output '
                             'fields; values that correspond to the '
                             'introduced grid points are calculated as mean '
                             'values of the arrays that correspond to the '
                             'highest latitude of the input fields',
                        type=bool, default=False)
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

    in_lat_var = in_ds.variables[args.lat_name]
    if 1 != len(in_lat_var.dimensions):
        raise Exception('\'%s\' is not 1D variable.' % args.lat_name)

    in_lon_var = in_ds.variables[args.lon_name]
    if 1 != len(in_lon_var.dimensions):
        raise Exception('\'%s\' is not 1D variable.' % args.lon_name)

    in_lat_dim_name = in_lat_var.dimensions[0]
    in_lon_dim_name = in_lon_var.dimensions[0]

    if in_lat_dim_name == in_lon_dim_name:
        raise Exception('Latitude and longitude dimension variables cannot '
                        'be specified along the same dimension.')

    scalar_vars = []
    vector_vars = []

    if args.var_names:
        for var_name_list in args.var_names:
            var_names = var_name_list.split('+')
            if len(var_names) == 1:
                scalar_vars.extend(var_names)
            elif len(var_names) == 2:
                vector_vars.append(var_names)
            elif len(var_names) > 2:
                raise Exception()

    out_ds = Dataset(args.output_file, 'w')

    # Add north pole.
    in_lat_list = in_lat_var[:]
    in_lat_dim = in_ds.dimensions[in_lat_dim_name]
    out_lat_list = in_lat_list

    prepend_north_pole = None
    if args.add_north_pole:
        prepend_north_pole = in_lat_list[0] > in_lat_list[1]
        if prepend_north_pole:
            out_lat_list = np.append([90.0], in_lat_list)
        else:
            out_lat_list = np.append(in_lat_list, [90.0])

    out_ds.createDimension(in_lat_dim_name,
                           None if in_lat_dim.isunlimited()
                           else out_lat_list.size)
    out_lat_var = out_ds.createVariable(args.lat_name,
                                        in_lat_var.dtype,
                                        dimensions=(in_lat_dim_name,))
    copy_nc_attributes(in_lat_var, out_lat_var)
    out_lat_var[:] = out_lat_list

    # Copy longitudes.
    in_lon_list = in_lon_var[:]
    in_lon_dim = in_ds.dimensions[in_lon_dim_name]
    out_lon_list = in_lon_list
    out_ds.createDimension(in_lon_dim_name,
                           None if in_lon_dim.isunlimited()
                           else out_lon_list.size)
    out_lon_var = out_ds.createVariable(args.lon_name,
                                        in_lon_var.dtype,
                                        dimensions=(in_lon_dim_name,))
    copy_nc_attributes(in_lon_var, out_lon_var)
    out_lon_var[:] = out_lon_list

    grid_ds = Dataset(args.grid_file, 'r')
    grid_proj_var = grid_ds.variables[names.VAR_PROJECTION]
    converter = init_converter_from_proj_var(grid_proj_var)

    out_lo, out_la = np.meshgrid(out_lon_list, out_lat_list)
    print 'Calculating coordinates of grid points:'
    xx, yy = convert_points(out_la, out_lo, converter, _progress)

    out_proj_var = out_ds.createVariable(names.VAR_PROJECTION,
                                         grid_proj_var.dtype)
    copy_nc_attributes(grid_proj_var, out_proj_var)

    out_x_var = out_ds.createVariable(args.x_name,
                                      xx.dtype,
                                      dimensions=(in_lat_dim_name,
                                                  in_lon_dim_name))
    copy_nc_attributes(grid_ds.variables[names.DIMVAR_X], out_x_var)
    out_x_var[:, :] = xx

    out_y_var = out_ds.createVariable(args.y_name,
                                      yy.dtype,
                                      dimensions=(in_lat_dim_name,
                                                  in_lon_dim_name))
    copy_nc_attributes(grid_ds.variables[names.DIMVAR_Y], out_y_var)
    out_y_var[:, :] = yy

    grid_ds.close()

    if len(scalar_vars) > 0:
        print 'Processing scalar fields:'

    for var_name in scalar_vars:
        print var_name
        in_var = in_ds.variables[var_name]
        _add_missing_dim_vars(in_ds, out_ds, in_var.dimensions)

        lat_idx, lon_idx = _find_dim_indices([args.lat_name, args.lon_name],
                                             in_var.dimensions)

        out_var = out_ds.createVariable(var_name,
                                        in_var.dtype,
                                        dimensions=in_var.dimensions)

        iter_mask = np.ones((len(in_var.shape, )), dtype=bool)

        if args.add_north_pole and lat_idx is not None:
            # In this case we add values for the North Pole.
            iter_mask[lat_idx] = False
            swap_axes = False
            if lon_idx is not None:
                iter_mask[lon_idx] = False
                swap_axes = lon_idx < lat_idx
            read_iter = DimIterator(in_var.shape, None, iter_mask)
            write_iter = DimIterator(out_var.shape, None, iter_mask)
            write_op_count = len(read_iter)
            for write_op_num, (read_slc, write_slc) in enumerate(
                    izip(read_iter.slice_tuples(), write_iter.slice_tuples())):
                _progress(write_op_num, write_op_count)
                in_field = in_var[read_slc]
                if swap_axes:
                    in_field = np.swapaxes(in_field, lat_idx, lon_idx)
                out_field = _add_north_pole(in_field, prepend_north_pole)
                if swap_axes:
                    out_field = np.swapaxes(out_field, lat_idx, lon_idx)
                out_var[write_slc] = out_field
            _progress(write_op_count, write_op_count)
        else:
            # Otherwise, we just copy the variable.
            iter_mask[-MAX_COPY_DIM_COUNT:] = False
            dim_iterator = DimIterator(in_var.shape, None, iter_mask)
            write_op_count = len(dim_iterator)
            for write_op_num, slc in enumerate(dim_iterator.slice_tuples()):
                _progress(write_op_num, write_op_count)
                out_var[slc] = in_var[slc]
            _progress(write_op_count, write_op_count)

    if len(vector_vars) > 0:
        print 'Processing vector fields:'
        in_lo, in_la = np.meshgrid(in_lon_list, in_lat_list)

    for var_name_pair in vector_vars:
        print var_name_pair
        in_u_var = in_ds.variables[var_name_pair[0]]
        in_v_var = in_ds.variables[var_name_pair[1]]

        if in_u_var.dimensions != in_v_var.dimensions:
            raise Exception()

        lat_idx, lon_idx = _find_dim_indices(
            [args.lat_name, args.lon_name],
            in_u_var.dimensions)

        if lat_idx is None or lon_idx is None:
            raise Exception()

        swap_axes = lon_idx < lat_idx

        _add_missing_dim_vars(in_ds, out_ds, in_u_var.dimensions)

        iter_mask = np.ones((len(in_u_var.shape, )), dtype=bool)
        iter_mask[lat_idx] = False
        iter_mask[lon_idx] = False

        read_iter = DimIterator(in_u_var.shape, None, iter_mask)
        write_op_count = len(read_iter)
        for write_op_num, read_slc in enumerate(read_iter.slice_tuples()):
            _progress(write_op_num, write_op_count)

            in_u_field = in_u_var[read_slc]
            in_v_field = in_v_var[read_slc]

            if swap_axes:
                in_u_field = np.swapaxes(in_u_field, lat_idx, lon_idx)
                in_v_field = np.swapaxes(in_v_field, lat_idx, lon_idx)

            out_x_field, out_y_field = convert_vectors(
                in_la, in_lo, in_u_field, in_v_field, converter)

            if args.add_north_pole:
                out_x_field = _add_north_pole(out_x_field,
                                              prepend_north_pole)
                out_y_field = _add_north_pole(out_y_field,
                                              prepend_north_pole)

            if swap_axes:
                out_x_field = np.swapaxes(out_x_field, lat_idx, lon_idx)
                out_y_field = np.swapaxes(out_y_field, lat_idx, lon_idx)

            if write_op_num == 0:
                out_x_var = out_ds.createVariable(
                    '_'.join(var_name_pair) + '_' + args.x_name,
                    out_x_field.dtype,
                    dimensions=in_u_var.dimensions)
                out_x_var.projection = out_proj_var.grid_mapping_name

                out_y_var = out_ds.createVariable(
                    '_'.join(var_name_pair) + '_' + args.y_name,
                    out_y_field.dtype,
                    dimensions=in_u_var.dimensions)
                out_y_var.projection = out_proj_var.grid_mapping_name

                write_iter_gen = DimIterator(out_x_var.shape, None,
                                             iter_mask).slice_tuples()

            write_slc = write_iter_gen.next()

            out_x_var[write_slc] = out_x_field
            out_y_var[write_slc] = out_y_field
        _progress(write_op_count, write_op_count)

    add_history(out_ds, get_history(in_ds))

    in_ds.close()
    out_ds.close()


def _progress(row_num, row_count):
    print '%.2f%%' % (float(row_num) / float(row_count) * 100.0)
    sys.stdout.flush()


def _add_north_pole(data, prepend):
    dim_num = len(data.shape)
    if dim_num == 1:
        return np.ma.append(data[0], data) \
            if prepend else np.ma.append(data, data[-1])
    elif dim_num == 2:
        if prepend:
            north_pole_value = data.dtype.type(np.ma.mean(data[0]))
            return np.ma.append(
                [np.ma.repeat(north_pole_value, data.shape[1])], data, axis=0)
        else:
            north_pole_value = data.dtype.type(np.ma.mean(data[-1]))
            return np.ma.append(
                data, [np.ma.repeat(north_pole_value, data.shape[1])], axis=0)
    else:
        raise Exception()


def _add_missing_dim_vars(src_ds, dst_ds, dim_names):
    for dim_name in dim_names:
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


def _find_dim_indices(names_to_find, dimensions):
    result = [None] * len(names_to_find)
    for result_idx, name_to_find in enumerate(names_to_find):
        for dim_idx, dim_name in enumerate(dimensions):
            if name_to_find == dim_name:
                if result[result_idx] is not None:
                    raise Exception()
                result[result_idx] = dim_idx
    return result


def _swap_values_in_tuple(tup, idx1, idx2):
    l = list(tup)
    l[idx1], l[idx2] = l[idx2], l[idx1]
    return tuple(l)
