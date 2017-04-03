import sys
from itertools import izip

import numpy as np
from netCDF4 import Dataset

from cmd.common.arg_processors import ListParser, split_scalar_and_vector_vars
from cmd.common.misc import create_dir_for_file
from cmd.common.nc_utils import copy_nc_attributes, add_missing_dim_vars, \
    find_dim_indices, DimIterator, add_history, get_history
from core.common import POLE_TOLERANCE, gen_rot_matrices_deg, apply_rot_matrices

description = 'extends input fields specified on a rectilinear lat/lon grid ' \
              'by adding grid points that correspond to the North (South) ' \
              'Pole; values for the new grid points are calculated as mean ' \
              'values along the highest (lowest) latitude of the fields'


def setup_parser(parser):
    mandatory_args = parser.add_argument_group('mandatory arguments')
    mandatory_args.add_argument('--input-file',
                                help='name of netcdf file that contains data '
                                     'that need to be modified',
                                required=True)
    mandatory_args.add_argument('--output-file',
                                help='output filename',
                                required=True)
    mandatory_args.add_argument('--lat-name',
                                help='name of 1D netcdf variable that '
                                     'contains latitudes',
                                required=True)

    parser.add_argument('--lon-name',
                        help='name of 1D netcdf variable that '
                             'contains longitudes')
    string_list_parser = ListParser()
    parser.add_argument('--var-names',
                        help='\'%s\'-separated list of names of netcdf '
                             'variables that will be extended and saved '
                             'to the output file; if the list is empty than '
                             'only lat/lon coordinates will be extended; if '
                             'a couple of variables contain vector field '
                             'components aligned with meridians and parallels '
                             'than they should appear as a single entry in '
                             'the list, joined with the symbol \'+\' '
                             '(e.g. uwnd+vwnd)'
                             % string_list_parser.separator,
                        type=string_list_parser)
    parser.add_argument('--add',
                        help='the Pole to add (default: %(default)s)',
                        choices=['north', 'south', 'both'], default='north')


def cmd(args):
    in_ds = Dataset(args.input_file, 'r')

    # Latitude variable is mandatory.
    in_lat_var = in_ds.variables[args.lat_name]
    if 1 != len(in_lat_var.dimensions):
        raise Exception('\'%s\' is not 1D variable.' % args.lat_name)

    in_lat_dim_name = in_lat_var.dimensions[0]

    # Longitude variable is optional (but only for scalar fields).
    if args.lon_name is not None:
        in_lon_var = in_ds.variables[args.lon_name]
        if 1 != len(in_lon_var.dimensions):
            raise Exception('\'%s\' is not 1D variable.' % args.lon_name)

        in_lon_dim_name = in_lon_var.dimensions[0]

        if in_lat_dim_name == in_lon_dim_name:
            raise Exception('Latitude and longitude dimension variables '
                            'can not be specified along the same dimension.')

    scalar_vars, vector_vars = split_scalar_and_vector_vars(args.var_names)

    if len(vector_vars) > 0 and args.lon_name is None:
        raise Exception('Vector fields cannot be processed without longitude '
                        'variable.')

    add_north_pole = (args.add == 'north' or args.add == 'both')
    add_south_pole = (args.add == 'south' or args.add == 'both')

    lat_list = in_lat_var[:]
    pole_tol = lat_list.dtype.type(POLE_TOLERANCE)
    np_lat = lat_list.dtype.type(90)

    min_lat_idx, max_lat_idx = 0, -1
    lat_list_ascending = True
    if lat_list.shape[0] > 1:
        if np.all(lat_list[1:] > lat_list[:-1]):
            pass
        elif np.all(lat_list[1:] < lat_list[:-1]):
            min_lat_idx, max_lat_idx = max_lat_idx, min_lat_idx
            lat_list_ascending = False
        else:
            raise Exception('Latitudes must be sorted and not be duplicated.')
    elif lat_list.shape[0] != 1:
        raise Exception('Array of latitudes must not be empty.')

    append_row = prepend_row = False
    if add_north_pole:
        if np.abs(lat_list[max_lat_idx] - np_lat) <= pole_tol:
            raise Exception('Input grid already contains grid points for the '
                            'North Pole.')
        if lat_list_ascending:
            append_row = True
        else:
            prepend_row = True

        lat_list = _extend_axis(lat_list, np_lat, not lat_list_ascending)

    if add_south_pole:
        if np.abs(lat_list[min_lat_idx] + np_lat) <= pole_tol:
            raise Exception('Input grid already contains grid points for the '
                            'South Pole.')
        if lat_list_ascending:
            prepend_row = True
        else:
            append_row = True

        lat_list = _extend_axis(lat_list, -np_lat, lat_list_ascending)

    create_dir_for_file(args.output_file)
    out_ds = Dataset(args.output_file, 'w')

    out_ds.createDimension(
        in_lat_dim_name,
        None if in_ds.dimensions[in_lat_dim_name].isunlimited()
        else lat_list.shape[0])
    out_lat_var = out_ds.createVariable(args.lat_name,
                                        in_lat_var.dtype,
                                        dimensions=(in_lat_dim_name,))
    copy_nc_attributes(in_lat_var, out_lat_var)
    out_lat_var[:] = lat_list

    if args.lon_name is not None:
        lon_list = in_lon_var[:]

        out_ds.createDimension(
            in_lon_dim_name,
            None if in_ds.dimensions[in_lon_dim_name].isunlimited()
            else lon_list.shape[0])
        out_lon_var = out_ds.createVariable(args.lon_name,
                                            in_lon_var.dtype,
                                            dimensions=(in_lon_dim_name,))
        copy_nc_attributes(in_lon_var, out_lon_var)
        out_lon_var[:] = lon_list

    if len(scalar_vars) > 0:
        print 'Processing scalar fields:'

    for var_name in scalar_vars:
        print var_name
        in_var = in_ds.variables[var_name]
        add_missing_dim_vars(in_ds, out_ds, in_var.dimensions)

        lat_idx, lon_idx = find_dim_indices(
            [args.lat_name, args.lon_name],
            in_var.dimensions)

        out_var = out_ds.createVariable(var_name,
                                        in_var.dtype,
                                        dimensions=in_var.dimensions)

        iter_mask = np.ones((len(in_var.shape, )), dtype=bool)

        if lat_idx is not None:
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

                out_field = in_field

                if prepend_row:
                    out_field = _extend_scalar_field(out_field, True)
                if append_row:
                    out_field = _extend_scalar_field(out_field, False)

                if swap_axes:
                    out_field = np.swapaxes(out_field, lat_idx, lon_idx)

                out_var[write_slc] = out_field

            _progress(write_op_count, write_op_count)

    if len(vector_vars) > 0:
        print 'Processing vector fields:'
        to_zero, from_zero = gen_rot_matrices_deg(lon_list, True)

    for var_name_pair in vector_vars:
        print var_name_pair

        in_u_var = in_ds.variables[var_name_pair[0]]
        in_v_var = in_ds.variables[var_name_pair[1]]

        if in_u_var.dimensions != in_v_var.dimensions:
            raise Exception()

        lat_idx, lon_idx = find_dim_indices(
            [args.lat_name, args.lon_name],
            in_u_var.dimensions)

        if lat_idx is None or lon_idx is None:
            raise Exception()

        add_missing_dim_vars(in_ds, out_ds, in_u_var.dimensions)

        out_u_var = out_ds.createVariable(var_name_pair[0],
                                          in_u_var.dtype,
                                          dimensions=in_u_var.dimensions)

        out_v_var = out_ds.createVariable(var_name_pair[1],
                                          in_v_var.dtype,
                                          dimensions=in_v_var.dimensions)

        swap_axes = lon_idx < lat_idx

        iter_mask = np.ones((len(in_u_var.shape, )), dtype=bool)
        iter_mask[lat_idx] = iter_mask[lon_idx] = False

        read_iter = DimIterator(in_u_var.shape, None, iter_mask)
        write_iter = DimIterator(out_u_var.shape, None, iter_mask)
        write_op_count = len(read_iter)
        for write_op_num, (read_slc, write_slc) in enumerate(
                izip(read_iter.slice_tuples(), write_iter.slice_tuples())):
            _progress(write_op_num, write_op_count)

            in_u_field = in_u_var[read_slc]
            in_v_field = in_v_var[read_slc]

            if swap_axes:
                in_u_field = np.swapaxes(in_u_field, lat_idx, lon_idx)
                in_v_field = np.swapaxes(in_v_field, lat_idx, lon_idx)

            out_u_field, out_v_field = in_u_field, in_v_field

            if prepend_row:
                out_u_field, out_v_field = \
                    _extend_vector_field(out_u_field, out_v_field,
                                         to_zero, from_zero,
                                         True)

            if append_row:
                out_u_field, out_v_field = \
                    _extend_vector_field(out_u_field, out_v_field,
                                         to_zero, from_zero,
                                         False)

            if swap_axes:
                out_u_field = np.swapaxes(out_u_field, lat_idx, lon_idx)
                out_v_field = np.swapaxes(out_v_field, lat_idx, lon_idx)

            out_u_var[write_slc] = out_u_field
            out_v_var[write_slc] = out_v_field

        _progress(write_op_count, write_op_count)

    add_history(out_ds, get_history(in_ds))

    in_ds.close()
    out_ds.close()


def _extend_axis(data, value, prepend):
    return np.insert(data, 0 if prepend else data.shape[0], value)


def _extend_scalar_field(data, prepend):
    dim_num = len(data.shape)

    if dim_num == 1:
        return np.ma.append(data[0], data) \
            if prepend else np.ma.append(data, data[-1])
    elif dim_num == 2:
        if prepend:
            pole_value = data.dtype.type(np.ma.mean(data[0]))
            return np.ma.append(
                [np.ma.repeat(pole_value, data.shape[1])], data, axis=0)
        else:
            pole_value = data.dtype.type(np.ma.mean(data[-1]))
            return np.ma.append(
                data, [np.ma.repeat(pole_value, data.shape[1])], axis=0)
    else:
        raise Exception()


def _extend_vector_field(u_data, v_data, to_zero_matrices, from_zero_matrices,
                         prepend):
    dim_num = len(u_data.shape)

    if dim_num == 2:
        row_idx = 0 if prepend else -1

        rot_u, rot_v = apply_rot_matrices(u_data[row_idx], v_data[row_idx],
                                          to_zero_matrices)
        u_pole_value = u_data.dtype.type(np.ma.mean(rot_u))
        v_pole_value = v_data.dtype.type(np.ma.mean(rot_v))
        new_u, new_v = \
            apply_rot_matrices(
                np.ma.repeat(u_pole_value, u_data.shape[1]),
                np.ma.repeat(v_pole_value, v_data.shape[1]),
                from_zero_matrices)

        if prepend:
            u_result = np.ma.append(
                [new_u], u_data, axis=0)
            v_result = np.ma.append(
                [new_v], v_data, axis=0)
        else:
            u_result = np.ma.append(
                u_data, [new_u], axis=0)
            v_result = np.ma.append(
                v_data, [new_v], axis=0)

        return u_result, v_result
    else:
        raise Exception()


def _progress(row_num, row_count):
    print '%.2f%%' % (float(row_num) / float(row_count) * 100.0)
    sys.stdout.flush()
