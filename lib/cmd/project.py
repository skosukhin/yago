import sys

import numpy as np
from netCDF4 import Dataset

import cmd.name_constants as names
from cmd.common import init_converter_from_proj_var, copy_nc_attributes, \
    check_preprocessed, get_history, add_history, copy_dim_var, ListParser
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

    parser.add_argument('--data-var-names',
                        help='semicolon-separated list of names of netcdf '
                             'variables that will be projected and/or saved '
                             'to the output file; if the list is empty than '
                             'only lat/lon coordinates will be projected; if '
                             'a couple of variables contain vector field '
                             'components aligned with meridians and parallels '
                             'than they should appear as a single entry in '
                             'the list joined by the symbol \'+\' '
                             '(e.g. uwnd+vwnd)',
                        type=ListParser())
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


def cmd(args):
    input_ds = Dataset(args.input_file, 'r')

    check_preprocessed(input_ds)

    scalar_vars = []
    vector_vars = []

    if args.data_var_names:
        for var_name_list in args.data_var_names:
            var_names = var_name_list.split(',')
            if len(var_names) == 1:
                scalar_vars.extend(var_names)
            elif len(var_names) == 2:
                vector_vars.append(var_names)
            elif len(var_names) > 2:
                raise Exception()

    output_ds = Dataset(args.output_file, 'w')

    # Add north pole.
    input_lat_var = input_ds.variables[names.DIMVAR_LAT]
    input_lat_list = input_lat_var[:]
    output_lat_list = input_lat_list

    prepend_north_pole = None
    if args.add_north_pole:
        prepend_north_pole = input_lat_list[0] > input_lat_list[1]
        if prepend_north_pole:
            output_lat_list = np.append([90.0], input_lat_list)
        else:
            output_lat_list = np.append(input_lat_list, [90.0])

    output_ds.createDimension(names.DIMVAR_LAT, output_lat_list.size)
    output_lat_var = output_ds.createVariable(names.DIMVAR_LAT,
                                              input_lat_var.dtype,
                                              dimensions=(names.DIMVAR_LAT,))
    output_lat_var[:] = output_lat_list
    copy_nc_attributes(input_lat_var, output_lat_var)

    # Copy longitudes.
    input_lon_list = copy_dim_var(input_ds, output_ds, names.DIMVAR_LON, True)
    output_lon_list = input_lon_list

    # Copy timestamps if any.
    if names.DIMVAR_TIME in input_ds.variables:
        copy_dim_var(input_ds, output_ds, names.DIMVAR_TIME)

    grid_ds = Dataset(args.grid_file, 'r')
    grid_proj_var = grid_ds.variables[names.VAR_PROJECTION]
    converter = init_converter_from_proj_var(grid_proj_var)

    output_lo, output_la = np.meshgrid(output_lon_list, output_lat_list)
    print 'Calculating coordinates of grid points:'
    xx, yy = convert_points(output_la, output_lo, converter, _progress)

    output_proj_var = output_ds.createVariable(names.VAR_PROJECTION,
                                               grid_proj_var.dtype)
    copy_nc_attributes(grid_proj_var, output_proj_var)

    output_x_var = output_ds.createVariable(names.DIMVAR_X, xx.dtype,
                                            dimensions=(names.DIMVAR_LAT,
                                                        names.DIMVAR_LON))
    copy_nc_attributes(grid_ds.variables[names.DIMVAR_X], output_x_var)
    output_x_var[:, :] = xx

    output_y_var = output_ds.createVariable(names.DIMVAR_Y, yy.dtype,
                                            dimensions=(names.DIMVAR_LAT,
                                                        names.DIMVAR_LON))
    copy_nc_attributes(grid_ds.variables[names.DIMVAR_Y], output_y_var)
    output_y_var[:, :] = yy

    grid_ds.close()

    nontemp_dim_tuple = (names.DIMVAR_LAT, names.DIMVAR_LON)
    temp_dim_tuple = (names.DIMVAR_TIME,) + nontemp_dim_tuple

    if len(scalar_vars) > 0:
        print 'Processing scalar fields:'

    for var_name in scalar_vars:
        print var_name
        input_var = input_ds.variables[var_name]
        if (input_var.dimensions != nontemp_dim_tuple and
                input_var.dimensions != temp_dim_tuple):
            raise Exception()
        output_var = output_ds.createVariable(var_name, input_var.dtype,
                                              dimensions=input_var.dimensions)

        copy_nc_attributes(input_var, output_var)

        total_fields = \
            (1
             if input_var.dimensions == nontemp_dim_tuple
             else input_var.shape[0])

        if input_var.dimensions == nontemp_dim_tuple:
            _progress(0, total_fields)
            if args.add_north_pole:
                output_var[:] = _add_north_pole(input_var[:],
                                                prepend_north_pole)
            else:
                output_var[:] = input_var[:]
            _progress(total_fields, total_fields)
        elif input_var.dimensions == temp_dim_tuple:
            for time_idx in xrange(0, total_fields):
                _progress(time_idx, total_fields)
                if args.add_north_pole:
                    output_var[time_idx, :] = _add_north_pole(
                        input_var[time_idx, :], prepend_north_pole)
                else:
                    output_var[time_idx, :] = input_var[:]
        else:
            raise Exception()

        _progress(total_fields, total_fields)

    if len(vector_vars) > 0:
        print 'Processing vector fields:'
        input_lo, input_la = np.meshgrid(input_lon_list, input_lat_list)

    for var_name_pair in vector_vars:
        print var_name_pair
        input_u_var = input_ds.variables[var_name_pair[0]]
        if (input_u_var.dimensions != nontemp_dim_tuple and
                input_u_var.dimensions != temp_dim_tuple):
            raise Exception()
        input_v_var = input_ds.variables[var_name_pair[1]]
        if input_u_var.dimensions != input_v_var.dimensions:
            raise Exception()

        input_u_data, input_v_data, total_fields = \
            ((input_u_var[:], input_v_var[:], 1)
             if input_u_var.dimensions == nontemp_dim_tuple
             else (input_u_var[0, :], input_v_var[0, :], input_u_var.shape[0]))

        _progress(0, total_fields)
        output_u_data, output_v_data = convert_vectors(
            input_la, input_lo, input_u_data, input_v_data,
            converter)

        if args.add_north_pole:
            output_u_data = _add_north_pole(output_u_data, prepend_north_pole)
            output_v_data = _add_north_pole(output_v_data, prepend_north_pole)

        output_u_var = output_ds.createVariable(
            var_name_pair[0],
            output_u_data.dtype,
            dimensions=input_u_var.dimensions)
        output_u_var.projection = output_proj_var.grid_mapping_name

        output_v_var = output_ds.createVariable(
            var_name_pair[1],
            output_v_data.dtype,
            dimensions=input_v_var.dimensions)
        output_v_var.projection = output_proj_var.grid_mapping_name

        if input_u_var.dimensions == nontemp_dim_tuple:
            output_u_var[:] = output_u_data
            output_v_var[:] = output_v_data
        elif input_u_var.dimensions == temp_dim_tuple:
            output_u_var[0, :] = output_u_data
            output_v_var[0, :] = output_v_data
            for time_idx in xrange(1, total_fields):
                _progress(time_idx, total_fields)
                input_u_data = input_u_var[time_idx, :]
                input_v_data = input_v_var[time_idx, :]
                output_u_data, output_v_data = convert_vectors(
                    input_la, input_lo, input_u_data, input_v_data,
                    converter)
                if args.add_north_pole:
                    output_u_data = _add_north_pole(output_u_data,
                                                    prepend_north_pole)
                    output_v_data = _add_north_pole(output_v_data,
                                                    prepend_north_pole)
                output_u_var[time_idx, :] = output_u_data
                output_v_var[time_idx, :] = output_v_data
        else:
            raise Exception()
        _progress(total_fields, total_fields)

    add_history(output_ds, get_history(input_ds))

    input_ds.close()


def _progress(row_num, row_count):
    print '%.2f%%' % (float(row_num) / float(row_count) * 100.0)
    sys.stdout.flush()


def _add_north_pole(data, prepend):
    dim_num = len(data.shape)
    if dim_num == 2:
        if prepend:
            north_pole_value = data.dtype.type(np.ma.mean(data[0]))
            return np.ma.append(
                [np.ma.repeat(north_pole_value, len(data[0]))], data, axis=0)
        else:
            north_pole_value = data.dtype.type(np.ma.mean(data[-1]))
            return np.ma.append(
                data, [np.ma.repeat(north_pole_value, len(data[-1]))], axis=0)
    else:
        raise Exception()
