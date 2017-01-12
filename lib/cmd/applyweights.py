import sys
import numpy as np
from netCDF4 import Dataset

import cmd.name_constants as names
from cmd.common import copy_nc_attributes, add_or_append_history, ListParser
from core.interpolation import apply_weights

description = 'applies weights to perform interpolation'


def setup_parser(parser):
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--output-file', required=True)
    parser.add_argument('--weight-file', required=True)
    list_parser = ListParser()
    parser.add_argument('--data-var-names', type=list_parser,
                        required=True)
    parser.add_argument('--data-var-types', type=list_parser)


def cmd(args):
    if not args.data_var_types or len(args.data_var_types) < 1:
        data_var_types = [None] * len(args.data_var_names)
    elif len(args.data_var_types) == 1:
        data_var_types = ([np.dtype(args.data_var_types[0])] *
                          len(args.data_var_names))
    elif len(args.data_var_types) != len(args.data_var_names):
        raise Exception()
    else:
        data_var_types = args.data_var_types

    weight_ds = Dataset(args.weight_file, 'r')

    output_ds = Dataset(args.output_file, 'w')

    weight_proj_var = weight_ds.variables[names.VAR_PROJECTION]
    output_proj_var = output_ds.createVariable(names.VAR_PROJECTION,
                                               weight_proj_var.dtype)
    copy_nc_attributes(weight_proj_var, output_proj_var)

    weight_x_var = weight_ds.variables[names.DIMVAR_X]
    output_ds.createDimension(names.DIMVAR_X, weight_x_var.size)
    output_x_var = output_ds.createVariable(names.DIMVAR_X,
                                            weight_x_var.dtype,
                                            dimensions=(names.DIMVAR_X,))
    copy_nc_attributes(weight_x_var, output_x_var)
    output_x_var[:] = weight_x_var[:]

    weight_y_var = weight_ds.variables[names.DIMVAR_Y]
    output_ds.createDimension(names.DIMVAR_Y, weight_y_var.size)
    output_y_var = output_ds.createVariable(names.DIMVAR_Y,
                                            weight_y_var.dtype,
                                            dimensions=(names.DIMVAR_Y,))
    copy_nc_attributes(weight_y_var, output_y_var)
    output_y_var[:] = weight_y_var[:]

    weight_lat_var = weight_ds.variables[names.DIMVAR_LAT]
    output_lat_var = output_ds.createVariable(names.DIMVAR_LAT,
                                              weight_lat_var.dtype,
                                              dimensions=(names.DIMVAR_Y,
                                                          names.DIMVAR_X))
    copy_nc_attributes(weight_lat_var, output_lat_var)
    output_lat_var[:] = weight_lat_var[:]

    weight_lon_var = weight_ds.variables[names.DIMVAR_LON]
    output_lon_var = output_ds.createVariable(names.DIMVAR_LON,
                                              weight_lon_var.dtype,
                                              dimensions=(names.DIMVAR_Y,
                                                          names.DIMVAR_X))
    copy_nc_attributes(weight_lon_var, output_lon_var)
    output_lon_var[:] = weight_lon_var[:]

    weight_var = weight_ds.variables[names.VAR_WEIGHTS]
    input_nontemp_dim_tuple = (names.DIMVAR_LAT, names.DIMVAR_LON)
    output_nontemp_dim_tuple = weight_var.dimensions[:-1]

    input_ds = Dataset(args.input_file, 'r')

    input_temp_dim_tuple, output_temp_dim_tuple = None, None
    if (names.DIMVAR_TIME in input_ds.variables and
            names.DIMVAR_TIME in input_ds.dimensions):
        input_temp_dim_tuple = \
            (names.DIMVAR_TIME,) + input_nontemp_dim_tuple
        output_temp_dim_tuple = (names.DIMVAR_TIME,) + output_nontemp_dim_tuple

        input_time_var = input_ds.variables[names.DIMVAR_TIME]
        output_ds.createDimension(names.DIMVAR_TIME, input_time_var.size)
        output_time_var = output_ds.createVariable(names.DIMVAR_TIME,
                                                   input_time_var.dtype,
                                                   dimensions=(
                                                       names.DIMVAR_TIME,))
        copy_nc_attributes(input_time_var, output_time_var)
        output_time_var[:] = input_time_var[:]

    quad_indices = weight_ds.variables[names.VAR_INDICES][:]
    weights = weight_var[:]
    weight_ds.close()

    for i, data_var_name in enumerate(args.data_var_names):
        print data_var_name
        input_data_var = input_ds.variables[data_var_name]
        output_data_type = data_var_types[i]

        if input_data_var.dimensions == input_nontemp_dim_tuple:
            output_data = apply_weights(input_data_var[:], quad_indices,
                                        weights, output_data_type)
            output_data_var = output_ds.createVariable(
                data_var_name, output_data.dtype,
                dimensions=output_nontemp_dim_tuple)
            output_data_var[:] = output_data
        elif (input_temp_dim_tuple and
                input_data_var.dimensions == input_temp_dim_tuple):
            _progress(0, input_data_var.shape[0])
            output_data = apply_weights(input_data_var[0, :], quad_indices,
                                        weights, output_data_type)
            output_data_var = output_ds.createVariable(
                data_var_name, output_data.dtype,
                dimensions=output_temp_dim_tuple)
            output_data_var[0, :] = output_data

            for time_idx in xrange(1, input_data_var.shape[0]):
                _progress(time_idx, input_data_var.shape[0])
                output_data = apply_weights(input_data_var[time_idx, :],
                                            quad_indices,
                                            weights,
                                            output_data_type)
                output_data_var[time_idx, :] = output_data

            _progress(input_data_var.shape[0], input_data_var.shape[0])
        else:
            raise Exception()

    input_ds.close()

    add_or_append_history(output_ds)
    output_ds.close()


def _progress(row_num, row_count):
    print '%.2f%%' % (float(row_num) / float(row_count) * 100.0)
    sys.stdout.flush()
