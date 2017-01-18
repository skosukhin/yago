import sys

import numpy as np
from netCDF4 import Dataset

import cmd.common.name_constants as names
from cmd.common.misc import create_dir_for_file
from cmd.common.nc_utils import set_generic_lat_attributes, \
    set_generic_lon_attributes, copy_nc_attributes, add_history, get_history
from cmd.common.arg_processors import ListParser

description = 'applies a number of hacks to prepare input file for the ' \
              'following processing (deprecated)'


def setup_parser(parser):
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--output-file', required=True)
    parser.add_argument('--lat-var-name', required=True)
    parser.add_argument('--lon-var-name', required=True)
    parser.add_argument('--time-var-name')
    parser.add_argument('--data-var-names', type=ListParser(),
                        required=True)
    parser.add_argument('--ignore-last-lat', type=bool, default=False)
    parser.add_argument('--ignore-last-lon', type=bool, default=False)
    parser.add_argument('--cut-latitude', type=np.float64,
                        default=np.float64(0.0))


def cmd(args):
    input_ds = Dataset(args.input_file, 'r')

    input_lat_var = input_ds.variables[args.lat_var_name]
    input_lat_list = input_lat_var[:]
    output_lat_indices = np.arange(len(input_lat_list))
    if args.ignore_last_lat:
        output_lat_indices = output_lat_indices[:-1]

    output_lat_indices = output_lat_indices[
        np.where(input_lat_list[output_lat_indices] > args.cut_latitude)]

    input_lon_var = input_ds.variables[args.lon_var_name]
    input_lon_list = input_lon_var[:]
    output_lon_indices = np.arange(len(input_lon_list))
    if args.ignore_last_lon:
        output_lon_indices = output_lon_indices[:-1]

    output_lat_list = input_lat_list[output_lat_indices]

    create_dir_for_file(args.output_file)
    output_ds = Dataset(args.output_file, 'w')

    output_ds.createDimension(names.DIMVAR_LAT, len(output_lat_list))
    output_lat_var = output_ds.createVariable(names.DIMVAR_LAT,
                                              input_lat_list.dtype,
                                              dimensions=(names.DIMVAR_LAT,))
    set_generic_lat_attributes(output_lat_var)
    output_lat_var[:] = output_lat_list

    output_ds.createDimension(names.DIMVAR_LON, len(output_lon_indices))
    output_lon_var = output_ds.createVariable(names.DIMVAR_LON,
                                              input_lon_list.dtype,
                                              dimensions=(names.DIMVAR_LON,))
    set_generic_lon_attributes(output_lon_var)
    output_lon_var[:] = input_lon_list[output_lon_indices]

    input_nontemp_dim_tuple = (args.lat_var_name, args.lon_var_name)
    output_nontemp_dim_tuple = (names.DIMVAR_LAT, names.DIMVAR_LON)

    process_temp_vars = args.time_var_name is not None

    input_temp_dim_tuple = None
    output_temp_dim_tuple = None
    if process_temp_vars:
        input_time_var = input_ds.variables[args.time_var_name]
        input_time_list = input_time_var[:]
        output_ds.createDimension(names.DIMVAR_TIME, len(input_time_list))
        output_time_var = output_ds.createVariable(
            names.DIMVAR_TIME, input_time_list.dtype,
            dimensions=(names.DIMVAR_TIME,))
        copy_nc_attributes(input_time_var, output_time_var)
        output_time_var[:] = input_time_list

        input_temp_dim_tuple = (args.time_var_name,) + input_nontemp_dim_tuple
        output_temp_dim_tuple = (names.DIMVAR_TIME,) + output_nontemp_dim_tuple

    for input_var_name in args.data_var_names:
        print input_var_name
        input_var = input_ds.variables[input_var_name]

        # Yes, this can be simplified. But later.
        if input_var.dimensions != input_nontemp_dim_tuple:
            if not process_temp_vars:
                raise Exception()
            else:
                if input_var.dimensions != input_temp_dim_tuple:
                    raise Exception()

        out_dims, total_fields = \
            ((output_nontemp_dim_tuple, 1)
             if input_var.dimensions == input_nontemp_dim_tuple
             else (output_temp_dim_tuple, input_var.shape[0]))

        output_var = output_ds.createVariable(input_var_name,
                                              input_var.dtype,
                                              dimensions=out_dims)

        copy_nc_attributes(input_var, output_var)

        if input_var.dimensions == input_nontemp_dim_tuple:
            _progress(0, total_fields)
            output_var[:] = input_var[output_lat_indices, output_lon_indices]
        elif (process_temp_vars and
                input_var.dimensions == input_temp_dim_tuple):
            for time_idx in xrange(0, total_fields):
                _progress(time_idx, total_fields)
                output_var[time_idx, :] = input_var[
                    time_idx, output_lat_indices, output_lon_indices]
        else:
            raise Exception()

        _progress(total_fields, total_fields)

    add_history(output_ds, get_history(input_ds))

    input_ds.close()
    output_ds.close()


def _progress(row_num, row_count):
    print '%.2f%%' % (float(row_num) / float(row_count) * 100.0)
    sys.stdout.flush()
