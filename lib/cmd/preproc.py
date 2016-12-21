import numpy as np
from netCDF4 import Dataset

from cmd.common import copy_nc_attributes, gen_hist_string, \
    parse_list_of_strings

description = 'applies a number of hacks to prepare input file for the ' \
              'following processing'


def setup_parser(parser):
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--output-file', required=True)
    parser.add_argument('--lat-var-name', required=True)
    parser.add_argument('--lon-var-name', required=True)
    parser.add_argument('--time-var-name')
    parser.add_argument('--data-var-names', type=parse_list_of_strings,
                        required=True)
    parser.add_argument('--ignore-last-lat', type=bool, required=True)
    parser.add_argument('--ignore-last-lon', type=bool, required=True)
    parser.add_argument('--add-north-pole', type=bool, required=True)
    parser.add_argument('--cut-latitude', type=np.float64,
                        default=np.float64(50.0))


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

    prepend_north_pole = None
    if args.add_north_pole:
        prepend_north_pole = (input_lat_list[output_lat_indices[0]] >
                              input_lat_list[output_lat_indices[1]])

    output_lat_list = input_lat_list[output_lat_indices]
    if args.add_north_pole:
        if prepend_north_pole:
            output_lat_list = np.append([90.0], output_lat_list)
        else:
            output_lat_list = np.append(output_lat_list, [90.0])

    output_ds = Dataset(args.output_file, 'w')
    copy_nc_attributes(input_ds, output_ds)

    output_ds.createDimension(args.lat_var_name, len(output_lat_list))
    output_lat_var = output_ds.createVariable(args.lat_var_name,
                                              input_lat_list.dtype,
                                              dimensions=(args.lat_var_name,))
    copy_nc_attributes(input_lat_var, output_lat_var)
    output_lat_var[:] = output_lat_list

    output_ds.createDimension(args.lon_var_name, len(output_lon_indices))
    output_lon_var = output_ds.createVariable(args.lon_var_name,
                                              input_lon_list.dtype,
                                              dimensions=(args.lon_var_name,))
    copy_nc_attributes(input_lon_var, output_lon_var)
    output_lon_var[:] = input_lon_list[output_lon_indices]

    nontemp_dim_tuple = (args.lat_var_name, args.lon_var_name)
    temp_dim_tuple = None

    if args.time_var_name:
        input_time_var = input_ds.variables[args.time_var_name]
        input_time_list = input_time_var[:]
        output_ds.createDimension(args.time_var_name, len(input_time_list))
        output_time_var = output_ds.createVariable(args.time_var_name,
                                                   input_time_list.dtype,
                                                   dimensions=(
                                                       args.time_var_name,))
        copy_nc_attributes(input_time_var, output_time_var)
        output_time_var[:] = input_time_list

        temp_dim_tuple = (
            args.time_var_name, args.lat_var_name, args.lon_var_name)

    for input_var_name in args.data_var_names:
        input_var = input_ds.variables[input_var_name]

        if input_var.dimensions == nontemp_dim_tuple:
            output_var_data = input_var[output_lat_indices, output_lon_indices]
        elif temp_dim_tuple and input_var.dimensions == temp_dim_tuple:
            output_var_data = \
                input_var[:, output_lat_indices, output_lon_indices]
        else:
            raise Exception()

        if args.add_north_pole:
            output_var_data = _add_north_pole(output_var_data,
                                              prepend_north_pole)

        output_var = output_ds.createVariable(input_var_name,
                                              output_var_data.dtype,
                                              dimensions=input_var.dimensions)
        copy_nc_attributes(input_var, output_var)
        output_var[:] = output_var_data

    input_ds.close()

    history = output_ds.history if hasattr(input_ds, 'history') else None
    history_update = gen_hist_string()
    history = [history_update + '\n', history] if history else history_update
    output_ds.history = history

    output_ds.close()


def _add_north_pole(data, prepend):
    dim_num = len(data.shape)
    if dim_num == 2:
        if prepend:
            north_pole_value = data.dtype.type(np.mean(data[0]))
            return np.append([np.repeat(north_pole_value, len(data[0]))], data,
                             axis=0)
        else:
            north_pole_value = data.dtype.type(np.mean(data[-1]))
            return np.append(data,
                             [np.repeat(north_pole_value, len(data[-1]))],
                             axis=0)
    elif dim_num == 3:
        return np.stack(_add_north_pole(dim2, prepend) for dim2 in data[:])
    else:
        raise Exception()
