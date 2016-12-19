import numpy as np
from netCDF4 import Dataset

from cmd.common import copy_nc_attributes, gen_hist_string

description = 'strips inputs'


def setup_parser(parser):
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--output-file', required=True)
    parser.add_argument('--lat-var-name', required=True)
    parser.add_argument('--lon-var-name', required=True)
    parser.add_argument('--time-var-name')
    parser.add_argument('--data-var-name', required=True)
    parser.add_argument('--ignore-last-lat', type=bool, required=True)
    parser.add_argument('--ignore-last-lon', type=bool, required=True)
    parser.add_argument('--cut-latitude', type=np.float64,
                        default=np.float64(50.0))


def cmd(args):
    input_ds = Dataset(args.input_file, 'r')

    input_data_var = input_ds.variables[args.data_var_name]

    if args.time_var_name:
        expected_dim_tuple = (
            args.time_var_name, args.lat_var_name, args.lon_var_name)
        input_time_var = input_ds.variables[args.time_var_name]
    else:
        expected_dim_tuple = (args.lat_var_name, args.lon_var_name)
        input_time_var = None

    if input_data_var.dimensions != expected_dim_tuple:
        raise Exception()

    input_lat_var = input_ds.variables[args.lat_var_name]
    input_lon_var = input_ds.variables[args.lon_var_name]

    field = _Field(input_data_var[:], input_lat_var[:], input_lon_var[:],
                   input_time_var[:] if input_time_var else None)

    if args.ignore_last_lat:
        field.exclude_last_lat()
    if args.ignore_last_lon:
        field.exclude_last_lon()
    field.cut_below_parallel(args.cut_latitude)

    output_ds = Dataset(args.output_file, 'w')
    copy_nc_attributes(input_ds, output_ds)

    output_ds.createDimension(args.lat_var_name, len(field.lat_list))
    output_lat_var = output_ds.createVariable(args.lat_var_name,
                                              field.lat_list.dtype,
                                              dimensions=(args.lat_var_name,))
    copy_nc_attributes(input_lat_var, output_lat_var)
    output_lat_var[:] = field.lat_list

    output_ds.createDimension(args.lon_var_name, len(field.lon_list))
    output_lon_var = output_ds.createVariable(args.lon_var_name,
                                              field.lon_list.dtype,
                                              dimensions=(args.lon_var_name,))
    copy_nc_attributes(input_lon_var, output_lon_var)
    output_lon_var[:] = field.lon_list

    if input_time_var:
        output_ds.createDimension(args.time_var_name, len(field.time_list))
        output_time_var = output_ds.createVariable(args.time_var_name,
                                                   field.time_list.dtype,
                                                   dimensions=(
                                                       args.time_var_name,))
        copy_nc_attributes(input_time_var, output_time_var)
        output_time_var[:] = field.time_list

        output_data_var = output_ds.createVariable(args.data_var_name,
                                                   field.data.dtype,
                                                   dimensions=(
                                                       args.time_var_name,
                                                       args.lat_var_name,
                                                       args.lon_var_name))
    else:
        output_data_var = output_ds.createVariable(args.data_var_name,
                                                   field.data.dtype,
                                                   dimensions=(
                                                       args.lat_var_name,
                                                       args.lon_var_name))

    copy_nc_attributes(input_data_var, output_data_var)
    output_data_var[:] = field.data

    history = input_ds.history if hasattr(input_ds, 'history') else None
    history_update = gen_hist_string()
    history = [history_update + '\n', history] if history else history_update
    output_ds.history = history

    input_ds.close()
    output_ds.close()


class _Field(object):
    def __init__(self, data, lat_list, lon_list, time_list):
        self.data = data
        self.lat_list = lat_list
        self.lon_list = lon_list
        self.time_list = time_list

    def exclude_last_lat(self):
        self.lat_list = self.lat_list[:-1]

        if self.time_list:
            self.data = self.data[:, :-1]
        else:
            self.data = self.data[:-1]

    def exclude_last_lon(self):
        self.lon_list = self.lon_list[:-1]

        if self.time_list:
            self.data = self.data[:, :, :-1]
        else:
            self.data = self.data[:, :-1]

    def cut_below_parallel(self, cut_lat):
        indices_to_keep = np.where(self.lat_list > cut_lat)
        self.lat_list = self.lat_list[indices_to_keep]

        if self.time_list:
            self.data = self.data[:, indices_to_keep]
        else:
            self.data = self.data[indices_to_keep]
