import numpy as np
import sys
from netCDF4 import Dataset

from cmd.common import init_converter_from_proj_var, copy_nc_attributes, \
    gen_hist_string
from core.converter import convert_points

description = 'projects geographical coordinates to a plane'


def setup_parser(parser):
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--grid-file', required=True)
    parser.add_argument('--output-file', required=True)
    parser.add_argument('--lat-var-name', required=True)
    parser.add_argument('--lon-var-name', required=True)


def cmd(args):
    input_ds = Dataset(args.input_file, 'r')
    input_lat_var = input_ds.variables[args.lat_var_name]
    input_lat_list = input_lat_var[:]
    input_lon_var = input_ds.variables[args.lon_var_name]
    input_lon_list = input_lon_var[:]

    grid_ds = Dataset(args.grid_file, 'r')
    grid_proj_var = grid_ds.variables['projection']
    converter = init_converter_from_proj_var(grid_proj_var)

    output_ds = Dataset(args.output_file, 'w')

    output_proj_var = \
        output_ds.createVariable(grid_proj_var.name, grid_proj_var.datatype)
    copy_nc_attributes(grid_proj_var, output_proj_var)

    output_ds.createDimension('lat', len(input_lat_list))
    output_lat_var = output_ds.createVariable('lat',
                                              input_lat_list.dtype,
                                              dimensions=('lat',))
    copy_nc_attributes(grid_ds.variables['lat'], output_lat_var)
    output_lat_var[:] = input_lat_list

    output_ds.createDimension('lon', len(input_lon_list))
    output_lon_var = output_ds.createVariable('lon',
                                              input_lon_list.dtype,
                                              dimensions=('lon',))
    copy_nc_attributes(grid_ds.variables['lon'], output_lon_var)
    output_lon_var[:] = input_lon_list

    lo, la = np.meshgrid(input_lon_list, input_lat_list)
    xx, yy = convert_points(la, lo, converter, _progress)

    output_x_var = output_ds.createVariable('x', xx.dtype,
                                            dimensions=('lat',
                                                        'lon'))
    copy_nc_attributes(grid_ds.variables['x'], output_x_var)
    output_x_var[:, :] = xx

    output_y_var = output_ds.createVariable('y', yy.dtype,
                                            dimensions=('lat',
                                                        'lon'))
    copy_nc_attributes(grid_ds.variables['y'], output_y_var)
    output_y_var[:, :] = yy

    output_ds.history = gen_hist_string()

    input_ds.close()
    grid_ds.close()
    output_ds.close()


def _progress(row_num, row_count):
    print '%.2f%%' % (float(row_num) / float(row_count) * 100.0)
    sys.stdout.flush()
