import sys

import numpy as np
from netCDF4 import Dataset

import cmd.name_constants as names
from cmd.common import init_converter_from_proj_var, copy_nc_attributes, \
    add_or_append_history, check_preprocessed
from core.converter import convert_points

description = 'projects geographical coordinates to a plane'


def setup_parser(parser):
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--grid-file', required=True)


def cmd(args):
    ds = Dataset(args.input_file, 'r+')

    check_preprocessed(ds)

    input_lat_list = ds.variables[names.DIMVAR_LAT][:]
    input_lon_list = ds.variables[names.DIMVAR_LON][:]

    grid_ds = Dataset(args.grid_file, 'r')
    grid_proj_var = grid_ds.variables[names.VAR_PROJECTION]
    converter = init_converter_from_proj_var(grid_proj_var)

    lo, la = np.meshgrid(input_lon_list, input_lat_list)
    xx, yy = convert_points(la, lo, converter, _progress)

    output_proj_var = ds.createVariable(names.VAR_PROJECTION,
                                        grid_proj_var.dtype)
    copy_nc_attributes(grid_proj_var, output_proj_var)

    output_x_var = ds.createVariable(names.DIMVAR_X, xx.dtype,
                                     dimensions=
                                     (names.DIMVAR_LAT, names.DIMVAR_LON))
    copy_nc_attributes(grid_ds.variables[names.DIMVAR_X], output_x_var)
    output_x_var[:, :] = xx

    output_y_var = ds.createVariable(names.DIMVAR_Y, yy.dtype,
                                     dimensions=
                                     (names.DIMVAR_LAT, names.DIMVAR_LON))
    copy_nc_attributes(grid_ds.variables[names.DIMVAR_Y], output_y_var)
    output_y_var[:, :] = yy

    grid_ds.close()

    add_or_append_history(ds)

    ds.close()


def _progress(row_num, row_count):
    print '%.2f%%' % (float(row_num) / float(row_count) * 100.0)
    sys.stdout.flush()
