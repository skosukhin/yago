import numpy as np
import sys
from netCDF4 import Dataset

import cmd.name_constants as names
from cmd.common import init_converter_from_proj_var, copy_nc_attributes, \
    gen_hist_string
from core.converter import convert_points

description = 'projects geographical coordinates to a plane'


def setup_parser(parser):
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--grid-file', required=True)


def cmd(args):
    ds = Dataset(args.input_file, 'r+')

    # Check that our input file was preprocessed.
    if (not hasattr(ds, names.ATTR_HISTORY) or
                'arctic preproc' not in ds.getncattr(names.ATTR_HISTORY)):
        raise Exception()

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

    history_update = gen_hist_string()
    history = [history_update + '\n',
               ds.getncattr(names.ATTR_HISTORY)] if hasattr(
        ds, names.ATTR_HISTORY) else history_update
    ds.setncattr(names.ATTR_HISTORY, history)

    ds.close()


def _progress(row_num, row_count):
    print '%.2f%%' % (float(row_num) / float(row_count) * 100.0)
    sys.stdout.flush()
