import sys

import numpy as np
from netCDF4 import Dataset

import cmd.name_constants as names

from cmd.common import check_preprocessed, init_converter_from_proj_var, \
    copy_nc_attributes, add_or_append_history
from core.converter import convert_points
from core.interpolation import RegularGrid, calc_weights

description = 'calculates weights for the following interpolation'


def setup_parser(parser):
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--output-file', required=True)
    parser.add_argument('--grid-file', required=True)
    parser.add_argument('--assume-lon-cycle', type=bool, required=True)


def cmd(args):
    input_ds = Dataset(args.input_file, 'r')

    check_preprocessed(input_ds)

    converter = None
    if names.VAR_PROJECTION in input_ds.variables:
        converter = init_converter_from_proj_var(
            input_ds.variables[names.VAR_PROJECTION])

    grid_ds = Dataset(args.grid_file, 'r')
    grid_proj_var = grid_ds.variables[names.VAR_PROJECTION]
    grid_converter = init_converter_from_proj_var(grid_proj_var)

    if converter:
        if converter != grid_converter:
            raise Exception('Input and grid files have different projection '
                            'descriptions.')
        input_xx = input_ds.variables[names.DIMVAR_X][:]
        input_yy = input_ds.variables[names.DIMVAR_Y][:]
    else:
        print 'Input file does not contain projection coordinates. We have ' \
              'to calculate them.'
        converter = grid_converter
        lo, la = np.meshgrid(input_ds.variables[names.DIMVAR_LON][:],
                             input_ds.variables[names.DIMVAR_LAT][:])
        input_xx, input_yy = convert_points(la, lo, converter, _progress)
        print 'Finished calculation of the projection coordinates.'

    input_ds.close()

    grid_x_var = grid_ds.variables[names.DIMVAR_X]
    grid_y_var = grid_ds.variables[names.DIMVAR_Y]

    grid = RegularGrid(grid_x_var[0], grid_ds.dimensions[names.DIMVAR_X].size,
                       grid_x_var.step,
                       grid_y_var[0], grid_ds.dimensions[names.DIMVAR_Y].size,
                       grid_y_var.step)

    quad_indices, weights = calc_weights(input_xx, input_yy,
                                         args.assume_lon_cycle, grid, _progress)

    output_ds = Dataset(args.output_file, 'w')

    grid_proj_var = grid_ds.variables[names.VAR_PROJECTION]
    output_proj_var = output_ds.createVariable(names.VAR_PROJECTION,
                                               grid_proj_var.dtype)
    copy_nc_attributes(grid_proj_var, output_proj_var)

    output_ds.createDimension(names.DIMVAR_X, grid_x_var.size)
    output_x_var = output_ds.createVariable(names.DIMVAR_X,
                                            grid_x_var.dtype,
                                            dimensions=(names.DIMVAR_X,))
    copy_nc_attributes(grid_x_var, output_x_var)
    output_x_var[:] = grid_x_var[:]

    output_ds.createDimension(names.DIMVAR_Y, grid_y_var.size)
    output_y_var = output_ds.createVariable(names.DIMVAR_Y,
                                            grid_y_var.dtype,
                                            dimensions=(names.DIMVAR_Y,))
    copy_nc_attributes(grid_y_var, output_y_var)
    output_y_var[:] = grid_y_var[:]

    grid_ds.close()

    output_ds.createDimension(names.DIM_DUO, 2)
    output_ds.createDimension(names.DIM_QUAD, 4)

    output_quad_var = output_ds.createVariable(names.VAR_INDICES,
                                               quad_indices.dtype,
                                               dimensions=(names.DIMVAR_Y,
                                                           names.DIMVAR_X,
                                                           names.DIM_DUO,
                                                           names.DIM_QUAD))
    output_quad_var[:] = quad_indices

    output_weight_var = output_ds.createVariable(names.VAR_WEIGHTS,
                                                 weights.dtype,
                                                 dimensions=(names.DIMVAR_Y,
                                                             names.DIMVAR_X,
                                                             names.DIM_QUAD))
    output_weight_var[:] = weights

    add_or_append_history(output_ds)

    output_ds.close()


def _progress(row_num, row_count):
    print '%.2f%%' % (float(row_num) / float(row_count) * 100.0)
    sys.stdout.flush()
