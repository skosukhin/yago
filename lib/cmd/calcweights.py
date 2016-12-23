import sys

import numpy as np
from netCDF4 import Dataset

import cmd.name_constants as names

from cmd.common import check_preprocessed, init_converter_from_proj_var
from core.converter import convert_points

description = 'calculates weights for the following interpolation'


def setup_parser(parser):
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--output-file', required=True)
    parser.add_argument('--grid-file', required=True)


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

    input_xx, input_yy = None, None
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

    grid_xx = grid_ds.variables[names.DIMVAR_X][:]
    grid_yy = grid_ds.variables[names.DIMVAR_Y][:]


def _progress(row_num, row_count):
    print '%.2f%%' % (float(row_num) / float(row_count) * 100.0)
    sys.stdout.flush()
