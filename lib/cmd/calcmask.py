import numpy as np
from netCDF4 import Dataset

import cmd.name_constants as names

from cmd.common import copy_dim_var

description = 'calculates ocean mask for the following input data adjustments'


def setup_parser(parser):
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--depth-data-var', required=True)
    parser.add_argument('--output-file', required=True)
    parser.add_argument('--depth-factor', type=np.float64,
                        default=np.float64(1.0))


def cmd(args):
    input_ds = Dataset(args.input_file, 'r')

    depth_var = input_ds.variables[args.depth_data_var]

    output_ds = Dataset(args.output_file, 'w')

    for dim_name in depth_var.dimensions:
        copy_dim_var(input_ds, output_ds, dim_name)

    depth_data = depth_var[:]
    depth_data *= args.depth_factor

    mask = np.ma.masked_where(depth_data <= 0.0, np.ones(depth_data.shape),
                              copy=False)

    output_var = output_ds.createVariable(names.VAR_MASK, 'u1',
                                          dimensions=depth_var.dimensions,
                                          fill_value=0)
    output_var[:] = mask

    input_ds.close()
    output_ds.close()
