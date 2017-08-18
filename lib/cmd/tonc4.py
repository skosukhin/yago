import sys

import numpy as np
from netCDF4 import Dataset

from cmd.common.misc import create_dir_for_file
from cmd.common.nc_utils import copy_nc_attributes, MAX_COPY_DIM_COUNT, \
    DimIterator, create_nc_var_like_other, add_or_append_history

description = 'convert files to NC4 format optionally applying compression'


def setup_parser(parser):
    mandatory_args = parser.add_argument_group('mandatory arguments')
    mandatory_args.add_argument('--input-file',
                                help='name of netcdf file that needs to be '
                                     'converted',
                                required=True)
    mandatory_args.add_argument('--output-file',
                                help='output filename',
                                required=True)

    compression_args = parser.add_argument_group('compression arguments')
    compression_args.add_argument('--no-compression',
                                  help='disables compression of the output '
                                       'file',
                                  action='store_false',
                                  dest='zlib')

    compression_args.add_argument('--complevel',
                                  help='sets compression level '
                                       '(default: \'%(default)s\')',
                                  default=4,
                                  type=np.intp,
                                  choices=np.arange(1, 10, dtype=np.intp))
    compression_args.add_argument('--no-shuffle',
                                  help='disables shuffle filter',
                                  action='store_false',
                                  dest='shuffle')
    compression_args.add_argument('--fletcher32',
                                  help='enables Fletcher32 HDF5 checksum '
                                       'algorithm',
                                  action='store_true')


def cmd(args):
    in_ds = Dataset(args.input_file, 'r')
    create_dir_for_file(args.output_file)
    out_ds = Dataset(args.output_file, 'w', format='NETCDF4')

    copy_nc_attributes(in_ds, out_ds)

    for in_dim in in_ds.dimensions.values():
        out_ds.createDimension(in_dim.name,
                               None if in_dim.isunlimited() else in_dim.size)

    all_args = vars(args)
    compression_args = {}
    for key in ['zlib', 'complevel', 'shuffle', 'fletcher32']:
        compression_args[key] = all_args[key]

    for in_var in in_ds.variables.values():
        print in_var.name
        out_var = create_nc_var_like_other(out_ds, in_var, **compression_args)

        iter_mask = np.ones((len(in_var.shape, )), dtype=bool)
        iter_mask[-MAX_COPY_DIM_COUNT:] = False
        dim_iterator = DimIterator(in_var.shape, None, iter_mask)
        write_op_count = len(dim_iterator)
        for write_op_num, slc in enumerate(dim_iterator.slice_tuples()):
            _progress(write_op_num, write_op_count)
            out_var[slc] = in_var[slc]

    add_or_append_history(out_ds)
    in_ds.close()
    out_ds.close()


def _progress(row_num, row_count):
    print '%.2f%%' % (float(row_num) / float(row_count) * 100.0)
    sys.stdout.flush()
