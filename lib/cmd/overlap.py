import numpy as np
import shutil
from netCDF4 import Dataset

from cmd.common.nc_utils import add_or_append_history
from cmd.common.arg_processors import ListParser

description = 'copies the first file in the list and fills masked values of ' \
              'a given field with values from fields of the following files'


def setup_parser(parser):
    list_parser = ListParser()
    parser.add_argument('--input-files', type=list_parser,
                        required=True)
    parser.add_argument('--output-file', required=True)
    parser.add_argument('--data-var-names', type=list_parser,
                        required=True)


def cmd(args):
    if len(args.data_var_names) == 1:
        data_var_names = ([args.data_var_names[0]] *
                          len(args.input_files))
    elif len(args.data_var_names) != len(args.input_files):
        raise Exception()
    else:
        data_var_names = args.data_var_names

    shutil.copyfile(args.input_files[0], args.output_file)

    output_ds = Dataset(args.output_file, 'r+')
    output_var = output_ds.variables[args.data_var_names[0]]
    output_var_data = output_var[:]

    for input_var_name, input_file in zip(data_var_names[1:],
                                          args.input_files[1:]):
        input_ds = Dataset(input_file, 'r')
        input_var_data = input_ds.variables[input_var_name][:]

        select_mask = np.logical_and(output_var_data.mask,
                                     ~input_var_data.mask)
        output_var_data[select_mask] = input_var_data[select_mask]

        input_ds.close()

    output_var[:] = output_var_data[:]

    add_or_append_history(output_ds)

    output_ds.close()
