import shutil

from netCDF4 import Dataset

from cmd.common.misc import create_dir_for_file
from cmd.common.nc_utils import copy_nc_attributes, add_or_append_history
from cmd.common.arg_processors import ListParser

description = 'copies input file and adds variables from another one to it'


def setup_parser(parser):
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--output-file', required=True)
    parser.add_argument('--appended-file', required=True)

    parser.add_argument('--appended-var-names', type=ListParser(),
                        required=True)


def cmd(args):
    create_dir_for_file(args.output_file)
    shutil.copyfile(args.input_file, args.output_file)

    output_ds = Dataset(args.output_file, 'r+')
    appended_ds = Dataset(args.appended_file, 'r')

    for appended_var_name in args.appended_var_names:
        appended_var = appended_ds.variables[appended_var_name]
        output_var = output_ds.createVariable(
            appended_var_name, appended_var.dtype,
            dimensions=appended_var.dimensions)
        output_var[:] = appended_var[:]
        copy_nc_attributes(appended_var, output_var)

    add_or_append_history(output_ds)

    output_ds.close()
