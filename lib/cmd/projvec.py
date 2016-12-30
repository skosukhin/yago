import sys

from netCDF4 import Dataset

import cmd.name_constants as names

from cmd.common import parse_list_of_strings, init_converter_from_proj_var, \
    add_or_append_history
from core.converter import convert_vectors

description = 'projects vector field to a plane'


def setup_parser(parser):
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--vector-var-names', type=parse_list_of_strings,
                        required=True)


def cmd(args):
    if len(args.vector_var_names) != 2:
        raise Exception()

    input_ds = Dataset(args.input_file, 'r+')

    converter = init_converter_from_proj_var(
        input_ds.variables[names.VAR_PROJECTION])

    input_lats = input_ds.variables[names.DIMVAR_LAT][:]
    input_lons = input_ds.variables[names.DIMVAR_LON][:]

    input_u_var = input_ds.variables[args.vector_var_names[0]]
    input_v_var = input_ds.variables[args.vector_var_names[1]]

    if input_u_var.dimensions != input_v_var.dimensions:
        raise Exception()

    # Hardcode to work with temporal vars only
    if (input_u_var.dimensions[0] != names.DIMVAR_TIME or
            len(input_u_var.dimensions) != 3):
        raise Exception()

    output_u_var = input_ds.createVariable('rot_' + args.vector_var_names[0],
                                           input_u_var.dtype,
                                           dimensions=input_u_var.dimensions)

    output_v_var = input_ds.createVariable('rot_' + args.vector_var_names[1],
                                           input_v_var.dtype,
                                           dimensions=input_v_var.dimensions)

    for time_idx in xrange(0, input_u_var.shape[0]):
        _progress(time_idx, input_u_var.shape[0])

        input_u_data = input_u_var[time_idx, :]
        input_v_data = input_v_var[time_idx, :]

        output_u_data, output_v_data = convert_vectors(
            input_lats, input_lons, input_u_data, input_v_data, converter)

        output_u_var[time_idx, :] = output_u_data
        output_v_var[time_idx, :] = output_v_data

    add_or_append_history(input_ds)

    input_ds.close()


def _progress(row_num, row_count):
    print '%.2f%%' % (float(row_num) / float(row_count) * 100.0)
    sys.stdout.flush()
