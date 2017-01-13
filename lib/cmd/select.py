import numpy as np
from netCDF4 import Dataset

from cmd.common import parse_slice, ListParser, DateTimeParser, PairParser, \
    get_history, add_history

description = 'selects subfields from input file'


def setup_parser(parser):
    mandatory_args = parser.add_argument_group('mandatory arguments')
    mandatory_args.add_argument('--input-file',
                                help='name of netcdf file that contains data '
                                     'that need to be projected',
                                required=True)
    mandatory_args.add_argument('--output-file',
                                help='output filename',
                                required=True)
    string_list_parser = ListParser()
    parser.add_argument('--var-names',
                        help='\'%s\'-separated list of names of netcdf '
                             'variables that will be saved '
                             'to the output file'
                             % string_list_parser.separator,
                        type=string_list_parser)
    string_slice_parser = PairParser(None, parse_slice)
    parser.add_argument('--slice',
                        help='\'%s\'-separated pair of strings where the '
                             'first string is a name of a netcdf dimension '
                             'and the second one is a python style slicing '
                             'string that will be applied to select a subset '
                             'of values along the dimension (e.g. if '
                             'dimension\'s name is \'lat\' and all but last '
                             'values along this dimension are required than '
                             'the value of this argument should be '
                             '\'lat::-1\'); this argument can be specified '
                             'multiple times to slice multiple dimensions'
                             % string_slice_parser.separator,
                        type=string_slice_parser,
                        action='append',
                        nargs='*')
    int_list_parser = ListParser(np.intp, ',')
    string_int_list_parser = PairParser(None, int_list_parser)
    parser.add_argument('--exclude-indices',
                        help='\'%s\'-separated pair of strings where the '
                             'first string is a name of a netcdf dimension '
                             'and the second one is a \'%s\'-separated list '
                             'of indices (starting from 0) that will NOT be '
                             'saved to the output file; this argument can be '
                             'specified multiple times to exclude indices '
                             'from multiple dimensions'
                             % (string_int_list_parser.separator,
                                int_list_parser.separator),
                        type=string_int_list_parser,
                        action='append',
                        nargs='*')
    string_int_parser = PairParser(None, np.int64)
    parser.add_argument('--min-int',
                        help='\'%s\'-separated pair of a string and an '
                             'integer number, where the string is a common '
                             'name for a netcdf dimension and a 1D netcdf '
                             'variable specified along this dimension, and '
                             'the number is the minimum threshold for the '
                             'variable\'s values that will be saved to the '
                             'output (indices that correspond to the values '
                             'that are less than the threshold will be '
                             'excluded from the dimension); this argument can '
                             'be specified multiple times to exclude indices '
                             'from multiple dimensions'
                             % string_int_parser.separator,
                        type=string_int_parser,
                        action='append',
                        nargs='*')
    parser.add_argument('--max-int',
                        help='\'%s\'-separated pair of a string and an '
                             'integer number, where the string is a common '
                             'name for a netcdf dimension and a 1D netcdf '
                             'variable specified along this dimension, and '
                             'the number is the maximum threshold for the '
                             'variable\'s values that will be saved to the '
                             'output (indices that correspond to the values '
                             'that are greater than the threshold will be '
                             'excluded from the dimension); this argument can '
                             'be specified multiple times to exclude indices '
                             'from multiple dimensions'
                             % string_int_parser.separator,
                        type=string_int_parser,
                        action='append',
                        nargs='*')
    string_float_parser = PairParser(None, np.float64)
    parser.add_argument('--min-float',
                        help='\'%s\'-separated pair of a string and a '
                             'floating point number, where the string is a '
                             'common name for a netcdf dimension and a 1D '
                             'netcdf variable specified along this dimension, '
                             'and the number is the minimum threshold for the '
                             'variable\'s values that will be saved to the '
                             'output (indices that correspond to the values '
                             'that are less than the threshold will be '
                             'excluded from the dimension); this argument can '
                             'be specified multiple times to exclude indices '
                             'from multiple dimensions'
                             % string_float_parser.separator,
                        type=string_float_parser,
                        action='append',
                        nargs='*')
    parser.add_argument('--max-float',
                        help='\'%s\'-separated pair of a string and a '
                             'floating point number, where the string is a '
                             'common name for a netcdf dimension and a 1D '
                             'netcdf variable specified along this dimension, '
                             'and the number is the maximum threshold for the '
                             'variable\'s values that will be saved to the '
                             'output (indices that correspond to the values '
                             'that are greater than the threshold will be '
                             'excluded from the dimension); this argument can '
                             'be specified multiple times to exclude indices '
                             'from multiple dimensions'
                             % string_float_parser.separator,
                        type=string_float_parser,
                        action='append',
                        nargs='*')
    time_parser = DateTimeParser()
    help_substring = time_parser.fmt.replace('%', '%%')
    string_time_parser = PairParser(None, time_parser)
    parser.add_argument('--min-time',
                        help='\'%s\'-separated pair of a string and a '
                             'timestamp in the format \'%s\', where the '
                             'string is a common name for a netcdf dimension '
                             'and a 1D netcdf variable specified along this '
                             'dimension, and the timestamp is the minimum '
                             'threshold for the variable\'s values that will '
                             'be saved to the output (indices that correspond '
                             'to the values that are less than the '
                             'threshold will be excluded from the dimension); '
                             'this argument can be specified multiple times '
                             'to exclude indices from multiple dimensions'
                             % (string_time_parser.separator, help_substring),
                        type=string_time_parser,
                        action='append',
                        nargs='*')
    parser.add_argument('--max-time',
                        help='\'%s\'-separated pair of a string and a '
                             'timestamp in the format \'%s\', where the '
                             'string is a common name for a netcdf dimension '
                             'and a 1D netcdf variable specified along this '
                             'dimension, and the timestamp is the maximum '
                             'threshold for the variable\'s values that will '
                             'be saved to the output (indices that correspond '
                             'to the values that are greater than the '
                             'threshold will be excluded from the dimension); '
                             'this argument can be specified multiple times '
                             'to exclude indices from multiple dimensions'
                             % (string_time_parser.separator, help_substring),
                        type=string_time_parser,
                        action='append',
                        nargs='*')


def cmd(args):
    in_ds = Dataset(args.input_file, 'r')
    out_ds = Dataset(args.output_file, 'w')

    modified_dimensions = {}

    for var_name in args.var_names:
        pass

    add_history(out_ds, get_history(in_ds))
    in_ds.close()
    out_ds.close()


def _get_dim_var_indices(ds, name, slice_, exclude, min_, max_):
    if name is None:
        return None
    indices = xrange(ds.dimensions[name].size)[slice_]
    if exclude:
        indices = [idx for idx in indices if idx not in exclude]
    val_list = None
    if min_:
        val_list = ds.variable[name][:]
        indices = [idx for idx in indices if val_list[idx] >= min_]
    if max_:
        if val_list is None:
            val_list = ds.variable[name][:]
        indices = [idx for idx in indices if val_list[idx] <= max_]

    return indices
