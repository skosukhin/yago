import argparse
from itertools import izip
from datetime import datetime

import numpy as np
from netCDF4 import Dataset

from cmd.common.misc import create_dir_for_file
from cmd.common.nc_utils import copy_nc_attributes, MAX_COPY_DIM_COUNT, \
    DimIterator, add_history, get_history, get_time_converter
from cmd.common.arg_processors import ListParser, PairParser, parse_slice, \
    DateTimeParser

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
    mandatory_args.add_argument('--var-names',
                                help='\'%s\'-separated list of names of '
                                     'netcdf variables that will be saved to '
                                     'the output file'
                                     % string_list_parser.separator,
                                type=string_list_parser,
                                action=AddToSetAction,
                                required=True)
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
                        action=AddToDictAction,
                        dest='slice_dict')
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
                        action=AddToDictAction,
                        dest='exclude_dict')
    string_int_parser = PairParser(None, np.int64)
    parser.add_argument('--min-val',
                        help='\'%s\'-separated pair of a string and a value, '
                             'where the string is a common name for a netcdf '
                             'dimension and a 1D netcdf variable specified '
                             'along this dimension, and the value is the '
                             'minimum threshold for the variable\'s values '
                             'that will be saved to the output (indices that '
                             'correspond to the values that are less than the '
                             'threshold will be excluded from the dimension); '
                             'this argument can be specified multiple times '
                             'to exclude indices from multiple dimensions'
                             % string_int_parser.separator,
                        type=string_int_parser,
                        action=AddToDictAction,
                        dest='min_dict')
    parser.add_argument('--max-val',
                        help='\'%s\'-separated pair of a string and a value, '
                             'where the string is a common name for a netcdf '
                             'dimension and a 1D netcdf variable specified '
                             'along this dimension, and the value is the '
                             'maximum threshold for the variable\'s values '
                             'that will be saved to the output (indices that '
                             'correspond to the values that are greater than '
                             'the threshold will be excluded from the '
                             'dimension); this argument can be specified '
                             'multiple times to exclude indices from multiple '
                             'dimensions'
                             % string_int_parser.separator,
                        type=string_int_parser,
                        action=AddToDictAction,
                        dest='max_dict')
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
                        action=AddToDictAction,
                        dest='min_dict')
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
                        action=AddToDictAction,
                        dest='max_dict')


def cmd(args):
    if len(args.var_names) == 0:
        raise Exception('Variable name list \'--var-names\' is empty.')

    in_ds = Dataset(args.input_file, 'r')
    create_dir_for_file(args.output_file)
    out_ds = Dataset(args.output_file, 'w')

    processed_requests = {}
    for var_name in args.var_names:
        in_var = in_ds.variables[var_name]
        var_request = []
        for dim_name in in_var.dimensions:
            if dim_name in processed_requests:
                dim_request = processed_requests[dim_name]
            else:
                dim = in_ds.dimensions[dim_name]

                in_dim_var = None
                if dim_name in in_ds.variables:
                    in_dim_var = in_ds.variables[dim_name]
                    if len(in_dim_var.dimensions) != 1 \
                            or in_dim_var.dimensions[0] != dim_name:
                        in_dim_var = None

                dim_request = _process_slice_request(dim_name,
                                                     args.slice_dict,
                                                     slice(None))
                dim_request = _process_exclude_request(dim_name,
                                                       args.exclude_dict,
                                                       dim.size,
                                                       dim_request)
                dim_request = _process_min_max_request(dim_name,
                                                       args.min_dict,
                                                       args.max_dict,
                                                       in_dim_var,
                                                       dim_request)

                processed_requests[dim_name] = dim_request
                out_ds.createDimension(dim_name,
                                       None if dim.isunlimited()
                                       else _calc_request_size(dim_request,
                                                               dim.size))

                if in_dim_var is not None:
                    out_dim_var = out_ds.createVariable(dim_name,
                                                        in_dim_var.dtype,
                                                        dimensions=(dim_name,))
                    copy_nc_attributes(in_dim_var, out_dim_var)
                    out_dim_var[:] = in_dim_var[dim_request]

            var_request.append(dim_request)

        out_var = out_ds.createVariable(var_name,
                                        in_var.dtype,
                                        dimensions=in_var.dimensions)
        copy_nc_attributes(in_var, out_var)

        iter_mask = np.ones((len(in_var.shape, )), dtype=bool)
        iter_mask[-MAX_COPY_DIM_COUNT:] = False

        read_iter = DimIterator(in_var.shape, var_request, iter_mask)
        write_iter = DimIterator(out_var.shape, None, iter_mask)

        for read_slc, write_slc in izip(read_iter.slice_tuples(),
                                        write_iter.slice_tuples()):
            out_var[write_slc] = in_var[read_slc]

    add_history(out_ds, get_history(in_ds))
    in_ds.close()
    out_ds.close()


def _process_slice_request(dim_name, request_dict, default_request):
    if request_dict is None or dim_name not in request_dict:
        return default_request
    return request_dict[dim_name]


def _process_exclude_request(dim_name, request_dict, dim_size,
                             current_request):
    if request_dict is None or dim_name not in request_dict:
        return current_request

    request_list = request_dict[dim_name]
    if len(request_list) == 0:
        return current_request

    exclude_indices = []
    max_pos_idx = dim_size - 1
    min_neg_idx = -dim_size - 1
    for idx in request_list:
        if 0 > idx >= min_neg_idx:
            exclude_indices.append(dim_size + idx)
        elif 0 <= idx <= max_pos_idx:
            exclude_indices.append(idx)

    exclude_indices = set(exclude_indices)

    if current_request is None:
        index_list = range(dim_size)
    elif isinstance(current_request, slice):
        index_list = range(dim_size)[current_request]
    else:
        index_list = current_request

    result = [idx for idx in index_list if idx not in exclude_indices]

    if len(result) == len(index_list):
        return current_request

    return result


def _process_min_max_request(dim_name, min_dict, max_dict, dim_var,
                             current_request):
    min_val = None
    if min_dict is not None and dim_name in min_dict:
        min_val = min_dict[dim_name]

    max_val = None
    if max_dict is not None and dim_name in max_dict:
        max_val = max_dict[dim_name]

    if min_val is None and max_val is None:
        return current_request

    if dim_var is None:
        raise Exception('Input file does not contain 1D variable \'%s\' '
                        'specified along the dimension with the same name.'
                        % dim_name)

    time_converter = None
    if isinstance(min_val, datetime):
        time_converter = get_time_converter(dim_var)
        min_val = time_converter.date2num(min_val)

    if isinstance(max_val, datetime):
        if time_converter is None:
            time_converter = get_time_converter(dim_var)
        max_val = time_converter.date2num(max_val)

    if current_request is None:
        index_list = range(dim_var.size)
        var_list = dim_var[:]
    elif isinstance(current_request, slice):
        index_list = range(dim_var.size)[current_request]
        var_list = dim_var[current_request]
    else:
        index_list = current_request
        var_list = dim_var[current_request]

    if min_val is None:
        result = [dim_idx for list_idx, dim_idx in enumerate(index_list) if
                  var_list[list_idx] <= max_val]
    elif max_val is None:
        result = [dim_idx for list_idx, dim_idx in enumerate(index_list) if
                  var_list[list_idx] >= min_val]
    else:
        result = [dim_idx for list_idx, dim_idx in enumerate(index_list) if
                  max_val >= var_list[list_idx] >= min_val]

    if len(result) == len(index_list):
        return current_request

    return result


def _calc_request_size(request, dim_size):
    if request is None:
        return dim_size
    elif isinstance(request, slice):
        return len(xrange(*request.indices(dim_size)))
    else:
        return len(request)


class AddToDictAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        d = getattr(namespace, self.dest)
        if d is None:
            d = {}
            setattr(namespace, self.dest, d)
        key = values[0]
        if key in d:
            raise Exception('Only one operation of the same type per variable '
                            'is allowed.')
        d[key] = values[1]


class AddToSetAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        s = getattr(namespace, self.dest)
        if s is None:
            s = set(values)
        else:
            s = s.union(values)
        setattr(namespace, self.dest, s)
