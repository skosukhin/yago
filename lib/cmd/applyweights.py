import numpy as np
from netCDF4 import Dataset

from cmd.common import parse_list_of_strings

description = 'applies weights to perform interpolation'


def setup_parser(parser):
    parser.add_argument('--input-files', type=parse_list_of_strings,
                        required=True)
    parser.add_argument('--output-file', required=True)
    parser.add_argument('--weight-file', required=True)
    parser.add_argument('--data-var-names', type=parse_list_of_strings,
                        required=True)


def cmd(args):
    if len(args.data) > 1 or len(args.input_files) > 1:
        raise Exception('TODO')

    weights_ds = Dataset(args.weight_file, 'r')
