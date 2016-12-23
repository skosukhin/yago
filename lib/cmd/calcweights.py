from netCDF4 import Dataset

description = 'calculates weights for the following interpolation'


def setup_parser(parser):
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--output-file', required=True)
    parser.add_argument('--grid-file', required=True)


def cmd(args):
    input_ds = Dataset(args.input_file, 'r')
