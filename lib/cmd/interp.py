from cmd.common.arg_processors import ListParser

description = 'interpolates fields specified on a structured grid in ' \
              'Cartesian coordinates'


def setup_parser(parser):
    mandatory_args = parser.add_argument_group('mandatory arguments')
    mandatory_args.add_argument('--input-file',
                                help='name of netcdf file that contains '
                                     'fields that need to be interpolated',
                                required=True)
    mandatory_args.add_argument('--grid-file',
                                help='name of netcdf file that contains '
                                     'coordinates of the target grid',
                                required=True)
    mandatory_args.add_argument('--weight-file',
                                help='name of netcdf file that interpolation '
                                     'weights will be saved to (if the file '
                                     'does not exist) or read from (if the '
                                     'file exists)')

    parser.add_argument('--output-file',
                        help='name of netcdf file that results of the '
                             'interpolation will be saved to; mandatory if '
                             'the weight file does not exist')

    string_list_parser = ListParser()
    parser.add_argument('--var-names',
                        help='\'%s\'-separated list of names of netcdf '
                             'variables that need to be interpolated; ignored '
                             'if output file is not specified'
                             % string_list_parser.separator,
                        type=string_list_parser)

    parser.add_argument('--cycled-dim',
                        help='name of netcdf dimension in the input file that '
                             'needs to be treated as if it is logically '
                             'cycled (e.g. dimension that corresponds to '
                             'longitudes of points of a global rectilinear '
                             'grid)')


def cmd(args):
    pass
