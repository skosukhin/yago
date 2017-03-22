from cmd.common.arg_processors import ListParser

description = 'extends input fields specified on a rectilinear lat/lon grid ' \
              'by adding grid points that correspond to the North (South) ' \
              'Pole; values for the new grid points are calculated as mean ' \
              'values along the highest (lowest) latitude of the fields; ' \
              'u-components of new vector values are given in the direction ' \
              'of the 90th meridian; v-components of new vector values are ' \
              'given in the direction of the 180th (for the North Pole) or ' \
              'zero (for the South Pole) meridian'


def setup_parser(parser):
    mandatory_args = parser.add_argument_group('mandatory arguments')
    mandatory_args.add_argument('--input-file',
                                help='name of netcdf file that contains data '
                                     'that need to be modified',
                                required=True)
    mandatory_args.add_argument('--output-file',
                                help='output filename',
                                required=True)
    mandatory_args.add_argument('--lat-name',
                                help='name of 1D netcdf variable that '
                                     'contains latitudes',
                                required=True)
    mandatory_args.add_argument('--lon-name',
                                help='name of 1D netcdf variable that '
                                     'contains longitudes',
                                required=True)

    string_list_parser = ListParser()
    parser.add_argument('--var-names',
                        help='\'%s\'-separated list of names of netcdf '
                             'variables that will be extended and saved '
                             'to the output file; if the list is empty than '
                             'only lat/lon coordinates will be extended; if '
                             'a couple of variables contain vector field '
                             'components aligned with meridians and parallels '
                             'than they should appear as a single entry in '
                             'the list, joined with the symbol \'+\' '
                             '(e.g. uwnd+vwnd)'
                             % string_list_parser.separator,
                        type=string_list_parser)
    parser.add_argument('--add',
                        help='the Pole to add (default: %(default)s)',
                        choices=['north', 'south', 'both'], default='north')


def cmd(args):
    pass
