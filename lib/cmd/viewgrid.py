import numpy as np
from netCDF4 import Dataset

import cmd.common.name_constants as names
from cmd.common.nc_utils import init_converter_from_proj_var

description = 'shows generated grids'


def setup_parser(parser):
    mandatory_args = parser.add_argument_group('mandatory arguments')
    mandatory_args.add_argument('--input-file',
                                help='path to a file generated with '
                                     '\'gengrid\' command',
                                required=True)
    parser.add_argument('--show-true-scale-lats',
                        help='show latitudes of true scale',
                        action='store_true',
                        default=False)


def cmd(args):
    try:
        # import matplotlib
        # matplotlib.use('TkAgg')
        from mpl_toolkits.basemap import Basemap
        import matplotlib.pyplot as plt
    except ImportError:
        raise Exception('Failed to load Basemap module.')

    in_ds = Dataset(args.input_file, 'r')
    lats = in_ds.variables[names.DIMVAR_LAT][:]
    lons = in_ds.variables[names.DIMVAR_LON][:]
    proj_var = in_ds.variables[names.VAR_PROJECTION]

    m = Basemap(projection='ortho',
                rsphere=proj_var.earth_radius,
                lon_0=proj_var.longitude_of_projection_origin,
                lat_0=proj_var.latitude_of_projection_origin,
                resolution='l')

    x, y = m(lons, lats)
    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='#cc9966', lake_color='#99ffff')
    m.scatter(x, y, 3, marker='o', color='k')

    if args.show_true_scale_lats:
        converter = init_converter_from_proj_var(proj_var)
        rot_lons, rot_lats = np.meshgrid(np.linspace(0.0, 360.0),
                                         converter.projection.true_scale_lats)
        true_lats_lats, true_lats_lons = \
            converter.rotor.restore_points(rot_lats, rot_lons)

        true_lats_xx, true_lats_yy = m(true_lats_lons, true_lats_lats)
        m.scatter(true_lats_xx, true_lats_yy, 3, marker='o', color='r')

    in_ds.close()
    plt.show()
