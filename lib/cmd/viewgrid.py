from netCDF4 import Dataset


import cmd.common.name_constants as names

description = 'shows generated a generated grid'


def setup_parser(parser):
    parser.add_argument('--input-file',
                        help='path to a file generated with \'gengrid\' '
                             'command',
                        required=True)


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

    in_ds.close()

    x, y = m(lons, lats)
    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='#cc9966', lake_color='#99ffff')
    m.scatter(x, y, 3, marker='o', color='k')
    plt.show()
