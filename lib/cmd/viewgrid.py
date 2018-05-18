import numpy as np
from netCDF4 import Dataset

import cmd.common.name_constants as names
from cmd.common.nc_utils import init_converter_from_proj_var

_cmd_disabled = False
_disabled_reason_msg = None
try:
    # import matplotlib
    # matplotlib.use('TkAgg')
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt

except ImportError as e:
    _cmd_disabled = True
    _disabled_reason_msg = e.message

description = 'shows generated grids'

if _cmd_disabled:
    description += ' (DISABLED: {0})'.format(_disabled_reason_msg)


def setup_parser(parser):
    mandatory_args = parser.add_argument_group('mandatory arguments')
    mandatory_args.add_argument('--input-file',
                                help='path to a file generated with '
                                     '\'gengrid\' command',
                                required=True)


def cmd(args):
    if _cmd_disabled:
        raise Exception(
            'The command is disabled: {0}'.format(_disabled_reason_msg))

    in_ds = Dataset(args.input_file, 'r')
    proj_var = in_ds.variables[names.VAR_PROJECTION]

    ref_lon = proj_var.longitude_of_projection_origin
    ref_lat = proj_var.latitude_of_projection_origin

    m = Basemap(projection='ortho',
                rsphere=proj_var.earth_radius,
                lon_0=ref_lon,
                lat_0=ref_lat,
                resolution='l')

    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='#cc9966', lake_color='#99ffff', alpha=0.3)
    m.drawmeridians([ref_lon], latmax=90)
    m.drawparallels([ref_lat], latmax=90)

    # Plot grid
    edge_lats, edge_lons = get_edge_values(in_ds.variables[names.DIMVAR_LAT],
                                           in_ds.variables[names.DIMVAR_LON])

    edge_map_xx, edge_map_yy = m(edge_lons, edge_lats)
    m.scatter(edge_map_xx[1:], edge_map_yy[1:], s=3, color='k',
              marker='o', label='Grid boundary points')
    m.scatter(edge_map_xx[0], edge_map_yy[0], s=10, color='r',
              marker='o', label='First grid point')

    # Plot reference point and the axes of the Cartesian grid
    converter = init_converter_from_proj_var(proj_var)
    delta_x = in_ds.variables[names.DIMVAR_X].step * 0.1
    delta_y = in_ds.variables[names.DIMVAR_Y].step * 0.1

    ref_map_xx, ref_map_yy = plot_grid_axes(m, converter,
                                            ref_lat, ref_lon,
                                            delta_x, delta_y)
    m.scatter(ref_map_xx, ref_map_yy, s=20, color='b',
              marker='o', label='Projection reference point (%.1f;%.1f)' %
                                (ref_lat, ref_lon))

    # Plot latitudes of true scale
    rot_lons, rot_lats = np.meshgrid(np.linspace(0.0, 360.0),
                                     converter.projection.true_scale_lats)
    true_lats_lats, true_lats_lons = converter.rotor.restore_points(rot_lats,
                                                                    rot_lons)
    true_lats_map_xx, true_lats_map_yy = m(true_lats_lons, true_lats_lats)
    m.scatter(true_lats_map_xx, true_lats_map_yy, s=3, color='g',
              marker='o', label='Latitude(s) of true scale')

    in_ds.close()

    plt.title(args.input_file)
    plt.legend(scatterpoints=1, loc=(0, 0), fontsize='x-small')
    plt.show()


def plot_grid_axes(basemap, converter, lat, lon, delta_x, delta_y):
    # Basemap can't handle axis directions at the Poles, so we use our own
    # implementation of rotations.
    grid_x, grid_y = converter.convert_points(lat, lon)
    axis_x_lat, axis_x_lon = converter.restore_points(grid_x + delta_x, grid_y)
    axis_y_lat, axis_y_lon = converter.restore_points(grid_x, grid_y + delta_y)
    map_xx, map_yy = basemap([lon, axis_x_lon, axis_y_lon],
                             [lat, axis_x_lat, axis_y_lat])
    q = basemap.quiver([map_xx[0], map_xx[0]], [map_yy[0], map_yy[0]],
                   [map_xx[1] - map_xx[0], map_xx[2] - map_xx[0]],
                   [map_yy[1] - map_yy[0], map_yy[2] - map_yy[0]],
                   width=0.003,
                   color='k')
    # Return map coordinates of the reference point
    return map_xx[0], map_yy[0]


def get_edge_values(lat_var, lon_var):
    edge_shape = (2 * (lat_var.shape[0] + lat_var.shape[1] - 2),)
    edge_lats = np.empty(edge_shape, dtype=lat_var.dtype)
    edge_lons = np.empty(edge_shape, dtype=lon_var.dtype)

    start, end = 0, lat_var.shape[1]
    edge_lats[start:end] = lat_var[0, :]
    edge_lons[start:end] = lon_var[0, :]

    start, end = end, end + lat_var.shape[0] - 2
    edge_lats[start:end] = lat_var[1:-1, -1]
    edge_lons[start:end] = lon_var[1:-1, -1]

    start, end = end, end + lat_var.shape[1]
    edge_lats[start:end] = lat_var[-1, -1::-1]
    edge_lons[start:end] = lon_var[-1, -1::-1]

    start, end = end, end + lat_var.shape[0] - 2
    edge_lats[start:end] = lat_var[-1:1:-1, 0]
    edge_lons[start:end] = lon_var[-1:1:-1, 0]

    return edge_lats, edge_lons

