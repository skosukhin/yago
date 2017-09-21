"""
Example parameters:
Common:
--x-start=-2695000.0
--x-count=386
--x-step=14000.0
--y-start=-2331000.0
--y-count=334
--y-step=14000.0
--orig-lat=88.9899731326
--orig-lon=-129.805571092
--adjust-angle=-39.805571092

For Polar stereographic:
--proj-name=stereo
--true-scale-lats=71.6577131288

For Mercator:
--proj-name=mercator
--true-scale-lats=10.6352550282

For Lambert:
--proj-name=lambert
--true-scale-lats=33.9172241958,54.4707286812

Example commands:
gengrid --x-start=-2695000.0 --x-count=386 --x-step=14000.0 --y-start=-2331000.0 --y-count=334 --y-step=14000.0 --orig-lat=88.9899731326 --orig-lon=-129.805571092 --adjust-angle=-39.805571092 --proj-name=stereo --true-scale-lats=71.6577131288 --output-file=grid_stereo.nc --output-format=nc
gengrid --x-start=-2695000.0 --x-count=386 --x-step=14000.0 --y-start=-2331000.0 --y-count=334 --y-step=14000.0 --orig-lat=88.9899731326 --orig-lon=-129.805571092 --adjust-angle=-39.805571092 --proj-name=mercator --true-scale-lats=10.6352550282 --output-file=grid_mercator.nc  --output-format=nc
gengrid --x-start=-2695000.0 --x-count=386 --x-step=14000.0 --y-start=-2331000.0 --y-count=334 --y-step=14000.0 --orig-lat=88.9899731326 --orig-lon=-129.805571092 --adjust-angle=-39.805571092 --proj-name=lambert --true-scale-lats=33.9172241958;54.4707286812 --output-file=grid_lambert.nc  --output-format=nc
gengrid --x-start=-2695000.0 --x-count=386 --x-step=14000.0 --y-start=-2331000.0 --y-count=334 --y-step=14000.0 --orig-lat=88.9899731326 --orig-lon=-129.805571092 --adjust-angle=-39.805571092 --proj-name=sinus --output-file=grid_sinus.nc  --output-format=nc
"""

import re

import numpy as np
from netCDF4 import Dataset

import cmd.common.name_constants as names
from cmd.common.arg_processors import ListParser, parse_pos_intp, \
    parse_pos_float, init_converter_from_args
from cmd.common.misc import create_dir_for_file
from cmd.common.nc_utils import add_or_append_history
from core.grids.axes import RegularAxis
from core.grids.rectilinear import RectilinearGrid
from core.projections import projections

description = 'generates a regular grid in Cartesian coordinates on a ' \
              'specified projection plane'

_DEFAULT_EARTH_RADIUS = np.float64(6370997.0)


def setup_parser(parser):
    mandatory_args = parser.add_argument_group('mandatory arguments')
    mandatory_args.add_argument('--proj-name',
                                help='name of projection to be used for grid '
                                     'generation',
                                choices=projections.keys(), required=True)
    mandatory_args.add_argument('--orig-lat',
                                help='latitude (in degrees) of the projection '
                                     'center',
                                type=np.float64, required=True)
    mandatory_args.add_argument('--orig-lon',
                                help='longitude (in degrees) of the '
                                     'projection center',
                                type=np.float64, required=True)
    mandatory_args.add_argument('--x-start',
                                help='x-coordinate (in meters) of the first '
                                     'grid point (after false easting)',
                                type=np.float64, required=True)
    mandatory_args.add_argument('--x-count',
                                help='number of grid points along x-axis',
                                type=parse_pos_intp, required=True)
    mandatory_args.add_argument('--x-step',
                                help='distance (in meters) between two '
                                     'consecutive grid points along x-axis',
                                type=parse_pos_float, required=True)
    mandatory_args.add_argument('--y-start',
                                help='y-coordinate (in meters) of the first '
                                     'grid point (after false northing)',
                                type=np.float64, required=True)
    mandatory_args.add_argument('--y-count',
                                help='number of grid points along y-axis',
                                type=parse_pos_intp, required=True)
    mandatory_args.add_argument('--y-step',
                                help='distance (in meters) between two '
                                     'consecutive grid points along y-axis',
                                type=parse_pos_float, required=True)
    mandatory_args.add_argument('--output-file',
                                help='output filename',
                                required=True)

    float_list_parser = ListParser(np.float64)
    true_scale_help = '\'%s\'-separated list of projection\'s latitudes (in ' \
                      'degrees) of true scale; number of values in the list ' \
                      'depend on the projection to be used; if the list is ' \
                      'empty than the default values that depend on the ' \
                      'projection are used (default: %%(default)s)' \
                      % float_list_parser.separator
    parser.add_argument('--true-scale-lats',
                        help=true_scale_help,
                        type=float_list_parser, default=[])
    parser.add_argument('--earth-radius',
                        help='earth radius (in meters) to be used for '
                             'projection (default: %(default)s)',
                        type=np.float64, default=_DEFAULT_EARTH_RADIUS)
    parser.add_argument('--adjust-angle',
                        help='optional angle of rotation for grid adjustment '
                             '(default: %(default)s)',
                        type=np.float64, default=np.float64(0.0))
    parser.add_argument('--easting',
                        help='value added to all x-coordinates '
                             '(false easting)',
                        type=np.float64, default=np.float64(0))
    parser.add_argument('--northing',
                        help='value added to all y-coordinates '
                             '(false northing)',
                        type=np.float64, default=np.float64(0))


def cmd(args):
    converter = init_converter_from_args(args)

    grid = RectilinearGrid(
        RegularAxis(args.x_start, args.x_count, args.x_step),
        RegularAxis(args.y_start, args.y_count, args.y_step),
        False)

    create_dir_for_file(args.output_file)
    grid_ds = Dataset(args.output_file, mode='w', format='NETCDF4')
    grid_ds.title = 'Geographic coordinates of points of a regular grid ' \
                    'defined in Cartesian coordinates on a ' + \
                    converter.projection.long_name + ' projection plane.'

    grid_ds.createDimension(names.DIMVAR_X, args.x_count)
    x_var = grid_ds.createVariable(names.DIMVAR_X, grid.dtype,
                                   dimensions=(names.DIMVAR_X,))
    x_var.long_name = 'x coordinate of projection'
    x_var.standard_name = 'projection_x_coordinate'
    x_var.axis = 'X'
    x_var.units = 'm'
    x_var.step = args.x_step
    x_var[:] = grid.x_axis

    grid_ds.createDimension(names.DIMVAR_Y, args.y_count)
    y_var = grid_ds.createVariable(names.DIMVAR_Y, grid.dtype,
                                   dimensions=(names.DIMVAR_Y,))
    y_var.long_name = 'y coordinate of projection'
    y_var.standard_name = 'projection_y_coordinate'
    y_var.axis = 'Y'
    y_var.units = 'm'
    y_var.step = args.y_step
    y_var[:] = grid.y_axis

    proj_var = grid_ds.createVariable(names.VAR_PROJECTION, 'c')
    proj_var.description = \
        re.sub(r'\s{2,}', ' ',
               converter.projection.__doc__.replace('\n', ' ')).strip() + \
        ' ' \
        'Before applying the projection, a series of rotations of the ' \
        'geographical coordinate system is performed to shift the point ' \
        '(origin_lat;origin_lon) to its center and to adjust the ' \
        'orientation of the axes of the projection plane with respect to ' \
        'the surface.'
    proj_var.grid_mapping_name = (converter.projection.standard_name +
                                  '+rotated_latitude_longitude')
    proj_var.earth_radius = args.earth_radius
    proj_var.latitude_of_projection_origin = args.orig_lat
    proj_var.longitude_of_projection_origin = args.orig_lon
    proj_var.standard_parallel = converter.projection.true_scale_lats
    proj_var.rot_axes = converter.rotor.rot_axes_ids
    proj_var.rot_angles_deg = converter.rotor.rot_angles_deg
    proj_var.short_name = converter.projection.short_name
    proj_var.false_easting = converter.translator.easting
    proj_var.false_northing = converter.translator.northing

    lats, lons = converter.restore_points(grid[:, :, 0], grid[:, :, 1])

    lats_var = grid_ds.createVariable(names.DIMVAR_LAT, lats.dtype,
                                      dimensions=(
                                          names.DIMVAR_Y,
                                          names.DIMVAR_X))
    lats_var.units = 'degrees_north'
    lats_var.long_name = 'latitude coordinate'
    lats_var.standard_name = 'latitude'
    lats_var[:] = lats

    lons_var = grid_ds.createVariable(names.DIMVAR_LON, lons.dtype,
                                      dimensions=(
                                          names.DIMVAR_Y,
                                          names.DIMVAR_X))
    lons_var.units = 'degrees_east'
    lons_var.long_name = 'longitude coordinate'
    lons_var.standard_name = 'longitude'
    lons_var[:] = lons

    add_or_append_history(grid_ds)

    grid_ds.close()
