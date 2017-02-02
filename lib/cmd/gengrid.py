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

import cmd.common.name_constants as names
from cmd.common.misc import create_dir_for_file
from cmd.common.nc_utils import set_generic_lat_attributes, \
    set_generic_lon_attributes, add_or_append_history
from cmd.common.arg_processors import ListParser, parse_pos_intp, \
    parse_pos_float, init_converter_from_args
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
                                     'grid point',
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
                                     'grid point',
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
                      'projection are used (default: %%(default)s)'\
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
    parser.add_argument('--x-offset',
                        help='value added to all x-coordinates '
                             '(false easting)',
                        type=np.float64, default=np.float64(0))
    parser.add_argument('--y-offset',
                        help='value added to all y-coordinates '
                             '(false northing)',
                        type=np.float64, default=np.float64(0))
    parser.add_argument('--output-format',
                        help='output file format (default: %(default)s)',
                        choices=['txt', 'nc'], default='txt')


def cmd(args):
    converter = init_converter_from_args(args)
    xx, yy = _generate_cartesian_grid(args.x_start, args.x_count, args.x_step,
                                      args.y_start, args.y_count, args.y_step)

    lats, lons = converter.restore_points(xx, yy)

    serializer = get_serializer_from_args(args)
    serializer.title = 'Geographic coordinates of points of a regular grid ' \
                       'defined in Cartesian coordinates on a ' + \
                       converter.projection.long_name + ' projection plane.'
    serializer.xx = xx
    serializer.x_step = args.x_step
    serializer.yy = yy
    serializer.y_step = args.y_step
    serializer.proj_description = \
        re.sub(r'\s{2,}', ' ',
               converter.projection.__doc__.replace('\n', ' ')).strip() + \
        ' ' \
        'To build a projection with a given point (origin_lat;origin_lon) ' \
        'in its origin, series of rotations of the spherical coordinate ' \
        'system are applied beforehand to move the origin point of the grid ' \
        'on the position of the North Pole. The rotations are performed ' \
        'around axes of a 3D right handed Cartesian coordinate system. The ' \
        'origin of the Cartesian system is in the center of the sphere that ' \
        'approximates the Earth. X-axis points to the intersection of the ' \
        'equator and the Greenwich Meridian; Y-axis points to the ' \
        'intersection of the equator and the 90th eastern meridian; Z-axis ' \
        'points to the North Pole.'
    serializer.mapping_name = (converter.projection.standard_name +
                               '+rotated_latitude_longitude')
    serializer.earth_radius = args.earth_radius
    serializer.orig_lat = args.orig_lat
    serializer.orig_lon = args.orig_lon
    serializer.standard_parallels = converter.projection.true_scale_lats
    serializer.rot_axes_ids = converter.rotor.rot_axes_ids
    serializer.rot_angles_deg = converter.rotor.rot_angles_deg
    serializer.proj_short_name = converter.projection.short_name
    serializer.lats = lats
    serializer.lons = lons
    serializer.x_offset = args.x_offset
    serializer.y_offset = args.y_offset

    serializer.save()


def get_serializer_from_args(args):
    format_lower = args.output_format.lower()

    if format_lower == 'txt':
        return TextSerializer(args.output_file)
    elif format_lower == 'nc':
        return NetCDFSerializer(args.output_file)
    else:
        raise Exception(
            'Unknown output format: \'' + args.output_format + '\'.')


def _generate_cartesian_grid(x_start, x_count, x_step,
                             y_start, y_count, y_step):
    x_array = np.linspace(x_start, x_start + (x_count - 1) * x_step,
                          num=x_count)
    y_array = np.linspace(y_start, y_start + (y_count - 1) * y_step,
                          num=y_count)

    return np.meshgrid(x_array, y_array)


class OutputSerializer(object):
    title = None
    xx = None
    x_step = None
    x_offset = None
    yy = None
    y_step = None
    y_offset = None
    proj_description = None
    mapping_name = None
    earth_radius = None
    orig_lat = None
    orig_lon = None
    standard_parallels = None
    rot_axes_ids = None
    rot_angles_deg = None
    proj_short_name = None
    lats = None
    lons = None

    def save(self):
        pass


class NetCDFSerializer(OutputSerializer):
    def __init__(self, filename):
        self.filename = filename

    def save(self):
        create_dir_for_file(self.filename)
        from netCDF4 import Dataset
        ds = Dataset(self.filename, mode='w', format='NETCDF4')
        ds.title = self.title

        ds.createDimension(names.DIMVAR_X, len(self.xx[0]))
        x_var = ds.createVariable(names.DIMVAR_X, self.xx.dtype,
                                  dimensions=(names.DIMVAR_X,))
        x_var.long_name = 'x coordinate of projection'
        x_var.standard_name = 'projection_x_coordinate'
        x_var.axis = 'X'
        x_var.units = 'm'
        x_var.step = self.x_step
        x_var[:] = self.xx[0] + self.x_offset

        ds.createDimension(names.DIMVAR_Y, len(self.yy))
        y_var = ds.createVariable(names.DIMVAR_Y, self.yy.dtype,
                                  dimensions=(names.DIMVAR_Y,))
        y_var.long_name = 'y coordinate of projection'
        y_var.standard_name = 'projection_y_coordinate'
        y_var.axis = 'Y'
        y_var.units = 'm'
        y_var.step = self.y_step
        y_var[:] = self.yy[:, 1] + self.y_offset

        proj_var = ds.createVariable(names.VAR_PROJECTION, 'c')
        proj_var.description = self.proj_description
        proj_var.grid_mapping_name = self.mapping_name
        proj_var.earth_radius = self.earth_radius
        proj_var.latitude_of_projection_origin = self.orig_lat
        proj_var.longitude_of_projection_origin = self.orig_lon
        proj_var.standard_parallel = self.standard_parallels
        proj_var.rot_axes = self.rot_axes_ids
        proj_var.rot_angles_deg = self.rot_angles_deg
        proj_var.short_name = self.proj_short_name
        proj_var.false_easting = self.x_offset
        proj_var.false_northing = self.x_offset

        lats_var = ds.createVariable(names.DIMVAR_LAT, self.lats.dtype,
                                     dimensions=(
                                         names.DIMVAR_Y,
                                         names.DIMVAR_X))
        set_generic_lat_attributes(lats_var)
        lats_var[:, :] = self.lats

        lons_var = ds.createVariable(names.DIMVAR_LON, self.lons.dtype,
                                     dimensions=(
                                         names.DIMVAR_Y,
                                         names.DIMVAR_X))
        set_generic_lon_attributes(lons_var)
        lons_var[:, :] = self.lons

        add_or_append_history(ds)

        ds.close()


class TextSerializer(OutputSerializer):
    def __init__(self, filename):
        self.filename = filename

    _COMMENT_PREFIX = '#'
    _MAX_WRAP_WIDTH = 79

    def save(self):
        create_dir_for_file(self.filename)
        with open(self.filename, 'w') as f:
            f.writelines(TextSerializer._comment('Title:'))
            f.writelines(
                TextSerializer._wrap_and_comment(self.title, indent=4))
            f.writelines(TextSerializer._comment('Projection:'))
            f.writelines(
                TextSerializer._comment('Name: ' + self.mapping_name,
                                        indent=4))
            f.writelines(TextSerializer._comment(
                'Earth radius: %s m' %
                TextSerializer._float_to_string(self.earth_radius), indent=4))

            f.writelines(TextSerializer._comment(
                'Origin: (%s; %s)' % (
                    TextSerializer._float_to_string(self.orig_lat),
                    TextSerializer._float_to_string(self.orig_lon)),
                indent=4))

            f.writelines(TextSerializer._comment(
                'Standard parallels: ' + '; '.join(
                    str(p) for p in self.standard_parallels), indent=4))

            f.writelines(TextSerializer._comment('Rotations: ' + '; '.join(
                '%s (%s)' % (i, TextSerializer._float_to_string(a)) for i, a in
                zip(self.rot_axes_ids, self.rot_angles_deg)), indent=4))

            f.writelines(TextSerializer._comment('Description:', indent=4))

            f.writelines(
                TextSerializer._wrap_and_comment(self.proj_description,
                                                 indent=8))

            f.writelines(TextSerializer._comment(
                'Grid size: %ix%i' % (len(self.xx), len(self.xx[0]))))

            f.writelines(TextSerializer._comment(
                'X offset: %s' % TextSerializer._float_to_string(
                    self.x_offset)))
            f.writelines(TextSerializer._comment(
                'Y offset: %s' % TextSerializer._float_to_string(
                    self.y_offset)))

            f.writelines(TextSerializer._comment(
                'X coordinates of projection (step %s m)' %
                TextSerializer._float_to_string(self.x_step)))
            TextSerializer._write_matrix(self.xx + self.x_offset, f)

            f.writelines(TextSerializer._comment(
                'Y coordinates of projection (step %s m)' %
                TextSerializer._float_to_string(self.x_step)))
            TextSerializer._write_matrix(self.yy + self.y_offset, f)

            f.writelines(TextSerializer._comment(
                'Latitude coordinates (degrees north)' % self.x_step))
            TextSerializer._write_matrix(self.lats, f)

            f.writelines(TextSerializer._comment(
                'Longitude coordinates (degrees east)' % self.x_step))
            TextSerializer._write_matrix(self.lons, f)

    @staticmethod
    def _wrap_and_comment(text, **kwargs):
        indent = kwargs.get('indent', 1)
        import textwrap
        return TextSerializer._comment(
            *textwrap.wrap(text, width=(
                TextSerializer._MAX_WRAP_WIDTH - indent - len(
                    TextSerializer._COMMENT_PREFIX))), indent=indent)

    @staticmethod
    def _comment(*lines, **kwargs):
        indent = kwargs.get('indent', 1)
        return [TextSerializer._COMMENT_PREFIX + (' ' * indent) + l + '\n' for
                l in lines]

    @staticmethod
    def _float_to_string(value):
        return ('%.15f' % value).rstrip('0').rstrip('.')

    @staticmethod
    def _write_matrix(matrix, stream):
        for row in matrix:
            stream.write(
                '\t'.join(
                    TextSerializer._float_to_string(v) for v in row) + '\n')
