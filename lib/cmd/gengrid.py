"""
Recommended parameters:
Common:
--x-count=386
--x-step=14000.0
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

gengrid --x-count=386 --x-step=14000.0 --y-count=334 --y-step=14000.0 --orig-lat=88.9899731326 --orig-lon=-129.805571092 --adjust-angle=-39.805571092 --proj-name=stereo --true-scale-lats=71.6577131288 --output-file=grid.nc --output-format=nc
gengrid --x-count=386 --x-step=14000.0 --y-count=334 --y-step=14000.0 --orig-lat=88.9899731326 --orig-lon=-129.805571092 --adjust-angle=-39.805571092 --proj-name=mercator --true-scale-lats=10.6352550282 --output-file=grid.nc  --output-format=nc
gengrid --x-count=386 --x-step=14000.0 --y-count=334 --y-step=14000.0 --orig-lat=88.9899731326 --orig-lon=-129.805571092 --adjust-angle=-39.805571092 --proj-name=lambert --true-scale-lats=33.9172241958,54.4707286812 --output-file=grid.nc  --output-format=nc
"""

import os
import re

import numpy as np
from core.converter import Converter, restore_points
from core.projections.mercator import MercatorProjector
from core.projections.polar_stereographic import PolarStereographicProjector

from cmd.common import parse_list_of_floats, \
    build_rotor_for_polar_stereographic, build_rotor_for_mercator, \
    build_rotor_for_lambert, generate_cartesian_grid
from core.projections.lambert import LambertConformalProjector

description = 'generates grids'


def setup_parser(parser):
    parser.add_argument('--proj-name', required=True)
    parser.add_argument('--output-file', required=True)
    parser.add_argument('--output-format', default='txt')
    parser.add_argument('--x-count', type=np.intp, required=True)
    parser.add_argument('--x-step', type=np.float64, required=True)
    parser.add_argument('--x-offset', type=np.float64, default=np.float64(0))
    parser.add_argument('--y-count', type=np.intp, required=True)
    parser.add_argument('--y-step', type=np.float64, required=True)
    parser.add_argument('--y-offset', type=np.float64, default=np.float64(0))
    parser.add_argument('--orig-lat', type=np.float64, required=True)
    parser.add_argument('--orig-lon', type=np.float64, required=True)
    parser.add_argument('--adjust-angle', type=np.float64, required=True)
    parser.add_argument('--true-scale-lats', type=parse_list_of_floats,
                        required=True)
    parser.add_argument('--earth-radius', type=np.float64,
                        default=np.float64(6370997.0))


def cmd(args):
    converter = get_converter_from_args(args)
    xx, yy = generate_cartesian_grid(args.x_count, args.x_step, args.y_count,
                                     args.y_step)
    xx += args.x_offset
    yy += args.y_offset
    lats, lons = restore_points(xx, yy, converter)

    serializer = get_serializer_from_args(args)
    serializer.title = 'Geographic coordinates of points of a regular grid ' \
                       'defined in Cartesian coordinates on a ' + \
                       converter.projector.long_name + ' projection plane.'
    serializer.xx = xx
    serializer.x_step = args.x_step
    serializer.yy = yy
    serializer.y_step = args.y_step
    serializer.proj_description = \
        re.sub(r'\s{2,}', ' ',
               converter.projector.__doc__.replace('\n', ' ')).strip() + \
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
    serializer.mapping_name = (converter.projector.standard_name +
                               '+rotated_latitude_longitude')
    serializer.earth_radius = args.earth_radius
    serializer.orig_lat = args.orig_lat
    serializer.orig_lon = args.orig_lon
    serializer.standard_parallels = args.true_scale_lats
    serializer.rot_axes_ids = converter.rotor.rot_axes_ids
    serializer.rot_angles_deg = converter.rotor.rot_angles_deg
    serializer.proj_short_name = converter.projector.short_name
    serializer.lats = lats
    serializer.lons = lons
    serializer.x_offset = args.x_offset
    serializer.y_offset = args.y_offset

    serializer.save()


def get_converter_from_args(args):
    if args.proj_name == 'stereo':
        converter = Converter(
            build_rotor_for_polar_stereographic(args.orig_lat, args.orig_lon,
                                                args.adjust_angle),
            PolarStereographicProjector(args.true_scale_lats[0],
                                        args.earth_radius))
    elif args.proj_name == 'mercator':
        converter = Converter(
            build_rotor_for_mercator(args.orig_lat, args.orig_lon,
                                     args.adjust_angle),
            MercatorProjector(args.true_scale_lats[0], args.earth_radius))
    elif args.proj_name == 'lambert':
        converter = Converter(
            build_rotor_for_lambert(args.orig_lat, args.orig_lon,
                                    args.adjust_angle),
            LambertConformalProjector(args.true_scale_lats[0],
                                      args.true_scale_lats[1],
                                      args.earth_radius))
    else:
        raise Exception('Unknown projection: \'' + args.proj_name + '\'.')

    return converter


def get_serializer_from_args(args):
    format_lower = args.output_format.lower()

    if format_lower == 'txt' or format_lower == 'text':
        return TextSerializer(args.output_file)
    elif format_lower == 'nc' or format_lower == 'netcdf':
        return NetCDFSerializer(args.output_file)
    else:
        raise Exception(
            'Unknown output format: \'' + args.output_format + '\'.')


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

    @staticmethod
    def _create_dir_for_file(filename):
        try:
            os.makedirs(os.path.dirname(filename))
        except:
            pass


class NetCDFSerializer(OutputSerializer):
    def __init__(self, filename):
        self.filename = filename

    def save(self):
        self._create_dir_for_file(self.filename)
        from netCDF4 import Dataset
        ds = Dataset(self.filename, mode='w', format='NETCDF4')
        ds.title = self.title

        ds.createDimension('x', len(self.xx[0]))
        x_var = ds.createVariable('x', 'f', dimensions=('x',))
        x_var.long_name = 'x coordinate of projection'
        x_var.standard_name = 'projection_x_coordinate'
        x_var.axis = 'X'
        x_var.units = 'm'
        x_var.step = self.x_step
        x_var[:] = self.xx[0]

        ds.createDimension('y', len(self.xx))
        y_var = ds.createVariable('y', 'f', dimensions=('y',))
        y_var.long_name = 'y coordinate of projection'
        y_var.standard_name = 'projection_y_coordinate'
        y_var.axis = 'Y'
        y_var.units = 'm'
        y_var.step = self.y_step
        y_var[:] = self.yy[:, 1]

        proj_var = ds.createVariable('projection', 'c')
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

        lats_var = ds.createVariable('lat', 'd', dimensions=('y', 'x'))
        lats_var.units = 'degrees_north'
        lats_var.long_name = 'latitude coordinate'
        lats_var.standard_name = 'latitude'
        lats_var[:, :] = self.lats

        lons_var = ds.createVariable('lon', 'd', dimensions=('y', 'x'))
        lons_var.units = 'degrees_east'
        lons_var.long_name = 'longitude coordinate'
        lons_var.standard_name = 'longitude'
        lons_var[:, :] = self.lons

        ds.close()


class TextSerializer(OutputSerializer):
    def __init__(self, filename):
        self.filename = filename

    _COMMENT_PREFIX = '#'
    _MAX_WRAP_WIDTH = 79

    def save(self):
        self._create_dir_for_file(self.filename)
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
            TextSerializer._write_matrix(self.xx, f)

            f.writelines(TextSerializer._comment(
                'Y coordinates of projection (step %s m)' %
                TextSerializer._float_to_string(self.x_step)))
            TextSerializer._write_matrix(self.yy, f)

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
