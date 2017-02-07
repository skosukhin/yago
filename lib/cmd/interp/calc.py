import numpy as np
from netCDF4 import Dataset

import cmd.common.name_constants as names
from cmd.common.misc import create_dir_for_file
from cmd.common.nc_utils import add_or_append_history, DimIterator
from core.grids.rectilinear import RectilinearGrid
from core.grids.structured import StructuredGrid

description = 'calculates weights for the following interpolation'


def setup_parser(parser):
    mandatory_args = parser.add_argument_group('mandatory arguments')
    mandatory_args.add_argument('--in-grid-file',
                                help='name of netcdf file that contains '
                                     'coordinates of the input grid',
                                required=True)
    mandatory_args.add_argument('--out-grid-file',
                                help='name of netcdf file that contains '
                                     'coordinates of the output grid',
                                required=True)
    mandatory_args.add_argument('--weight-file',
                                help='name of netcdf file that interpolation '
                                     'weights will be saved to',
                                required=True)
    parser.add_argument('--in-x-name',
                        help='name of netcdf variable that contains '
                             'x-coordinates of the input grid (default: '
                             '\'%(default)s\')',
                        default=names.DIMVAR_X)
    parser.add_argument('--in-y-name',
                        help='name of netcdf variable that contains '
                             'y-coordinates of the input grid (default: '
                             '\'%(default)s\')',
                        default=names.DIMVAR_Y)
    parser.add_argument('--out-x-name',
                        help='name of netcdf variable that contains '
                             'x-coordinates of the output grid (default: '
                             '\'%(default)s\')',
                        default=names.DIMVAR_X)
    parser.add_argument('--out-y-name',
                        help='name of netcdf variable that contains '
                             'y-coordinates of the output grid (default: '
                             '\'%(default)s\')',
                        default=names.DIMVAR_Y)


def cmd(args):
    in_grid_ds = Dataset(args.in_grid_file, 'r')
    in_grid, in_grid_dim_names = \
        _init_grid(in_grid_ds.variables[args.in_x_name],
                   in_grid_ds.variables[args.in_y_name])
    in_grid_ds.close()

    out_grid_ds = Dataset(args.out_grid_file, 'r')
    out_grid, out_grid_dim_names = \
        _init_grid(out_grid_ds.variables[args.out_x_name],
                   out_grid_ds.variables[args.out_y_name])
    out_grid_ds.close()

    if len(in_grid.shape) != len(out_grid.shape):
        raise Exception()

    indices = np.ma.masked_all(out_grid.shape +
                               (len(in_grid.shape),
                                in_grid.cell_vert_count),
                               dtype=np.intp)

    weights = np.ma.masked_all(out_grid.shape +
                               (in_grid.cell_vert_count,),
                               dtype=np.float64)

    op_iter = DimIterator(out_grid.shape)
    for slc in op_iter.slice_tuples():
        cell_indices, cell_weights = in_grid.calc_weights(*out_grid[slc])
        if cell_indices is not None:
            indices[slc] = cell_indices
            weights[slc] = cell_weights

    create_dir_for_file(args.weight_file)
    out_ds = Dataset(args.weight_file, 'w')

    for dim_idx, dim_name in enumerate(out_grid_dim_names):
        out_ds.createDimension(dim_name, out_grid.shape[dim_idx])

    out_ds.createDimension(names.DIM_DIM_IDX, len(in_grid.shape))
    out_ds.createDimension(names.DIM_VERTEX_IDX,
                           in_grid.cell_vert_count)

    out_grid_shape_var = out_ds.createVariable(names.VAR_INPUT_SHAPE,
                                               np.intp,
                                               dimensions=
                                               (names.DIM_DIM_IDX,))
    out_grid_shape_var[:] = in_grid.shape

    out_grid_dim_names_var = out_ds.createVariable(names.VAR_INPUT_DIMS,
                                                   str,
                                                   dimensions=
                                                   (names.DIM_DIM_IDX,))
    for dim_idx, dim_name in enumerate(in_grid_dim_names):
        out_grid_dim_names_var[dim_idx] = dim_name

    out_indices_var = \
        out_ds.createVariable(names.VAR_INDICES,
                              indices.dtype,
                              dimensions=
                              out_grid_dim_names +
                              (names.DIM_DIM_IDX,
                               names.DIM_VERTEX_IDX))

    out_indices_var[:] = indices

    out_weights_var = \
        out_ds.createVariable(names.VAR_WEIGHTS,
                              weights.dtype,
                              dimensions=
                              out_grid_dim_names +
                              (names.DIM_VERTEX_IDX,))

    out_weights_var[:] = weights

    add_or_append_history(out_ds)
    out_ds.close()


def _init_grid(x_var, y_var):
    if len(x_var.shape) != len(y_var.shape) or len(x_var.shape) > 2:
        raise Exception()

    if len(x_var.shape) == 2:
        if x_var.dimensions != y_var.dimensions:
            raise Exception()
        return StructuredGrid(xx=x_var[:], yy=y_var[:]), x_var.dimensions
    else:
        return RectilinearGrid(x=x_var[:], y=y_var[:]), \
               y_var.dimensions + x_var.dimensions
