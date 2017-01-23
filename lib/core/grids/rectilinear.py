import numpy as np

from core.grids.axes import RegularAxis, RectilinearAxis
from core.grids.grid import Grid
from core.grids.structured import StructuredGrid


class RectilinearGrid(StructuredGrid):
    def __init__(self, x, y):
        if len(x.shape) != 1 or len(y.shape) != 1:
            raise Exception('Rectilinear grid requires 1D axis arrays.')

        self._x_axis = RectilinearGrid._build_axis(x)
        self._y_axis = RectilinearGrid._build_axis(y)
        xx, yy = np.meshgrid(self._x_axis, self._y_axis)
        super(RectilinearGrid, self).__init__(xx, yy)

    @property
    def x_axis(self):
        return self._x_axis

    @property
    def y_axis(self):
        return self._y_axis

    @staticmethod
    def _build_axis(axis_values):
        count = len(axis_values)
        first = axis_values[0]
        step = axis_values[1] - axis_values[0]

        result = RegularAxis(first, count, step)

        try:
            eps = np.finfo(step).eps
        except:
            eps = 0

        if np.allclose(axis_values, result, atol=eps):
            return result
        else:
            return RectilinearAxis(axis_values)

    def calc_weights(self, x, y):
        col_indices, col_weights = self._x_axis.calc_weights(x)

        if col_indices is None:
            return Grid._EMPTY_WEIGHT_RESULT

        row_indices, row_weights = self._y_axis.calc_weights(y)

        if row_indices is None:
            return Grid._EMPTY_WEIGHT_RESULT

        i_indices = [row_indices[0], row_indices[1],
                     row_indices[1], row_indices[0]]

        j_indices = [col_indices[0], col_indices[0],
                     col_indices[1], col_indices[1]]

        weights = [row_weights[0] * col_weights[0],
                   row_weights[1] * col_weights[0],
                   row_weights[1] * col_weights[1],
                   row_weights[0] * col_weights[1]]

        return [i_indices, j_indices], weights