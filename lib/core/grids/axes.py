import copy

import numpy as np

from core.grids.grid import Grid


class RectilinearAxis(Grid):
    cell_vert_count = 2

    def __init__(self, values, copy_values=True):

        if len(values.shape) != 1:
            raise Exception('Axis requires 1D array.')

        if len(values) < 2:
            raise Exception('Axis requires at least two values.')

        self._ascending = RectilinearAxis._check_ascending(values)

        if self._ascending is None:
            raise Exception('Axis values must be sorted and unique.')

        self._orig_order = np.copy(values) if copy_values else values
        self._sorted_order = self._orig_order \
            if self._ascending \
            else self._orig_order[::-1]

    @property
    def shape(self):
        return self._orig_order.shape

    @property
    def dtype(self):
        return self._orig_order.dtype

    def __len__(self):
        return len(self._orig_order)

    def __getitem__(self, item):
        return self._orig_order[item]

    def __str__(self):
        return self._orig_order.__str__()

    def init_cell_locator(self, no_gap_along_axis=None):
        if no_gap_along_axis is not None:
            raise NotImplementedError()

    def calc_weights(self, x):
        r_idx_sorted = np.searchsorted(self._sorted_order, x, side='right')
        if r_idx_sorted == 0 or r_idx_sorted == len(self._orig_order):
            return Grid._EMPTY_WEIGHT_RESULT
        l_idx_sorted = r_idx_sorted - 1
        l_val_sorted = self._sorted_order[l_idx_sorted]
        r_val_sorted = self._sorted_order[r_idx_sorted]
        w_l_sorted = (r_val_sorted - x) / (r_val_sorted - l_val_sorted)
        if self._ascending:
            return [l_idx_sorted, r_idx_sorted], [w_l_sorted, 1.0 - w_l_sorted]
        else:
            r_idx = len(self._orig_order) - r_idx_sorted
            return [r_idx - 1, r_idx], [1.0 - w_l_sorted, w_l_sorted]

    @staticmethod
    def _check_ascending(values):
        if all(values[idx] < values[idx + 1]
               for idx in xrange(len(values) - 1)):
            return True
        elif all(values[idx] > values[idx + 1]
                 for idx in xrange(len(values) - 1)):
            return False
        else:
            return None


class RegularAxis(RectilinearAxis):
    def __init__(self, first, count, step):
        if count < 2:
            raise Exception('Axis requires at least two values.')

        values, step = np.linspace(first, first + step * (count - 1),
                                   num=count, retstep=True)
        super(RegularAxis, self).__init__(values, False)

        self.count = count
        self.step = step
        self.first = values[0]

    def init_cell_locator(self, no_gap_along_axis=None):
        if no_gap_along_axis is not None:
            raise NotImplementedError()

    def calc_weights(self, x):
        if x < self._sorted_order[0] or x >= self._sorted_order[-1]:
            return Grid._EMPTY_WEIGHT_RESULT

        if self.step < 0:
            l_idx = np.intp(np.ceil((x - self.first) / self.step)) - 1
        else:
            l_idx = np.intp(np.floor((x - self.first) / self.step))

        w_r = (x - self[l_idx]) / self.step
        return [l_idx, l_idx + 1], [1 - w_r, w_r]


def build_axis(axis_values, copy_values=True):
    if isinstance(axis_values, RegularAxis):
        return copy.deepcopy(axis_values) if copy_values else axis_values

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
        return axis_values \
            if (not copy_values and isinstance(axis_values, RectilinearAxis)) \
            else RectilinearAxis(axis_values, copy_values)
