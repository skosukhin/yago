import numpy as np

from core.grids.grid import Grid


class RectilinearAxis(Grid):
    cell_vert_count = 2

    def __init__(self, values):

        if len(values.shape) != 1:
            raise Exception('Axis requires 1D array.')

        if len(values) < 2:
            raise Exception('Axis requires at least two values.')

        self._ascending = RectilinearAxis._check_ascending(values)

        if self._ascending is None:
            raise Exception('Axis values must be sorted and unique.')

        self._orig_order = np.copy(values)
        self._sorted_order = self._orig_order \
            if self._ascending \
            else self._orig_order[::-1]

        self._values = np.copy(values)

    @property
    def shape(self):
        return (len(self._orig_order),)

    def __len__(self):
        return len(self._orig_order)

    def __getitem__(self, item):
        return self._orig_order[item]

    def __str__(self):
        return self._values.__str__()

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
        super(RegularAxis, self).__init__(values)

        self.count = count
        self.step = step
        self.first = values[0]

    def calc_weights(self, x):
        if x < self._sorted_order[0] or x >= self._sorted_order[-1]:
            return Grid._EMPTY_WEIGHT_RESULT

        if self.step < 0:
            l_idx = np.intp(np.ceil((x - self.first) / self.step)) - 1
        else:
            l_idx = np.intp(np.floor((x - self.first) / self.step))

        w_r = (x - self[l_idx]) / self.step
        return [l_idx, l_idx + 1], [1 - w_r, w_r]
