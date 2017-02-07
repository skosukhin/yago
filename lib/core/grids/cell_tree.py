from itertools import izip

import numpy as np


class CellTree(object):
    def __init__(self, coords, max_leaf_size, max_depth):
        if coords is None or max_leaf_size < 1 or max_depth < 0:
            raise Exception()

        if coords.shape[-1] != len(coords.shape[:-1]):
            raise Exception()

        if coords.shape[-1] != 2:
            raise NotImplementedError('Only 2D grids are supported.')

        self.coords = coords
        self.root = CellTree._build(self.coords, max_leaf_size, max_depth)

    def find_cell(self, v):
        return self.root.find_cell(self.coords, v)

    @staticmethod
    def _build(coords, max_leaf_size, max_depth):
        vertex_slices = [[slice(0, -1) if dim_idx == 0 else slice(1, None)
                          for dim_idx in idx]
                         for idx in np.ndindex((2,) * (len(coords.shape) - 1))]
        maximums = np.maximum.reduce(
            [coords[slc] for slc in vertex_slices])
        indices = np.rollaxis(
            np.indices(maximums.shape[:-1]), 0,
            len(maximums.shape)).reshape((-1, maximums.shape[-1]))
        maximums = maximums.reshape((-1, coords.shape[-1]))
        minimums = np.minimum.reduce(
            [coords[slc] for slc in vertex_slices]).reshape(
            (-1, coords.shape[-1]))
        centers = (maximums + minimums) / 2.0

        return _CellNode.build(indices, maximums, minimums, centers,
                               max_leaf_size, 0, max_depth)


class _CellNode(object):
    _CELL_POLY_I = np.array([0, 1, 1, 0, 0])
    _CELL_POLY_J = np.array([0, 0, 1, 1, 0])

    def __init__(self, left=None, right=None, max_left=None, min_right=None,
                 axis=None, leaf=None):
        self.left = left
        self.right = right
        self.max_left = max_left
        self.min_right = min_right
        self.axis = axis
        self.leaf = leaf

    def find_cell(self, coords, v):
        if self.leaf is None:
            if v[self.axis] < self.max_left:
                result = self.left.find_cell(coords, v)
                if result[0] is not None:
                    return result
            if v[self.axis] >= self.min_right:
                return self.right.find_cell(coords, v)
        else:
            indices, minimums, maximums = self.leaf
            result_mask = np.logical_and.reduce(
                [np.logical_and(minimums[:, axis] <= v[axis],
                                maximums[:, axis] > v[axis])
                 for axis in xrange(indices.shape[-1])])

            indices = indices[result_mask]
            if indices.shape[0] > 0:
                for i, j in indices:
                    cel_ii = _CellNode._CELL_POLY_I + i
                    cel_jj = _CellNode._CELL_POLY_J + j
                    cell_poly = coords[cel_ii, cel_jj]
                    if self._point_inside_cell(v[0], v[1], cell_poly):
                        return [cel_ii[:-1], cel_jj[:-1]], cell_poly[:-1]
        return None, None

    @staticmethod
    def build(indices, maximums, minimums, centers,
              max_leaf_size, depth, max_depth):
        if indices.shape[0] <= max_leaf_size or depth >= max_depth:
            return _CellNode(leaf=(indices, minimums, maximums))

        split_axis = depth % indices.shape[-1]
        axis_centers = centers[:, split_axis]
        split_value = (np.max(axis_centers) + np.min(axis_centers)) / 2.0

        left_mask = axis_centers < split_value
        right_mask = ~left_mask
        left_indices = indices[left_mask]
        split_idx = len(left_indices)

        if split_idx == 0 or split_idx == indices.shape[0]:
            return _CellNode(leaf=(indices, minimums, maximums))

        indices[:] = np.concatenate((left_indices, indices[right_mask]))
        del left_indices

        mins_right = minimums[right_mask]
        min_right = np.min(mins_right[:, split_axis])
        minimums[:] = np.concatenate((minimums[left_mask], mins_right))
        del mins_right

        maxs_left = maximums[left_mask]
        max_left = np.max(maxs_left[:, split_axis])
        maximums[:] = np.concatenate((maxs_left, maximums[right_mask]))
        del maxs_left

        centers[:] = np.concatenate((centers[left_mask], centers[right_mask]))

        return _CellNode(
            _CellNode.build(indices[:split_idx], maximums[:split_idx],
                            minimums[:split_idx], centers[:split_idx],
                            max_leaf_size, depth + 1, max_depth),
            _CellNode.build(indices[split_idx:], maximums[split_idx:],
                            minimums[split_idx:], centers[split_idx:],
                            max_leaf_size, depth + 1, max_depth),
            max_left,
            min_right,
            split_axis)

    @staticmethod
    def _point_inside_cell(x, y, cell_poly):
        result = False
        for v1, v2 in izip(cell_poly, cell_poly[1:]):
            if (v1[1] <= y < v2[1]) or (v1[1] > y >= v2[1]):
                vt = (y - v1[1]) / (v2[1] - v1[1])
                if x < v1[0] + vt * (v2[0] - v1[0]):
                    result = not result
        return result
