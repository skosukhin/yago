from itertools import izip

import numpy as np

from core.grids.common import KDTree, BoundingBox
from core.grids.grid import Grid


class StructuredGrid(Grid):
    _MAX_LOC_TREE_DEPTH = 100
    _MAX_LOC_TREE_LEAF_LENGTH = 100
    _QUAD_VERT_WEIGHTS = np.array(
        [[1.0, 0.0, 0.0, 0.0],
         [1.0, 1.0, 0.0, 0.0],
         [1.0, 1.0, 1.0, 1.0],
         [1.0, 0.0, 1.0, 0.0]])

    cell_vert_count = 4

    def __init__(self, xx, yy):
        if xx.shape != yy.shape or len(xx.shape) != 2:
            raise Exception()

        self._data = np.dstack((xx, yy))

        self._loc_tree = None

    @property
    def shape(self):
        return self._data.shape[:-1]

    def __getitem__(self, item):
        return self._data[item]

    def calc_weights(self, x, y):
        cell = self._find_cell(x, y)
        if cell is None:
            return Grid._EMPTY_WEIGHT_RESULT
        poly = np.array(cell.poly[:-1])
        alphas = np.linalg.solve(StructuredGrid._QUAD_VERT_WEIGHTS, poly[:, 0])
        betas = np.linalg.solve(StructuredGrid._QUAD_VERT_WEIGHTS, poly[:, 1])
        aa = alphas[3] * betas[2] - alphas[2] * betas[3]
        bb = (alphas[3] * betas[0] - alphas[0] * betas[3]
              + alphas[1] * betas[2] - alphas[2] * betas[1]
              + x * betas[3] - y * alphas[3])
        cc = (alphas[1] * betas[0] - alphas[0] * betas[1]
              + x * betas[1] - y * alphas[1])

        roots = np.roots([aa, bb, cc])

        m = roots[np.logical_and(roots >= 0.0, roots <= 1.0)][0]
        l = (y - betas[0] - betas[2] * m) / (betas[1] + betas[3] * m)

        return ([cell.i_indices,
                cell.j_indices],
                [(1.0 - l) * (1.0 - m), l * (1.0 - m), l * m, (1.0 - l) * m])

    def _build_loc_tree(self):
        cell_list = []

        for i in xrange(self.shape[0] - 1):
            for j in xrange(self.shape[1] - 1):
                i_indices = [i, i + 1, i + 1, i]
                j_indices = [j, j, j + 1, j + 1]
                cell = self._data[i_indices, j_indices]
                poly = np.append(cell, cell[0, None], axis=0).tolist()
                cell_list.append(CellBoundingBox(i_indices, j_indices, poly,
                                                 cell.min(axis=0),
                                                 cell.max(axis=0)))

        return KDTree.build(cell_list,
                            StructuredGrid._MAX_LOC_TREE_LEAF_LENGTH,
                            StructuredGrid._MAX_LOC_TREE_DEPTH)

    def _find_cell(self, x, y):
        if self._loc_tree is None:
            self._loc_tree = self._build_loc_tree()

        for cell in self._loc_tree.get_boxes([x, y]):
            if self._point_inside_cell(x, y, cell.poly):
                return cell

        return None

    @staticmethod
    def _point_inside_cell(x, y, cell_poly):
        result = False
        for v1, v2 in izip(cell_poly, cell_poly[1:]):
            if (v1[1] <= y < v2[1]) or (v1[1] > y >= v2[1]):
                vt = (y - v1[1]) / (v2[1] - v1[1])
                if x < v1[0] + vt * (v2[0] - v1[0]):
                    result = not result
        return result


class CellBoundingBox(BoundingBox):
    axis_count = 2

    def __init__(self, i_indices, j_indices, poly, min_arr, max_arr):
        self.i_indices, self.j_indices = i_indices, j_indices
        self.poly = poly
        self.min_arr = min_arr
        self.max_arr = max_arr

    def get_bounds(self, axis):
        return self.min_arr[axis], self.max_arr[axis]
