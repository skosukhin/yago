import numpy as np

from core.grids.cell_tree import CellTree
from core.grids.grid import Grid


class StructuredGrid(Grid):
    _MAX_LOC_TREE_DEPTH = 20
    _MAX_LOC_TREE_LEAF_SIZE = 500
    _QUAD_VERT_WEIGHTS = np.array(
        [[1.0, 0.0, 0.0, 0.0],
         [1.0, 1.0, 0.0, 0.0],
         [1.0, 1.0, 1.0, 1.0],
         [1.0, 0.0, 1.0, 0.0]])

    cell_vert_count = 4

    def __init__(self, xx, yy):
        if xx.shape != yy.shape or len(xx.shape) != 2:
            raise Exception()

        self._data = np.stack((xx, yy), axis=-1)

        self._loc_tree = None

    @property
    def shape(self):
        return self._data.shape[:-1]

    def __getitem__(self, item):
        return self._data[item]

    def init_cell_locator(self, no_gap_along_axis=None):
        self._loc_tree = CellTree(self._data, no_gap_along_axis,
                                  StructuredGrid._MAX_LOC_TREE_LEAF_SIZE,
                                  StructuredGrid._MAX_LOC_TREE_DEPTH)

    def calc_weights(self, x, y):
        cell = self._find_cell(x, y)
        if cell[0] is None:
            return Grid._EMPTY_WEIGHT_RESULT
        poly = cell[1]
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

        return (cell[0],
                [(1.0 - l) * (1.0 - m), l * (1.0 - m), l * m, (1.0 - l) * m])

    def _find_cell(self, x, y):
        if self._loc_tree is None:
            raise Exception('Cell locator is not initialized.')
        return self._loc_tree.find_cell([x, y])
