class Grid(object):
    _EMPTY_WEIGHT_RESULT = None, None

    shape = None
    cell_vert_count = None

    def init_cell_locator(self, no_gap_along_axis=None):
        raise NotImplementedError()

    def calc_weights(self, *args):
        raise NotImplementedError()
