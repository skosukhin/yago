class Grid(object):
    _EMPTY_WEIGHT_RESULT = None, None

    shape = None
    cell_vert_count = None

    def calc_weights(self, *args):
        raise NotImplementedError()
