class BoundingBox(object):
    def get_bounds(self, axis):
        raise NotImplementedError()

    @property
    def axis_count(self):
        raise NotImplementedError()


class KDTree(object):
    def __init__(self, left, right, split_val, axis, leaf_list):
        self.left = left
        self.right = right
        self.split_val = split_val
        self.axis = axis
        self.leaf_list = leaf_list

    def get_boxes(self, v):
        if self.leaf_list is None:
            if v[self.axis] < self.split_val:
                return self.left.get_boxes(v)
            else:
                return self.right.get_boxes(v)
        else:
            return self.leaf_list

    @staticmethod
    def build(box_list, max_leaf_list_len=None, max_depth=None):
        if max_leaf_list_len is None:
            max_leaf_list_len = 1

        if max_depth is None:
            max_depth = 0

        if max_leaf_list_len < 1 or max_depth < 0:
            raise Exception()

        return KDTree._build(box_list, max_leaf_list_len, 0, max_depth, None)

    @staticmethod
    def _build(box_list, max_leaf_list_len, depth, max_depth, split_value):
        if len(box_list) <= max_leaf_list_len or depth >= max_depth:
            return KDTree(None, None, None, None, box_list)

        split_axis = depth % box_list[0].axis_count

        if split_value is None:
            mi, ma = box_list[0].get_bounds(split_axis)
            for box in box_list[1:]:
                box_mi, box_ma = box.get_bounds(split_axis)
                if box_mi < mi:
                    mi = box_mi
                if box_ma > ma:
                    ma = box_ma

            split_value = (ma + mi) / 2.0

        next_depth = depth + 1
        next_split_axis = next_depth % box_list[0].axis_count
        left_list, right_list = [], []
        left_min, left_max, right_min, right_max = None, None, None, None

        for box in box_list:
            box_mi, box_ma = box.get_bounds(split_axis)
            box_mi_next, box_ma_next = box.get_bounds(next_split_axis)
            if box_mi < split_value:
                left_list.append(box)
                if left_min is None:
                    left_min, left_max = box_mi_next, box_ma_next
                else:
                    if box_mi_next < left_min:
                        left_min = box_mi_next
                    if box_ma_next > left_max:
                        left_max = box_ma_next

            if box_ma > split_value:
                right_list.append(box)
                if right_min is None:
                    right_min, right_max = box_mi_next, box_ma_next
                else:
                    if box_mi_next < right_min:
                        right_min = box_mi_next
                    if box_ma_next > right_max:
                        right_max = box_ma_next

        left_split_val = \
            None if left_min is None else (left_min + left_max) / 2.0
        right_split_val = \
            None if right_min is None else (right_min + right_max) / 2.0

        return KDTree(
            KDTree._build(left_list, max_leaf_list_len, next_depth, max_depth,
                          left_split_val),
            KDTree._build(right_list, max_leaf_list_len, next_depth, max_depth,
                          right_split_val),
            split_value,
            split_axis,
            None)
