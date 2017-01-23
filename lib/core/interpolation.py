import numpy as np

from core.common import min_max


class RegularGrid(object):
    def __init__(self,
                 x_start, x_count, x_step,
                 y_start, y_count, y_step):
        self.x_start = x_start
        self.x_count = x_count
        self.x_step = x_step
        self.y_start = y_start
        self.y_count = y_count
        self.y_step = y_step
        self.x_end = x_start + (x_count - 1) * x_step
        self.y_end = y_start + (y_count - 1) * y_step

    def get_x(self, col_num):
        return self.x_start + col_num * self.x_step

    def get_y(self, row_num):
        return self.y_start + row_num * self.y_step

    def get_right_x_idx(self, x):
        return np.intp(np.ceil((x - self.x_start) / self.x_step))

    def get_top_y_idx(self, y):
        return np.intp(np.ceil((y - self.y_start) / self.y_step))


def calc_weights(quad_xx, quad_yy, xx_cycled, regular_grid,
                 progress_callback=None):
    quad_indices = -1 * np.ones(
        (regular_grid.y_count, regular_grid.x_count, 2, 4), dtype=np.intp)

    weights = np.zeros(
        (regular_grid.y_count, regular_grid.x_count, 4), dtype=np.float64)

    row_count = len(quad_xx) - 1
    for m in xrange(row_count):
        if progress_callback:
            progress_callback(m, row_count)
        m_next = m + 1
        quad_col_idx_num = len(quad_xx[m])
        if not xx_cycled:
            quad_col_idx_num -= 1
        m_indices = m, m, m_next, m_next

        for n in xrange(quad_col_idx_num):
            n_next = n + 1
            if xx_cycled:
                n_next %= quad_col_idx_num
            n_indices = n, n_next, n_next, n
            poly_xx = quad_xx[m_indices, n_indices]
            mi_x, ma_x = min_max(poly_xx)
            if mi_x <= regular_grid.x_end and ma_x >= regular_grid.x_start:
                poly_yy = quad_yy[m_indices, n_indices]
                mi_y, ma_y = min_max(poly_yy)
                if mi_y <= regular_grid.y_end and ma_y >= regular_grid.y_start:
                    if mi_x > regular_grid.x_start:
                        mi_x_idx = regular_grid.get_right_x_idx(mi_x)
                    else:
                        mi_x_idx = 0

                    if ma_x < regular_grid.x_end:
                        ma_x_idx = regular_grid.get_right_x_idx(ma_x)
                    else:
                        ma_x_idx = regular_grid.x_count

                    if mi_y > regular_grid.y_start:
                        mi_y_idx = regular_grid.get_top_y_idx(mi_y)
                    else:
                        mi_y_idx = 0

                    if ma_y < regular_grid.y_end:
                        ma_y_idx = regular_grid.get_top_y_idx(ma_y)
                    else:
                        ma_y_idx = regular_grid.y_count

                    for i in xrange(mi_y_idx, ma_y_idx):
                        y = regular_grid.get_y(i)
                        for j in xrange(mi_x_idx, ma_x_idx):
                            x = regular_grid.get_x(j)
                            if point_inside_poly(x, y, poly_xx, poly_yy):
                                if quad_indices[i, j, 0, 0] >= 0:
                                    raise Exception("Quad overlap.")
                                quad_indices[i, j, 0, :] = m_indices
                                quad_indices[i, j, 1, :] = n_indices
                                weights[i, j, :] = calc_quad_weights(
                                    poly_xx, poly_yy, x, y)

    if progress_callback:
        progress_callback(row_count, row_count)

    return quad_indices, weights


def apply_weights(input_field, quad_indices, weights, dtype=None):
    result = np.ma.masked_all(weights.shape[:-1],
                              dtype=dtype if dtype is not None
                              else weights.dtype)

    for i in xrange(weights.shape[0]):
        for j in xrange(weights.shape[1]):
            if quad_indices[i, j, 0, 0] >= 0:
                input_data = input_field[quad_indices[i, j, 0, :],
                                         quad_indices[i, j, 1, :]]
                if not np.ma.is_masked(input_data):
                    result[i, j] = np.dot(input_data, weights[i, j, :])

    return result


def point_inside_poly(p_x, p_y, v_xx, v_yy):
    result = False
    v_len = len(v_xx)
    for i in range(v_len):
        j = (i + 1) % v_len
        if (v_yy[i] <= p_y < v_yy[j]) or (v_yy[i] > p_y >= v_yy[j]):
            vt = (p_y - v_yy[i]) / (v_yy[j] - v_yy[i])
            if p_x < v_xx[i] + vt * (v_xx[j] - v_xx[i]):
                result = not result
    return result


_QUAD_VERT_WEIGHTS = np.array(
    [[1.0, 0.0, 0.0, 0.0],
     [1.0, 1.0, 0.0, 0.0],
     [1.0, 1.0, 1.0, 1.0],
     [1.0, 0.0, 1.0, 0.0]])


def calc_quad_weights(poly_xx, poly_yy, p_x, p_y):
    alphas = np.linalg.solve(_QUAD_VERT_WEIGHTS, poly_xx)
    betas = np.linalg.solve(_QUAD_VERT_WEIGHTS, poly_yy)
    aa = (alphas[3] * betas[2] - alphas[2] * betas[3])
    bb = (alphas[3] * betas[0] - alphas[0] * betas[3] + alphas[1] * betas[2] -
          alphas[2] * betas[1] + p_x * betas[3] - p_y * alphas[3])
    cc = (alphas[1] * betas[0] - alphas[0] * betas[1] +
          p_x * betas[1] - p_y * alphas[1])

    m, l = None, None

    if np.fabs(aa) < np.finfo(float).eps:
        m = -cc / bb
    else:
        det = bb * bb - 4.0 * aa * cc
        if det < 0:
            raise Exception('Negative det!!!')
        sqrt_det = np.sqrt(det)
        m_1 = (-bb + sqrt_det) / (2.0 * aa)
        m_2 = (-bb - sqrt_det) / (2.0 * aa)
        if np.isnan(m_2) or m_2 < 0 or m_2 > 1.0:
            if np.isnan(m_1) or m_1 < 0 or m_1 > 1.0:
                raise Exception('Wrong m')
            else:
                m = m_1
        else:
            if np.isnan(m_1) or m_1 < 0 or m_1 > 1.0:
                m = m_2
            else:
                m = m_1
                l_1 = ((p_x - alphas[0] - alphas[2] * m) /
                       (alphas[1] + alphas[3] * m))
                if np.isnan(l_1) or l_1 < 0 or l_1 > 1.0:
                    m = m_2
                else:
                    l_2 = ((p_x - alphas[0] - alphas[2] * m_2) /
                           (alphas[1] + alphas[3] * m_2))
                    if np.isnan(l_2) or l_2 < 0 or l_2 > 1.0:
                        l = l_1
                    else:
                        raise Exception('Ambiguous l')

    if np.isnan(m) or m < 0 or m > 1.0:
        raise Exception('Wrong m')

    if not l:
        l = (p_x - alphas[0] - alphas[2] * m) / (alphas[1] + alphas[3] * m)

    if np.isnan(l) or l < 0 or l > 1.0:
        raise Exception('Wrong l')

    w1 = (1 - l) * (1 - m)
    w2 = l * (1 - m)
    w3 = l * m
    w4 = (1 - l) * m

    return w1, w2, w3, w4
