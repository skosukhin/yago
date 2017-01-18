import shutil
import sys

import numpy as np
from netCDF4 import Dataset

import cmd.common.name_constants as names
from cmd.common.misc import create_dir_for_file
from cmd.common.nc_utils import add_or_append_history
from cmd.common.arg_processors import ListParser

description = 'adjusts fields to a given mask'

_MAX_NEIGHBOUR_COUNT = 9
_DEFAULT_BACKUP_VALUE = 0.0
_DEFAULT_MAX_ITER_COUNT = 1000


def setup_parser(parser):
    list_parser = ListParser()
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--data-var-names', required=True,
                        type=list_parser)
    parser.add_argument('--mask-file', required=True)
    parser.add_argument('--output-file', required=True)
    parser.add_argument('--backup-values', type=list_parser)
    parser.add_argument('--max-iter-count', type=np.intp,
                        default=_DEFAULT_MAX_ITER_COUNT)


def cmd(args):
    if args.backup_values is not None and \
                    len(args.backup_values) != len(args.data_var_names):
        raise Exception()

    mask_ds = Dataset(args.mask_file, 'r')
    mask = np.ma.getmaskarray(mask_ds.variables[names.VAR_MASK][:])
    mask_ds.close()

    create_dir_for_file(args.output_file)
    shutil.copyfile(args.input_file, args.output_file)

    ds = Dataset(args.output_file, 'r+')

    for i, data_var_name in enumerate(args.data_var_names):
        print data_var_name
        data_var = ds.variables[data_var_name]

        if args.backup_values is None:
            backup_value = _DEFAULT_BACKUP_VALUE
        else:
            backup_value = args.backup_values[i]

        backup_value = data_var.dtype.type(backup_value)

        if len(data_var.dimensions) == 2:
            _progress(0, 1)
            data_var[:] = _apply_mask_fast(data_var[:], mask, backup_value)
            _progress(1, 1)
        # Assume that time is the first dimension
        elif len(data_var.dimensions) == 3:
            for time_idx in xrange(0, data_var.shape[0]):
                _progress(time_idx, data_var.shape[0])
                data_var[time_idx, :] = _apply_mask_fast(data_var[time_idx, :],
                                                         mask, backup_value)
            _progress(data_var.shape[0], data_var.shape[0])
        else:
            raise Exception()

    add_or_append_history(ds)
    ds.close()


def _progress(row_num, row_count):
    print '%.2f%%' % (float(row_num) / float(row_count) * 100.0)
    sys.stdout.flush()


def _apply_mask_fast(field, mask, backup_value,
                     max_iter_count=_DEFAULT_MAX_ITER_COUNT, copy=False):
    """
        Walk order dependent algorithm for mask adjustment.
    """
    result = np.ma.masked_where(mask, field, copy=copy)
    result.unshare_mask()
    result_mask = np.ma.getmaskarray(result)

    iter_count = 0
    while iter_count < max_iter_count:
        none_fixed = True
        for i in xrange(result.shape[0]):
            for j in xrange(result.shape[1]):
                if not mask[i, j] and result_mask[i, j]:
                    left_i = i - 1
                    if left_i < 0:
                        left_i = 0

                    right_i = i + 2
                    if right_i > result.shape[0]:
                        right_i = result.shape[0]

                    bottom_j = j - 1
                    if bottom_j < 0:
                        bottom_j = 0

                    top_j = j + 2
                    if top_j > result.shape[1]:
                        top_j = result.shape[1]

                    new_value = \
                        np.ma.mean(result[left_i:right_i, bottom_j:top_j])

                    if new_value is not np.ma.masked:
                        result[i, j] = new_value
                        none_fixed = False

        if none_fixed:
            break

        iter_count += 1

    for i in xrange(result.shape[0]):
        for j in xrange(result.shape[1]):
            if not mask[i, j] and result_mask[i, j]:
                result[i, j] = backup_value

    return result


def _apply_mask(field, mask, backup_value,
                max_iter_count=_DEFAULT_MAX_ITER_COUNT, copy=False):
    """
        Failed attempt to implement walk order independent algorithm. WIP.
    """
    result = np.ma.masked_where(mask, field, copy=copy)
    result.unshare_mask()

    result_mask = np.ma.getmaskarray(result)

    min_masked_count = _MAX_NEIGHBOUR_COUNT
    min_masked_indices = []
    ignore_mask = np.zeros(result.shape, dtype=bool)

    for k in xrange(max_iter_count):
        for i in xrange(result.shape[0]):
            for j in xrange(result.shape[1]):
                if not mask[i, j] and \
                        result_mask[i, j] and \
                        not ignore_mask[i, j]:
                    left_i = i - 1
                    if left_i < 0:
                        left_i = 0

                    right_i = i + 2
                    if right_i > result.shape[0]:
                        right_i = result.shape[0]

                    bottom_j = j - 1
                    if bottom_j < 0:
                        bottom_j = 0

                    top_j = j + 2
                    if top_j > result.shape[1]:
                        top_j = result.shape[1]

                    masked_count = np.ma.count_masked(
                        result[left_i:right_i, bottom_j:top_j])

                    if masked_count < min_masked_count:
                        min_masked_count = masked_count
                        min_masked_indices = []

                    if min_masked_count == masked_count:
                        min_masked_indices.append(
                            (i, j, left_i, right_i, bottom_j, top_j))

        if len(min_masked_indices) == 0:
            break

        for i, j, left_i, right_i, bottom_j, top_j in min_masked_indices:
            new_value = np.ma.mean(result[left_i:right_i, bottom_j:top_j])

            if new_value is not np.ma.masked:
                result[i, j] = new_value
                result_mask[i, j] = False
            else:
                ignore_mask[i, j] = True

        min_masked_count = _MAX_NEIGHBOUR_COUNT
        min_masked_indices = []

    for i in xrange(result.shape[0]):
        for j in xrange(result.shape[1]):
            if not mask[i, j] and result_mask[i, j]:
                result[i, j] = backup_value

    return result
