
import pathlib

import h5py
import numpy as np

from tools.file_helpers import get_path_matching_regex_with_groups

basepath_ref = pathlib.Path("rd_test/ref/")
basepath_cmp = pathlib.Path("rd_test_258/ref/")

patterns = ['bindata_Alf_comp_trans_1-{}.h5',
            'bindata_Alf_dir_beam_trans_1-{}.h5',
            'bindata_Alf_interpTnP_1-{}.h5',
            'bindata_Alf_int_flx_1-{}.h5',
            'bindata_Alf_pop_spec_flx_thomas_1-{}.h5',
            'bindata_Alf_prep_flx_1-{}.h5',
            'bindata_Alf_prep_II_1-{}.h5']


def compare_data(base_ref, base_cmp, fname, keys_in=None):
    keys = []
    file_ref = base_ref / fname.format("merged")
    file_cmp = base_cmp / fname.format("merged")
    with h5py.File(file_ref, "r") as f_ref:
        with h5py.File(file_cmp, "r") as f_cmp:
            if keys_in is not None:
                keys = keys_in
            else:
                keys = f_ref.keys()

            for k in keys:
                data_ref = f_ref[k][...]
                data_cmp = f_cmp[k][...]

                l_r = len(data_ref)
                l_c = len(data_cmp)

                min_l = min(l_r, l_c)
                diffs = 0
                diff_zones = []
                last_diff = -5
                diff_start = -5

                for i in range(min_l):
                    difference = data_ref[i] != data_cmp[i]
                    if difference:
                        if last_diff != i-1:
                            if diff_start != -5:
                                diff_zones.append((diff_start, last_diff))
                            diff_start = i
                        last_diff = i

                        diffs += 1
                        #                        if diffs < 100:
                        #                            print(
                        #                                f"difference[{i}]: {data_ref[i]} != {data_cmp[i]}")
                        #                        break
#                for start, end in diff_zones:
#                    print(f"\tzone: ({start}, {end}) d: {end - start}")
                print(f"dataset[{k}] - length ({l_r}, {l_c}) - diffs: {diffs}")


for p in patterns:
    print(f"comparing: {p}")
    compare_data(basepath_ref, basepath_cmp, p)
    print()

#compare_data(basepath_ref, basepath_cmp, patterns[0], keys_in=['w0u', 'w0l'])
