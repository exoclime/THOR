"""quick and dirty debugging tool to merge multiple h5 files that get dumped in a loop over columns in Alfrodull"""

import pathlib

import h5py
import numpy as np

from tools.file_helpers import get_path_matching_regex_with_groups


import argparse

parser = argparse.ArgumentParser(description='merge multiple h5 dump files.')
parser.add_argument('path', metavar='PATH', type=str,
                    help='path of data')

args = parser.parse_args()


basepath = pathlib.Path(args.path)


def merge_pattern(basepath, pattern):
    print('looking for files with pattern ' + pattern.format("(\d+)"))
    files = get_path_matching_regex_with_groups(
        basepath, pattern.format("(\d+)"))

    files = sorted(files, key=lambda d: int(d['groups'][0]))

    for f in files:
        print(f)

    # get keys from first file
    keys = []
    with h5py.File(files[0]['path'], "r") as f:
        for k in f.keys():
            keys.append((k, f[k].dtype))

    outfile = h5py.File(basepath / pattern.format("merged"), "w")

    for k, dtype in keys:

        data = []
        for f in files:

            p = f['path']
            inpt = h5py.File(p, "r")
            data.append(inpt[k][...])

        dataset = np.concatenate(data)

        dset = outfile.create_dataset(k, (len(dataset),), dtype=dtype)
        dset[...] = dataset

    outfile.close()


patterns = ['bindata_Alf_comp_trans_1-{}.h5',
            'bindata_Alf_dir_beam_trans_1-{}.h5',
            'bindata_Alf_interpTnP_1-{}.h5',
            'bindata_Alf_int_flx_1-{}.h5',
            'bindata_Alf_pop_spec_flx_thomas_1-{}.h5',
            'bindata_Alf_prep_flx_1-{}.h5',
            'bindata_Alf_prep_II_1-{}.h5']

for p in patterns:
    merge_pattern(basepath, p)
