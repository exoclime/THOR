import numpy as np
import h5py
import subprocess as spr
import argparse

# code to test the state of thor by running and checking outputs of config files
parser = argparse.ArgumentParser()
parser.add_argument("-n","--num_steps",nargs=1,default=[1920],help='number of time steps')
parser.add_argument("-o","--outdir",nargs=1,default=['testing'],help='output directory for testing model benchmarks and features'])
args = parser.parse_args()

write = [{'config':'ifile/earth_hstest.thr','num_steps':args.num_steps[0],
                'outdir':args.outdir+'/write/earth_hs_rest'},
         {'config':'ifile/earth_hstest.thr','num_steps':args.num_steps[0],
                'outdir':args.outdir+'/write/earth_hs_norest'},
         {'config':'ifile/earth_sync.thr','num_steps':args.num_steps[0],
                'outdir':args.outdir+'/write/earth_sync'},
         {'config':'ifile/shallowhj.thr','num_steps':args.num_steps[0],
                'outdir':args.outdir+'/write/shallowhj'},
         {'config':'ifile/deephj.thr','num_steps':args.num_steps[0],
                'outdir':args.outdir+'/write/deephj'},
         {'config':'ifile/wasp43b_ex.thr','num_steps':args.num_steps[0],
                'outdir':args.outdir+'/write/wasp43'}]
