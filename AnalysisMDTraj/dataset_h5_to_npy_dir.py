#!/usr/bin/env python
from msmbuilder.dataset import dataset
from msmbuilder.io import save_trajs, load_meta
import argparse
parser = argparse.ArgumentParser(prog='dataset_h5_to_npy_dir.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''version1''')

parser.add_argument("dataset", help="""An HDF5 dataset""", type=str)
parser.add_argument("meta", help="A metadata pickl file", type=str)
parser.add_argument("trajs", help="The folder in which to store the trajs",
                    type=str, default='trajs')

if __name__ == '__main__':
    args = parser.parse_args()
    meta = load_meta(args.meta)
    ds = dataset(args.dataset)
    trajs = {}
    for k, v in ds.items():
        trajs[k] = v
    save_trajs(trajs, args.trajs, meta)
