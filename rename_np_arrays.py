#!/usr/bin/env python
from os import rename
from glob import glob

if __name__ == "__main__":
    fnames = glob("*.npy")
    for fname in fnames:
        new_name = fname.split(".")[0].zfill(8) + "." + fname.split(".")[1]
        rename(fname, new_name)
