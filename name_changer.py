#!/usr/bin/env python

"""
Use to rename MD trajectory files
to 4 digit numbering, to account
for trajectories that go past 1000 ns
"""

import re
import os
import argparse
import sys
parser = argparse.ArgumentParser(
                usage="{} Trajectories*.nc".format(sys.argv[0]),
                epilog="Changes the naming of the selected trajectory files from 3-digit numbering to 4-digit numbering")

parser.add_argument("Trajectories", help="An indefinite amount of AMBER\
                    trajectories", nargs="+")
parser.add_argument("-c", "--change", help="Actually perform the name change. Default is false.",action="store_true")
args = parser.parse_args()



def namechange(files, pattern):
    for pathname in files:
        basename = os.path.basename(pathname)
        final_name = []
        new_filename = basename
        match = pattern.match(str(basename))
        if match is not None:
            dict = match.groupdict()
            if len(dict['B']) == 4:
                dict['B'] = "0" + dict['B']
            if len(dict['C']) == 3:
                dict['C'] = "0" + dict['C']
            for k, v in sorted(dict.items()):
                final_name.append(v)
            new_filename = ''.join(final_name)
            print(new_filename)
        if (new_filename != basename) and args.change:
            os.rename(
                pathname,
                os.path.join(os.path.dirname(pathname), new_filename))
        else:
            print(basename, "\t", new_filename)


def main():
    if args:
        pattern = re.compile(r'(?P<A>.*_)(?P<B>[0-9]{3}-)(?P<C>[0-9]{3,4})(?P<D>.*)')
        namechange(args.Trajectories, pattern)

if __name__ == "__main__":
    main()
