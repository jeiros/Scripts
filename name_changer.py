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
                epilog="""Changes the naming of the selected trajectory files
                 from 3-digit numbering to 4-digit numbering""")

parser.add_argument("Trajectories", help="An indefinite amount of AMBER\
                    trajectories", nargs="+")

parser.add_argument("-c", "--change", help="""Actually perform the name 
    change. Default is false.""",action="store_true")

args = parser.parse_args()



def namechange(files, pattern):
    for pathname in files:
        basename = os.path.basename(pathname)
        final_name = []
        new_filename = basename
        match = pattern.match(str(basename))
        if match is not None:
            # The pattern is split into 4 groups.
            # The groupdict() method returns a dictionary that
            # contains the value of the match for every group
            dict = match.groupdict()
            if len(dict['B']) == 4:
                # The 'B' key stores the first number with the trailing
                # - character. That's why we check if it has len(4)
                dict['B'] = dict['B'].zfill(5)
            if len(dict['C']) == 3:
                dict['C'] = dict['C'].zfill(4)
            for k, v in sorted(dict.items()):
                # Looping trough the dictionary items we can 
                # create a string concatenating the values
                # in a list and then join it to make a string
                final_name.append(v)
            new_filename = ''.join(final_name)
        if (new_filename != basename) and args.change:
            os.rename(
                pathname,
                os.path.join(os.path.dirname(pathname), new_filename))
        print(basename, "\t", "----->", "\t", new_filename, "\n")


def main():
    if args:
        # This pattern matches files that have the following name
        # ****_XXX-XXXns****
        # 'A' stores the initial characters
        # 'B' stores the first number plus the hyphen
        # 'C' stores the second number
        # 'D' stores the rest of the file name
        pattern = re.compile(r'(?P<A>.*_)(?P<B>[0-9]{3}-)(?P<C>[0-9]{3,4})(?P<D>.*)')
        namechange(args.Trajectories, pattern)

if __name__ == "__main__":
    main()
