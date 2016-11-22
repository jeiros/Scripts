#!/usr/bin/env python
import re
import os
import argparse
import sys

parser = argparse.ArgumentParser(
    usage="{} Files*[0-9].dat".format(sys.argv[0]),
    epilog="""Changes the naming of the selected files from
         1,2-digit numbering to 3-digit numbering""")

parser.add_argument("Files", help="An indefinite amount of files", nargs="+")

parser.add_argument("-c", "--change", help="""Actually perform the name
                    change. Default is false.""", action="store_true")

args = parser.parse_args()


def namechange(files, pattern):
    for pathname in files:
        basename = os.path.basename(pathname)
        final_name = []
        new_filename = basename
        match = pattern.match(str(basename))
        if match is not None:
            # The pattern is split into 3 groups.
            # The groupdict() method returns a dictionary that
            # contains the value of the match for every group
            dict = match.groupdict()
            dict['B'] = dict['B'].zfill(3)  # Convert to 3-digit numbers
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
        else:
            print(basename, "\t", "----->", "\t", new_filename, "\n")


def main():
    if args:
        # This pattern matches files that have the following name
        # ****[0-9].dat
        # 'A' stores the initial characters
        # 'B' stores the number
        # 'C' stores the .dat string
        pattern = re.compile(
            r'(?P<A>.*?(?=[0-9]+.dat))(?P<B>[0-9]*?(?=.dat))(?P<C>.dat)')
        namechange(args.Files, pattern)
if __name__ == "__main__":
    main()
