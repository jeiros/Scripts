#!/usr/bin/env python

"""
Use to rename MD trajectory files
to 4 digit numbering, to account
for trajectories that go past 1000 ns

"""




import re, glob, os


filenames = sorted(glob.glob("05*.nc"))


pattern = re.compile(r'(?P<A>.*_)(?P<B>[0-9]{3}-)(?P<C>[0-9]{3,4})(?P<D>.*)')

for pathname in filenames:
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
    print(basename, "\t", new_filename)
    # if new_filename != basename:
    #     os.rename(
    #         pathname,
    #         os.path.join(os.path.dirname(pathname), new_filename))