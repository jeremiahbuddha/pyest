#! /usr/bin/env python


# Read in text file with tracking data
with open('tracking_data.txt') as f:
    tdata = f.readlines()

# Reformat data as nested python lists
out = 'TRACKING_DATA = ['
for line in tdata:
    out+='\n    [{0},{1},{2},{3}],'.format(line.split()[0],line.split()[1],
                                           line.split()[2],line.split()[3] )
out += '\n    ]'

# Write tracking data to importable module
with open('tracking_data.py','w') as f:
    f.writelines(out)
