#! /usr/bin/env python

import sys

my_file = open(sys.argv[1])

write_out = []

for line in my_file:
    write_out += [ "[ {0} ],\n".format(','.join(line.split())) ]

with open(sys.argv[1].replace('txt','py'),'w') as f:
    f.writelines(write_out)


