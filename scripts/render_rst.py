#! /usr/bin/env python

import sys
from os import remove
from subprocess import Popen, PIPE, STDOUT, call
from glob import glob

rst_file = sys.argv[1]
tex_file = rst_file.replace('rst','tex')
pdf_file = rst_file.replace('rst','pdf')

rst_cmd = 'rst2latex.py {0} {1} > /dev/null'.format(rst_file, tex_file)
call(rst_cmd,shell=True)
#rst_prcss = Popen(rst_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT)
#rst_prcss.wait()

tex_cmd = 'pdflatex {0} -interaction=nonstopmode > /dev/null'.format(tex_file)
call(tex_cmd,shell=True)
#tex_prcss = Popen(tex_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT)
#tex_prcss.wait()

gv_cmd = 'gv {0} > dev_null'.format(pdf_file)
call(gv_cmd,shell=True)

raw_input()

offending_extensions = ['tex', 'pdf','out','log','aux','synctex.gz']

chaff = []
for ext in offending_extensions:
    chaff += glob('*.'+ext)

for garbage in chaff:
    remove(garbage)


