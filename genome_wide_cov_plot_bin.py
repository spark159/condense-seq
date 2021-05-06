import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
from sklearn import linear_model
import random

def read_bincountfile (fname, chr_list=None):
    First = True
    for line in open(fname):
        if First:
            cols = line.strip().split()
            names = [name.rsplit('.')[-2] for name in cols[4:-1]]
            chr_binID_counts = [{} for i in range(len(names))]
            chr_binID_range = {}
            chr_binID_GC = {}
            First = False
            continue
        line = line.strip()
        if not line:
            continue
        cols = line.strip().split()
        ID, chr, st, ed = cols[:4]
        if chr_list != None and chr not in chr_list:
            continue
        ID = int(ID)
        st, ed = int(st), int(ed)
        GC = float(cols[-1])
        if chr not in chr_binID_range:
            chr_binID_range[chr] = []
        chr_binID_range[chr].append((st, ed))
        if chr not in chr_binID_GC:
            chr_binID_GC[chr] = []
        chr_binID_GC[chr].append(GC)
        datas = cols[4:-1]
        for i in range(len(datas)):
            data = float(datas[i])
            chr_binID_count = chr_binID_counts[i]
            if chr not in chr_binID_count:
                chr_binID_count[chr] = []
            chr_binID_count[chr].append(data)
    return names, chr_binID_counts, chr_binID_range, chr_binID_GC

path = ""
#fname = "H1_NCP-new_spd_10kb_bin.cn"
fname = "NCP_Spermine(4+)_10kb_bin.cn"
names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname, chr_list=['chr1'])
chr_binID_control = chr_binID_counts[-1]

ID_pos = {}
ID_count = {}
for i in range(len(chr_binID_control['chr1'])):
    start, end = chr_binID_range['chr1'][i]
    pos = (start+end)/2
    ID_pos[i] = pos
    ID_count[i] = chr_binID_control['chr1'][i]

ID_counts = chr_binID_control['chr1']
graphics.draw_along_genome (ID_pos, [ID_count], win=1000, labels=['Input'], ylabel='Read counts per 10kb', ylim=[-5,100], title='Chromosome 1', scatt=True, note="_readcounts")

