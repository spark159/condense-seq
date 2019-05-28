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

def tuple_cmp (a, b):
    if a[0] < b[0]:
        return -1
    else:
        return 1

def read_Ncov(fname):
    ID_pos = {}
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            label = cols[3:]
            ID_sig_list = [{} for i in range(len(label))]
            First = False
            continue
        ID, chr, pos = cols[:3]
        ID_pos[ID] = int(pos)
        sig_list = cols[3:]
        for i in range(len(sig_list)):
            sig = sig_list[i]
            ID_sig_list[i][ID] = float(sig)/167.0
    return ID_pos, ID_sig_list

#ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("data/hg19_chr1_167win25step_anot.cn", num_max=1000000)
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("data/hg19_chr1_167win25step_new_anot.cn", num_max=10000)
#ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("hg19_chr1_171_everything_anot.cn")
print "reading done"

ID_value = name_ID_value["work/condense_seq/sp10_hg19_chr1"]
for i in range(5):
    graphics.draw_along_genome (ID_pos, [ID_value], 10**i, labels=["Win_size:" + str(10**i)], ylabel="", scatt=True)
