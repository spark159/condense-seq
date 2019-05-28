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

def read_data (fname, chr_target, score_col=6):
    ID = 0
    ID_interval, ID_score = {}, {}
    for line in open(fname):
        cols = line.strip().split()
        chr, st, ed = cols[:3]
        if chr != chr_target:
            continue
        score = cols[score_col]
        st, ed = int(st), int(ed)
        score = float(score)
        interval = (st, ed)
        assert ID not in ID_interval
        ID_interval[ID] = interval
        assert ID not in ID_score
        ID_score[ID] = score
        ID +=1
    return ID_interval, ID_score

dID_interval1, dID_DNase1 = read_data("data/ENCFF706BCL_peaks.bed", "chr1")
dID_interval2, dID_DNase2 = read_data("data/ENCFF554RWR_hotspots.bed", "chr1", score_col=4)
dinterval_dict1 = Interval_dict.double_hash(dID_interval1, 10000, 250000000)
dinterval_dict2 = Interval_dict.double_hash(dID_interval2, 10000, 250000000)

ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("data/hg19_chr1_171_everything_anot.cn")

ID_DNase1 = {}
ID_DNase2 = {}
for ID in ID_pos:
    pos = ID_pos[ID]
    dIDs1 = dinterval_dict1.find(pos)
    if len(dIDs1) <= 0:
        DNase1 = 0
    else:
        DNase1 = np.mean([dID_DNase1[dID] for dID in dIDs1])
    ID_DNase1[ID] = DNase1
    dIDs2 = dinterval_dict2.find(pos)
    if len(dIDs2) <= 0:
        DNase2 = 0
    else:
        DNase2 = np.mean([dID_DNase2[dID] for dID in dIDs2])
    ID_DNase2[ID] = DNase2

frac = [4**i for i in range(10)]
for name in name_ID_value:
    ID_value = name_ID_value[name]
    graphics.PartitionBoxplot(ID_value, ID_DNase1, frac, xlabel='DNase-seq score1', ylabel=name.split('/')[-1], note='DNase1')
    graphics.PartitionBoxplot(ID_value, ID_DNase2, frac, xlabel='DNase-seq score2', ylabel=name.split('/')[-1], note='DNase2')
