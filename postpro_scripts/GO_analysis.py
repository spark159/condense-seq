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

def tuple_cmp(a, b):
    if a[0] < b[0]:
        return -1
    else:
        return 1

# read gene intervals
gID_field_values, field_gID_values = load_file.read_GTF ("data/Homo_sapiens.GRCh37.87.gtf", "chr1", mode="both")

gID_ginterval = {}
for gID in gID_field_values:
    TSS = gID_field_values[gID]['TSS']
    TTS = gID_field_values[gID]['TTS']
    strand = gID_field_values[gID]['strand']
    if strand == '+':
        interval = (TSS, TTS)
    else:
        interval = (TTS, TSS)
    gID_ginterval[gID] = interval

ginterval_dict = Interval_dict.double_hash(gID_ginterval, 10000, 250000000)


ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("data/hg19_chr1_171_everything_anot.cn")
ID_score2 = name_ID_value['work/condense_seq/sp10_hg19_chr1']
ID_AT = name_ID_value['ATcontent']

#new_ID_score2 = statis.neutralize_score_by_target(ID_score2, ID_AT)
#fluffy_IDs, middle_IDs, sticky_IDs = statis.quantile (new_ID_score2, 3, [5, 90, 5])

name_gID_values = {}
for ID in ID_pos:
    pos = ID_pos[ID]
    gIDs = ginterval_dict.find(pos)
    if len(gIDs) <= 0:
        continue
    for name in name_ID_value:
        if name not in name_gID_values:
            name_gID_values[name] = {}
        value = name_ID_value[name][ID]
        for gID in gIDs:
            if gID not in name_gID_values[name]:
                name_gID_values[name][gID] = []
            name_gID_values[name][gID].append(value)

name_gID_mean = {}
for name in name_gID_values:
    if name not in name_gID_mean:
        name_gID_mean[name] = {}
    for gID in name_gID_values[name]:
        assert gID not in name_gID_mean[name]
        name_gID_mean[name][gID] = np.mean(name_gID_values[name][gID])

gID_score2 = name_gID_mean['work/condense_seq/sp10_hg19_chr1']
fluffy_gIDs, middle_gIDs, sticky_gIDs = statis.quantile (gID_score2, 3, [5, 90, 5])

for ID in sticky_gIDs:
    print ID
