import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
import Interval_dict


#GC/me/PTM analysis
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("data/hg19_chr1_171_everything_anot.cn")

"""
#divide into domains
gID_field_values, field_gID_values = load_file.read_GTF ("data/Homo_sapiens.GRCh37.87.gtf", "chr1", mode="both")
x
gID_ginterval = {}
for gID in gID_field_values.keys():
    TSS = gID_field_values[gID]['TSS']
    TTS = gID_field_values[gID]['TTS']
    strand = gID_field_values[gID]['strand']
    interval = (TSS-500, TSS+500)
    #if strand == '+':
        #interval = (TSS, TTS)
        #interval = (TSS-250, TSS)
        #interval = (TSS, TSS+2500)
    #else:
        #interval = (TTS, TSS)
        #interval = (TSS, TSS+250)
        #interval = (TSS-2500, TSS)
    gID_ginterval[gID] = interval

ginterval_dict = Interval_dict.double_hash(gID_ginterval, 10000, 250000000)

temp_name_ID_value = {}
for ID in ID_pos:
    pos = ID_pos[ID]
    gIDs = ginterval_dict.find(pos)
    if not gIDs:
        continue
    for name in name_ID_value:
        if name not in temp_name_ID_value:
            temp_name_ID_value[name] = {}
        temp_name_ID_value[name][ID] = name_ID_value[name][ID]

name_ID_value = temp_name_ID_value
"""

#ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("hg19_chr1_171_anot.cn")
ID_score1 = name_ID_value['work/condense_seq/sp9_hg19_chr1']
ID_score2 = name_ID_value['work/condense_seq/sp10_hg19_chr1']
ID_AT = name_ID_value['ATcontent']
ID_CpG = name_ID_value['CpGNumber']
ID_me = name_ID_value['meGCNumber']

ID_chip_list, chip_names = [], []
for name in name_ID_value:
    if name.startswith('k'):
        ID_value = name_ID_value[name]
        chip_names.append(name)
        ID_chip_list.append(ID_value)

for ID in ID_AT:
    ID_AT[ID] = ID_AT[ID]*100

"""
ID_meCpGfrac = {}
ID_meGCfrac = {}
for ID in ID_CpG:
    CpG = ID_CpG[ID]
    me = ID_me[ID]
    GC = 100-ID_AT[ID]
    if CpG <= 0:
        value = 0.0
    else:
        float(me)/(171*CpG)
    ID_meCpGfrac[ID] = value
    ID_meGCfrac[ID] = float(me)/GC
"""

new_ID_score2 = statis.neutralize_score_by_target(ID_score2, ID_AT)

#graphics.Scatter_plot(ID_AT, ID_score1, note='sp9')
graphics.Scatter_plot(ID_AT, ID_score2, note='sp10')

frac = [(4**i) for i in range(1,11)]
frac = frac[::-1]


#graphics.PartitionMeanplot(ID_AT, ID_score1, ID_me, frac, note="sp9_me")
graphics.PartitionMeanplot(ID_AT, ID_score2, ID_me, frac, note="sp10_me")
graphics.PartitionScatterplot(ID_AT, ID_score2, ID_me, frac, note="sp10_me")
graphics.PartitionBoxplot(new_ID_score2, ID_me, frac, note="sp10_me")

for i in range(len(ID_chip_list)):
    ID_chip = ID_chip_list[i]
    name = chip_names[i]
    #graphics.PartitionMeanplot(ID_AT, ID_score1, ID_chip, frac, note="sp9_" + name)
    graphics.PartitionMeanplot(ID_AT, ID_score2, ID_chip, frac, note="sp10_" + name)
    graphics.PartitionScatterplot(ID_AT, ID_score2, ID_chip, frac, note="sp10_" + name)
    graphics.PartitionBoxplot(new_ID_score2, ID_chip, frac, note="sp10_" + name)
