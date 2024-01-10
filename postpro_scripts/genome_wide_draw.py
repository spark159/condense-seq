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
#ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("data/hg19_chr1_167win25step_new_anot.cn", num_max=1000000)

#target_names = ["work/condense_seq/sp10_hg19_chr1", "ATcontent"]
#target_names = ["work/condense_seq/sp10_hg19_chr1", "CpGNumber"]
#target_names = ["work/condense_seq/sp10_hg19_chr1", "meGCNumber"]
#target_names = ["work/condense_seq/sp10_hg19_chr1", "k27ac"]
#target_names = ["work/condense_seq/sp10_hg19_chr1", "k27me3"]
#target_names = ["work/condense_seq/sp10_hg19_chr1", "k36me3"]
#target_names = ["work/condense_seq/sp10_hg19_chr1", "k4me3"]
#target_names = ["work/condense_seq/sp10_hg19_chr1", "k9ac"]
#target_names = ["work/condense_seq/sp10_hg19_chr1", "k9me2"]
#target_names = ["work/condense_seq/sp10_hg19_chr1", "k9me3"]
target_names = ["work/condense_seq/sp1_hg19_chr1", "work/condense_seq/sp10_hg19_chr1"]
target_names = ["data/sp_spd_tests_detail/sp7", "data/sp_spd_tests_detail/sp8",	"data/sp_spd_tests_detail/sp1"]
#target_names = ["ATcontent"]
#ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("data/hg19_chr1_167win25step_new_anot.cn", target_names=target_names)
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("/home/spark159/../../media/spark159/sw/sp_spd_tests_detail/hg19_chr1_167win25step_cov_Bsig.cn", target_names=target_names)
print "reading done"

#posID = [[pos, ID] for ID, pos in ID_pos.items()]
#posID = sorted(posID, cmp=tuple_cmp)
#ID_list = [ value[1] for value in posID ]

#ID_score2 = name_ID_value["work/condense_seq/sp10_hg19_chr1"]
#graphics.draw_along_genome_pair (ID_pos, ID_score2, name_ID_value[target_names[1]], 100, ylabel1='Condensability (A.U.)', ylabel2=target_names[1], note="Condensability_" + target_names[1])
graphics.draw_along_genome_pair (ID_pos, name_ID_value[target_names[0]], name_ID_value[target_names[1]], 1000,  ylabel1=target_names[0], ylabel2=target_names[1], note=target_names[0] + "_" + target_names[1])

#graphics.draw_along_genome (ID_pos, [name_ID_value[target_names[0]]], 100, target_names, ylabel="", note="")

#for name in name_ID_value:
#    ID_value = name_ID_value[name]
#    graphics.draw_along_genome (ID_pos, [ID_value], 10000, labels = [name], ylabel="", note=name)


"""
fig = plt.figure()
for name in name_ID_value:
    if name == "work/condense_seq/sp9_hg19_chr1":
        continue
    ID_value = name_ID_value[name]
    value_list = [ID_value[ID] for ID in ID_list]
    acf_list = statis.acf(value_list)
    X = np.asarray([i*25 for i in range(len(acf_list))])
    if name == "work/condense_seq/sp10_hg19_chr1":
        name = "Condensability"
    plt.plot(np.log10(X+1), acf_list, label=name, alpha=0.8)
plt.title("chr1")
plt.axhline(y=0, color='k', linestyle='--')
plt.xlabel("Lag (log(bp+1))")
plt.ylabel("Autocorrelation")
plt.legend()
plt.show()
plt.close()



ID_chip_list, chip_names = [], []
for name in name_ID_value:
    if name.startswith('k'):
        ID_value = name_ID_value[name]
        chip_names.append(name)
        ID_chip_list.append(ID_value)

ID_AT = {}
for ID in ID_AT_raw:
    ID_AT[ID] = ID_AT_raw[ID] - 0.5

ID_me_r = {}
for ID in ID_me:
    ID_me_r[ID] = ID_me[ID]/20.0

ID_meden = {}
ID_meden2 = {}
for ID in ID_CpG:
    CpG = ID_CpG[ID]
    me = ID_me[ID]
    if CpG <= 0:
        value = 0.0 -0.5
    else:
        value = (float(me) / CpG) - 0.5
    ID_meden[ID] = value
    ID_meden2[ID] = me/(10*(1-ID_AT[ID]))

ID_k9me3 = {}
for ID in ID_k9me3_raw:
    k9me3 = ID_k9me3_raw[ID]
    if k9me3 > 100:
        value = 1
    else:
        value = 0
    ID_k9me3[ID] = value
"""
