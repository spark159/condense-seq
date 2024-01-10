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
from pyliftover import LiftOver

def read_ATAC (fname, chr_choice):
    lID_score = {}
    lID_interval = {}
    for line in open(fname):
        cols = line.strip().split()
        chr, st, ed, id, _, _, counts, score, _ = cols
        if chr != chr_choice:
            continue
        st, ed = int(st), int(ed)
        counts = int(counts)
        #score = float(score)
        score = np.log2(counts+1.0)
        lID_score[id] = score
        lID_interval[id] = (st, ed)
    return lID_score, lID_interval

path = "./data/"
    
lID_score, lID_interval = read_ATAC (path + "GM12878_ATACseq_50k_AllReps_ZINBA_pp08.bed", chr_choice="chr1")
linterval_dict = Interval_dict.double_hash(lID_interval, 100000, 250000000)

ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path + "hg19_chr1_NCP_anot.cn")

ID_score = name_ID_value['data/sp_spd_tests_detail/sp7']

ID_AT = name_ID_value['ATcontent']
for ID in ID_AT:
    ID_AT[ID] = ID_AT[ID]*100

ID_CG = name_ID_value['CpGNumber']
ID_me = name_ID_value['meGCNumber']
ID_mefrac = {}
for ID in ID_CG:
    CG = ID_CG[ID]
    if CG <= 0:
        continue
    me = ID_me[ID]
    mefrac = float(me) / (2*CG)
    ID_mefrac[ID] = mefrac
name_ID_value['meCpG density'] = ID_mefrac
#names.append('meCpG density')

#names = name_ID_value.keys()
names = ['meGCNumber', 'k27me3a', 'k27ac', 'k9me3', 'k9me2', 'ATcontent', 'k36me3', 'meCpG density', 'k4me3', 'k9ac', 'CpGNumber', 'data/sp_spd_tests_detail/sp7']

ID_ATAC = {}
for ID in ID_pos:
    pos = ID_pos[ID]
    lIDs = linterval_dict.find(pos)
    if len(lIDs) <= 0:
        ATAC = 0.0
    else:
        ATAC = np.mean([lID_score[lID] for lID in lIDs])
    ID_ATAC[ID] = ATAC


frac=[(4**i) for i in range(1,11)]
#graphics.PartitionScatterplot (ID_AT, ID_score, ID_ATAC, frac=frac, xlim=[0, 100], ylim=[-3.5, 3.5], xlabel="AT content (%)", ylabel = 'Condensability (A.U.)', title="", note="")

X = [ ID_ATAC[ID] for ID in sorted(ID_pos.keys()) ]
Y_list = []
for i in range(len(names)):
    name = names[i]
    ID_value = name_ID_value[name]
    Y = []
    for ID in sorted(ID_pos.keys()):
        try:
            value = ID_value[ID]
        except:
            value = np.nan
        Y.append(value)
    Y_list.append(Y)

corr_list = []
for i in range(len(names)):
    name = names[i].split('/')[-1]
    Y = Y_list[i]
    corr = statis.get_corr(X, Y)
    corr_list.append(corr)
    print name, corr
    
names[-1] = "Condensability"
fig = plt.figure()
plt.bar(range(len(corr_list)), corr_list, width=0.5, color='g')
plt.xticks(range(len(corr_list)), names, rotation=90)
plt.ylabel("Pearson Correlation")
plt.axhline(y=0, color='k', linestyle='--')
plt.title("Correlation with ATAC-seq counts")
plt.savefig("bar_corr.png",bbox_inches='tight')
#plt.show()
plt.close()
