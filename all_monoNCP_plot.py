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
from scipy.stats import gaussian_kde

#path = "/home/spark159/../../media/spark159/sw/sp_spd_tests_detail/"
path='./data/'

#GC/me/PTM analysis
#ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+ "hg19_chr1_1001win501step_anot.cn")
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+ "hg19_chr1_NCP_anot.cn")

ID_score1 = name_ID_value["data/sp_spd_tests_detail/sp7"]
ID_score2 = name_ID_value["data/sp_spd_tests_detail/sp8"]
ID_AT = name_ID_value['ATcontent']

for ID in ID_AT:
    ID_AT[ID] = 100*ID_AT[ID]

X, Y = [], []
xvalue_scores = {}
for ID in ID_score1.keys():
    xvalue = ID_AT[ID]
    score = ID_score1[ID]
    X.append(xvalue)
    Y.append(score)
    if xvalue not in xvalue_scores:
        xvalue_scores[xvalue] = []
    xvalue_scores[xvalue].append(score)
    
Xmean, Ymean, Yerr = [], [], []
for xvalue in xvalue_scores:
    if len(xvalue_scores[xvalue]) <= 1:
        continue
    Xmean.append(xvalue)
    Ymean.append(np.mean(xvalue_scores[xvalue]))
    Yerr.append(np.std(xvalue_scores[xvalue]/np.sqrt(len(xvalue_scores[xvalue]))))
 
# Calculate the point density
X, Y = np.asarray(X), np.asarray(Y)
XY = np.vstack([X,Y])
Z = gaussian_kde(XY)(XY)

# Sort the points by density, so that the densest points are plotted last
order = np.argsort(Z)
X, Y, Z = X[order], Y[order], Z[order]

fig = plt.figure()
plt.scatter(X, Y, c=Z, s=5, cmap='jet', edgecolor='', alpha=0.1)
plt.errorbar(Xmean, Ymean, yerr=Yerr, fmt='k.')
plt.plot(Xmean, Ymean,'k.')
plt.xlabel("AT content (%)")
plt.ylabel("Condensability (A.U.)")
plt.ylim([-3, 3.5])
plt.show()
plt.close()

#with open("X.pickle", "wb") as f:
#    pickle.dump(X, f)
#with open("Y.pickle", "wb") as f:
#    pickle.dump(Y, f)
#with open("Z.pickle", "wb") as f:
#    pickle.dump(Z, f)


sys.exit(1)

#graphics.Scatter_plot (ID_AT, ID_score1, ylim=[-3, 3.5], note='raw')

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
ID_score1 = name_ID_value['data/sp_spd_tests_detail/sp7']
ID_score2 = name_ID_value['data/sp_spd_tests_detail/sp8']
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

#new_ID_score2 = statis.neutralize_score_by_target(ID_score2, ID_AT)

#graphics.Scatter_plot(ID_AT, ID_score1, note='sp6')
#graphics.Scatter_plot(ID_AT, ID_score2, note='sp10')

#frac=[(4**i) for i in range(1,11)]
frac=[1 for i in range(10)]
frac = frac[::-1]


#graphics.PartitionMeanplot(ID_AT, ID_score1, ID_me, frac, note="sp6_me")
#graphics.PartitionMeanplot(ID_AT, ID_score2, ID_me, frac, note="sp7_me")
#graphics.PartitionScatterplot(ID_AT, ID_score1, ID_me, frac, note="sp6_me")
graphics.PartitionBoxplot(ID_score1, ID_CpG, frac, xlabel='CpG number', note="sp6_CpG")

for i in range(len(ID_chip_list)):
    ID_chip = ID_chip_list[i]
    name = chip_names[i]
    #graphics.PartitionMeanplot(ID_AT, ID_score1, ID_chip, frac, note="sp6_" + name)
    #graphics.PartitionMeanplot(ID_AT, ID_score2, ID_chip, frac, note="sp7_" + name)
    #graphics.PartitionScatterplot(ID_AT, ID_score1, ID_chip, frac, note="sp6_" + name)
    graphics.PartitionBoxplot(ID_score1, ID_chip, frac, xlabel=name, note="sp6_" + name)
