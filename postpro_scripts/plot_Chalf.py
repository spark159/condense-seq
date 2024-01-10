import sys
import copy
import matplotlib.pyplot as plt
import numpy as np
import math
import random
from sklearn import linear_model
import matplotlib.cm as cm
import statis

def read_bin_Chalf (fname, chr_choices=None):
    chr_binID_value = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        if First:
            First = False
            continue
        cols = line.split('\t')
        binID, chr, st, ed, Chalf = cols
        if chr_choices and chr not in chr_choices:
            continue
        st, ed = int(st), int(ed)
        binID = (st, ed)
        Chalf = float(Chalf)
        if chr not in chr_binID_value:
            chr_binID_value[chr] = {}
        chr_binID_value[chr][binID] = Chalf
    return chr_binID_value

# check Chalf
chr_binID_Chalf1 = read_bin_Chalf("mCD8T_WT-NCP_sp_10kb_Chalf.cn")
chr_binID_Chalf2 = read_bin_Chalf("mCD8T_inht-NCP_sp_10kb_Chalf.cn")
chr_binID_Chalf3 = read_bin_Chalf("mCD8T_KO-NCP_sp_10kb_Chalf.cn")

chr_field_binID_Chalf = {}
for chr in chr_binID_Chalf1:
    binID_Chalf1 = chr_binID_Chalf1[chr]
    binID_Chalf2 = chr_binID_Chalf2[chr]
    binID_Chalf3 = chr_binID_Chalf3[chr]
    binIDs = list(set(binID_Chalf1) & set(binID_Chalf2) & set(binID_Chalf3))
    field_binID_Chalf = {'WT-Chalf':{}, 'inht-Chalf':{}, 'KO-Chalf':{}}
    for binID in binIDs:
        Chalf1 = binID_Chalf1[binID]
        Chalf2 = binID_Chalf2[binID]
        Chalf3 = binID_Chalf3[binID]
        field_binID_Chalf['WT-Chalf'][binID] = Chalf1
        field_binID_Chalf['inht-Chalf'][binID] = Chalf2
        field_binID_Chalf['KO-Chalf'][binID] = Chalf3
    chr_field_binID_Chalf[chr] = field_binID_Chalf

del chr_binID_Chalf1
del chr_binID_Chalf2
del chr_binID_Chalf3


# compute acf
# select the longest data list with constant step size
chr_stID_length = {}
for chr in chr_field_binID_Chalf:
    binID_Chalf = chr_field_binID_Chalf[chr]['WT-Chalf']
    binIDs = sorted(binID_Chalf.keys())
    i = 0
    while i < len(binIDs) - 1:
        if binIDs[i][1] == binIDs[i+1][0]:
            j = 1
            while i + j < len(binIDs) - 1:
                if binIDs[i+j][1] != binIDs[i+j+1][0]:
                    if chr not in chr_stID_length:
                        chr_stID_length[chr] = {}
                    chr_stID_length[chr][binIDs[i]] = j + 1
                    break
                j +=1
            if i + j >= len(binIDs) - 1:
                if chr not in chr_stID_length:
                    chr_stID_length[chr] = {}
                chr_stID_length[chr][binIDs[i]] = j + 1
            i +=j
        i +=1

# select longest data stretch
chr_field_stretch = {}
for chr in chr_stID_length:
    stID_length = chr_stID_length[chr]
    length, stID = sorted([(length, stID) for stID, length in stID_length.items()], reverse=True)[0]
    st, ed = stID
    binstep = ed - st
    
    # sanity check
    binID = (st - binstep, ed - binstep)
    assert binID not in chr_field_binID_Chalf[chr]['WT-Chalf']
    for i in range(length):
        binID = (st + i*binstep, ed + i*binstep)
        assert binID in chr_field_binID_Chalf[chr]['WT-Chalf']
    binID = (st + length*binstep, ed + length*binstep)
    assert binID not in chr_field_binID_Chalf[chr]['WT-Chalf']

    # get data stretch
    for field in chr_field_binID_Chalf[chr]:
        binID_Chalf = chr_field_binID_Chalf[chr][field]
        for i in range(length):
            binID = (st + i*binstep, ed + i*binstep)
            Chalf = binID_Chalf[binID]
            if chr not in chr_field_stretch:
                chr_field_stretch[chr] = {}
            if field not in chr_field_stretch[chr]:
                chr_field_stretch[chr][field] = []
            chr_field_stretch[chr][field].append(Chalf)

# plot acfs
for chr in ['chr1']:
    fig = plt.figure()
    for field in chr_field_stretch[chr]:
        acf = statis.acf(chr_field_stretch[chr][field])
        X = np.asarray([i*binstep for i in range(len(acf))])
        plt.plot(X[1:], acf[1:], label=field, alpha=0.7, lw=2)
    plt.xlabel("Distance (bp)", fontsize=8)
    plt.ylabel("Auto-correlation", fontsize=8)
    plt.xscale("log")
    #plt.show()
    plt.close()


# Chalf histogram
data1, data2, data3 = [], [], []
for chr in chr_field_binID_Chalf:
    field_binID_Chalf = chr_field_binID_Chalf[chr]
    data1 += field_binID_Chalf['WT-Chalf'].values()
    data2 += field_binID_Chalf['inht-Chalf'].values()
    data3 += field_binID_Chalf['KO-Chalf'].values()

fig = plt.figure()
#plt.violinplot([data1, data2, data3], range(3), quantiles=[[0.05, 0.1, 0.8, 0.9]]*3, showmedians=True, widths=0.3)
plt.boxplot([data1, data2, data3], range(3))

#plt.hist(data1, bins=100, density=True, stacked=True, histtype='step', label='WT-Chalf')
#plt.hist(data2, bins=100, density=True, stacked=True, histtype='step', label='inht-Chalf')
#plt.hist(data3, bins=100, density=True, stacked=True, histtype='step', label='KO-Chalf')
plt.legend()
#plt.show()
plt.close()
