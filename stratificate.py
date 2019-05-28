import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math
import Interval_dict
from sklearn import linear_model
from scipy.stats import norm
from matplotlib_venn import venn3
import itertools

def all_path (N, states='ATCG'):
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output

def poly_score (seq, nts='AT', pos=False):
    num = []
    num_pos = {}
    i = 0
    while i < len(seq):
        if seq[i] in nts:
            nt = seq[i]
            count = 1
            j = i + 1
            while j < len(seq):
                if seq[j] != nt:
                    break
                count +=1
                j +=1
            num.append(count)
            if count not in num_pos:
                num_pos[count] = []
            num_pos[count].append(i)
            i = j
        else:
            i +=1
    if pos:
        return num_pos
    if len(num) == 0:
        return 0
    return max(num)

def get_dincount(seq, din=None):
    if din:
        count = 0
        for i in range(len(seq)-1):
            if seq[i:i+2].upper() == din.upper():
                count +=1
        return count
    din_count = {}
    seq = seq.upper()
    for i in range(len(seq)-1):
        din = seq[i:i+2]
        if 'N' in din:
            continue
        if din not in din_count:
            din_count[din] = 0
        din_count[din] += 1
    return din_count

def categorize (ID_value_list):
    values_IDs = {}
    IDs = ID_value_list[0].keys()
    for ID in IDs:
        values = tuple([ID_value[ID] for ID_value in ID_value_list])
        if values not in values_IDs:
            values_IDs[values] = []
        values_IDs[values].append(ID)
    return values_IDs

ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("data/hg19_chr1_171_icdseq_anot.cn")
ID_score2 = name_ID_value['work/condense_seq/sp10_hg19_chr1']
ID_seq = name_ID_value['Sequence']
ID_AT = name_ID_value['ATcontent']
ID_CpG = name_ID_value['CpGNumber']
ID_me = name_ID_value['meGCNumber']

for ID in ID_AT:
    ID_AT[ID] = ID_AT[ID]*100

# sequence dependence
ID_din_count = {}
for ID in ID_seq:
    seq = ID_seq[ID]
    din_count = get_dincount(seq)
    ID_din_count[ID] = din_count

ID_TA = {}
for ID in ID_seq:
    seq = ID_seq[ID]
    count = get_dincount(seq, din="TA")
    ID_TA[ID] = count

ID_polyAT, ID_polyGC = {}, {}
for ID in ID_seq:
    seq = ID_seq[ID]
    num_pos = poly_score(seq, nts='AT', pos=True)
    score = 0.0
    for num, pos in num_pos.items():
        if num < 3:
            continue
        score += len(pos)*(num**2)
    ID_polyAT[ID] = score
    num_pos = poly_score(seq, nts='GC', pos=True)
    score = 0.0
    for num, pos in num_pos.items():
        if num < 3:
            continue
        score += len(pos)*(num**2)
    ID_polyGC[ID] = score

print "complete to record sequence features"

AT_IDs = categorize ([ID_AT])
all_din = all_path(2)

#alpha_list = np.linspace(0.02, 1, num=len(AT_IDs.keys()))
color_list = np.linspace(0.01, 1, num=len(all_din))
marker = itertools.cycle(('.', '<', 'x', 'd', '*'))
cmap = mpl.cm.get_cmap("jet")                  

fig = plt.figure()
plt.axhline(y=0, color='k', linestyle='--')

for i in range(len(all_din)):
    din = all_din[i]
    X_axis, Y_axis = [], []
    data_count = []
    for AT in AT_IDs:
        IDs = AT_IDs[AT]
        X, Y = [], []
        for ID in IDs:
            din_count = ID_din_count[ID]
            try:
                count = din_count[din]
            except:
                count = 0
            X.append(count)
            Y.append(ID_score2[ID])
        corr = statis.get_corr (X, Y)
        X_axis.append(AT[0])
        Y_axis.append(corr)
        data_count.append(len(X))
    temp = list(cmap(color_list[i]))
    rgba_colors = []
    for j in range(len(X_axis)):
        temp[-1] = 0.01 + float(data_count[j]-min(data_count)) / (max(data_count)-min(data_count))
        rgba_colors.append(temp)
    plt.scatter(X_axis, Y_axis, s=5, color=cmap(color_list[i]), marker = marker.next(), label=din)

    
plt.xlabel("AT content (%)")
plt.ylabel("Correlation")
plt.title("Condensability VS Dinucleotide count")
plt.legend(markerscale=2, ncol=3)
plt.savefig("corr_strait_din.png", bbox_inches='tight')
#plt.colorbar()
#plt.show()
plt.close()

# Epigenetic markers
names = ["CpG me", "k4me3", "k27ac", "k9ac", "k36me3", "k9me2", "k9me3", "k27me3"]
ID_value_list = [ID_me, name_ID_value['k4me3'], name_ID_value['k27ac'], name_ID_value['k9ac'], name_ID_value['k36me3_2'], name_ID_value['k9me2_2'], name_ID_value['k9me3_2'], name_ID_value['k27me3a_2']]

color_list = np.linspace(0.01, 1, num=len(names))
marker = itertools.cycle(('.', '<', 'x', 'd', '*'))
cmap = mpl.cm.get_cmap("jet")                  

fig = plt.figure()
plt.axhline(y=0, color='k', linestyle='--')

for i in range(len(ID_value_list)):
    name = names[i]
    ID_value = ID_value_list[i]
    X_axis, Y_axis = [], []
    data_count = []
    for AT in AT_IDs:
        IDs = AT_IDs[AT]
        X, Y = [], []
        for ID in IDs:
            try:
                value = ID_value[ID]
            except:
                value = 0
            X.append(value)
            Y.append(ID_score2[ID])
        corr = statis.get_corr (X, Y)
        X_axis.append(AT[0])
        Y_axis.append(corr)
        data_count.append(len(X))
    temp = list(cmap(color_list[i]))
    rgba_colors = []
    for j in range(len(X_axis)):
        temp[-1] = 0.01 + float(data_count[j]-min(data_count)) / (max(data_count)-min(data_count))
        rgba_colors.append(temp)
    plt.scatter(X_axis, Y_axis, s=5, color=cmap(color_list[i]), marker=marker.next(), label=name)

plt.xlabel("AT content (%)")
plt.ylabel("Correlation")
plt.title("Condensability VS Epigenetic markers")
plt.legend(markerscale=2, ncol=3)
#plt.colorbar()
plt.savefig("corr_strait_epi.png", bbox_inches='tight')
#plt.show()
plt.close()

# group by same AT content and CG count

ATCG_IDs = categorize ([ID_AT, ID_CpG])

X_axis, Y_axis, Z_axis = [], [], []
data_count = []
for AT_CG in ATCG_IDs:
    IDs = ATCG_IDs[AT_CG]
    X = [ ID_me[ID] for ID in IDs ]
    Y = [ ID_score2[ID] for ID in IDs ]
    corr = statis.get_corr(X, Y)
    X_axis.append(AT_CG[0])
    Y_axis.append(AT_CG[1])
    Z_axis.append(corr)
    data_count.append(len(X))

size_list = []
for i in range(len(data_count)):
    count = data_count[i]
    size = 1000*(float(count-min(data_count)) / (max(data_count)-min(data_count)))
    size_list.append(size)
    
fig = plt.figure()
#plt.scatter(X_axis, Y_axis, s=size_list, c=Z_axis, edgecolor='k', linewidth=0.05, cmap='bwr')
plt.scatter(X_axis, Y_axis, s=5, c=Z_axis, cmap='bwr')
plt.title("Condensability VS CpG methylation")
plt.xlabel("AT content (%)")
plt.ylabel("CG count")
plt.colorbar()
plt.savefig("corr_strait_CpGme.png", bbox_inches='tight')
#plt.show()
plt.close()

fig = plt.figure()
#plt.scatter(X_axis, Y_axis, s=size_list, c=Z_axis, edgecolor='k', linewidth=0.05, cmap='bwr')
plt.scatter(X_axis, Y_axis, s=5, c=np.log10(data_count), cmap='jet')
plt.title("Sample size ($10^x$)")
plt.xlabel("AT content (%)")
plt.ylabel("CG count")
plt.colorbar()
plt.savefig("count_strait_CpGme.png", bbox_inches='tight')
#plt.show()
plt.close()
