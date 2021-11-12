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
from sklearn import linear_model
from scipy.stats import norm
import scipy.stats
from matplotlib import colors
#from matplotlib_venn import venn3
#from npeet import entropy_estimators as ee

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

path = ''
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+"H1_NCP_sp_chr1_anot.cn")
ID_score2 = name_ID_value['work/2021_06_07_H1_sp_detail/H1-NCP-sp-8']
ID_seq = name_ID_value['Sequence']
ID_AT = name_ID_value['ATcontent']
ID_CpG = name_ID_value['CNumber(CpG)']
ID_me = name_ID_value['meCNumber(CpG)']

for ID in ID_AT:
    ID_AT[ID] = ID_AT[ID]*100

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

ID_TA = {}
for ID in ID_seq:
    seq = ID_seq[ID]
    count = get_dincount(seq, din="TA")
    ID_TA[ID] = count

ID_mefrac = {}
for ID in ID_CpG:
    CpG = ID_CpG[ID]
    if CpG <= 0:
        ID_mefrac[ID] = np.NaN
        continue
    me = ID_me[ID]
    #mefrac = float(me) / (2*CpG)
    mefrac = float(me) / (CpG)
    ID_mefrac[ID] = mefrac

# collect all features
names = ["AT content", "TpA count", "CpG count", "Poly-G", "meCpG density", "H2AZ", "H3k4me3", "H3k27ac", "H3k9ac", "H3k36me3", "H3k9me3", "H3k27me3"]
ID_value_list = [ID_AT, ID_TA, ID_CpG, ID_polyGC, ID_mefrac, name_ID_value['H2AFZ'], name_ID_value['H3K4me3'], name_ID_value['H3k27ac'], name_ID_value['H3K9ac'], name_ID_value['H3K36me3'], name_ID_value['H3K9me3'], name_ID_value['H3K27me3']]


ID_state = {}
for ID in ID_seq:
    state = [ID_value[ID] for ID_value in ID_value_list]
    ID_state[ID] = state

"""
# multivariate linear regression
IDs = ID_state.keys()
feature_list = [ ID_state[ID] for ID in IDs]
test_list = [[ID_score1[ID]] for ID in IDs]
reg = linear_model.Ridge(alpha=0.5)
reg.fit (feature_list, test_list)
coef_list = reg.coef_[0]
print "Linear regression"
for i in range(len(names)):
    print names[i], coef_list[i]

fig = plt.figure()
plt.bar(range(len(coef_list)), coef_list, width=0.5, color='g')
plt.xticks(range(len(coef_list)), names, rotation=90)
plt.ylabel("Coefficient")
plt.axhline(y=0, color='k', linestyle='--')
plt.title("Multivariate linear regression with Condensability")
plt.savefig("bar_linear.png",bbox_inches='tight')
plt.show()
plt.close()

# partial correlation
print "Partial correlation"
pcorr_list = []
for i in range(len(names)):
    Zs = ID_value_list[:i] + ID_value_list[i+1:]
    ID_newvalue = statis.neutralize (ID_value_list[i], Zs)
    ID_newscore1 = statis.neutralize (ID_score1, Zs)
    A, B = [], []
    for ID in ID_newvalue:
        newvalue = ID_newvalue[ID]
        newscore1 = ID_newscore1[ID]
        A.append(newvalue)
        B.append(newscore1)
    pcorr = statis.get_corr(A, B)
    pcorr_list.append(pcorr)
    name = names[i]
    print name, pcorr

fig = plt.figure()
plt.bar(range(len(pcorr_list)), pcorr_list, width=0.5, color='g')
plt.xticks(range(len(pcorr_list)), names, rotation=90)
plt.ylabel("Pearson correlation")
plt.axhline(y=0, color='k', linestyle='--')
plt.title("Partial correlation with Condensability")
plt.savefig("bar_pcorr.png",bbox_inches='tight')
plt.show()
plt.close()
"""

# conditinoal correlation
print "Conditional correlation"
cdcorr_list = []
weights_list, corrs_list = [], []
for i in range(len(names)):
    rstate_IDs = {}
    for ID in ID_state:
        state = ID_state[ID]
        rstate = tuple(state[:i] + state[i+1:])
        if rstate not in rstate_IDs:
            rstate_IDs[rstate] = []
        rstate_IDs[rstate].append(ID)
    total = 0
    cdcorr = 0.0
    weights, corrs = [], []
    for IDs in rstate_IDs.values():
        if len(IDs) < 5:
            continue
        X = [ID_value_list[i][ID] for ID in IDs]
        Y = [ID_score2[ID] for ID in IDs]
        
        #corr = statis.get_corr(X, Y)
        #corr = scipy.stats.spearmanr(X, Y)[0]
        corr = statis.get_spearman_corr(X, Y)
        if np.isnan(corr):
            continue
        weights.append(len(IDs))
        corrs.append(corr)
        total += len(IDs)
        cdcorr += len(IDs)*corr
    if total > 0:
        cdcorr = cdcorr/float(total)
        #weights = [value/float(total) for value in weights]
    name = names[i]
    cdcorr_list.append(cdcorr)
    weights_list.append(weights)
    corrs_list.append(corrs)
    print name, cdcorr


# state vs corr list with weight (new)
color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'lime', 'salmon']

fig = plt.figure(figsize=(6,8))
for i in range(len(names)):
    corrs = np.asarray(corrs_list[i])
    weights = np.asarray(weights_list[i])
    frac_weights = weights / float(sum(weights))
    rgb_colors = np.zeros((len(corrs),4))
    rgb_colors[:,:3] = colors.to_rgba(color_list[i])[:3]
    rgb_colors[:,3] = 0.15 * weights / float(max(weights))
    order = np.argsort(weights)
    plt.scatter(corrs[order], [-i]*len(corrs), s=50000*frac_weights[order], color=rgb_colors[order])
    plt.annotate('x', (cdcorr_list[i], -i), ha='center', va='center')
plt.axvline(x=0, linestyle='--', color='k')
plt.xlim([-0.2, 0.2])
plt.xlabel("Spearman correlation")
plt.yticks([-i for i in range(len(names))], names)
plt.title("Data Stratification")
plt.tight_layout()
plt.savefig("Data_strat.png", bbox_inches='tight')
plt.show()
plt.close()


yset1, dataset1 = [], []
yset2, dataset2 = [], []
for i in range(len(names)):
    cdcorr = cdcorr_list[i]
    if cdcorr < 0:
        yset1.append(-i)
        dataset1.append(cdcorr)
    else:
        yset2.append(-i)
        dataset2.append(cdcorr)
        
fig = plt.figure(figsize=(6,8))
plt.barh(yset1, dataset1, align='center', color='tab:red', height=0.5, edgecolor='k')
plt.barh(yset2, dataset2, align='center', color='tab:blue', height=0.5, edgecolor='k')
plt.axvline(x=0, linestyle='--', color='k')
plt.xlabel("Averaged correlation")
plt.yticks([-i for i in range(len(names))], names)
plt.title("Conditional Correlation")
plt.tight_layout()
plt.savefig("Conditional_corr.png", bbox_inches='tight')
plt.show()
plt.close()


#fig = plt.figure(figsize=(5,8))
#for i in range(len(names)):
#    corrs = np.asarray(corrs_list[i])
#    weights = np.asarray(weights_list[i])
#    frac_weights = weights / float(sum(weights))
#    order = np.argsort(weights)
#    plt.scatter(corrs[order], [-i]*len(corrs), s=50000*frac_weights[order], alpha=0.2)
#    plt.annotate('x', (cdcorr_list[i], -i), ha='center', va='center')
#plt.axvline(x=0, linestyle='--', color='k')
#plt.yticks([-i for i in range(len(names))], names)
#plt.xlabel("Spearman correlation")
#plt.tight_layout()
#plt.show()
#plt.close()


## state vs corr list with weight
#fig = plt.figure()
#mweight_list = [ max(weights_list[i]) for i in range(len(weights_list)) ]
#max_weight = max(mweight_list)
#for i in range(len(names)):
#    corrs = np.asarray(corrs_list[i])
#    weights = np.asarray(weights_list[i]) #/ float(max_weight)
#    order = np.argsort(weights)
#    #plt.plot([i]*len(corrs), corrs, '.', alpha=weights)
#    plt.scatter([i]*len(corrs), corrs[order], s=100*weights[order]/max_weight, c=weights[order], cmap='#hot_r', alpha=0.2)
#    #plt.scatter([i]*len(corrs), corrs, s=100*weights, c=corrs, cmap='seismic', vmin=-1, vmax=1, alpha=#0.5)
#plt.axhline(y=0, linestyle='--', color='k')
#plt.xticks(range(len(names)), names, rotation=90)
##plt.ylabel("Pearson correlation")
#plt.ylabel("Spearman correlation")
#cbar = plt.colorbar()
#cbar.ax.set_ylabel('stratified sample size', rotation=-90, va="bottom")
#plt.tight_layout()
#plt.show()
#plt.close()

    
# state vs conditional corr
#fig = plt.figure()
#plt.bar(range(len(cdcorr_list)), cdcorr_list, width=0.5, color='g')
#plt.xticks(range(len(cdcorr_list)), names, rotation=90)
##plt.ylabel("Pearson correlation")
#plt.ylabel("Spearman correlation")
#plt.axhline(y=0, color='k', linestyle='--')
##plt.title("Conditional correlation with Condensability")
#plt.tight_layout()
##plt.savefig("bar_cdcorr.png",bbox_inches='tight')
#plt.show()
#plt.close()


"""
# conditional mutual information
print "Conditional mutual information"

state_IDs = {}
for ID in ID_state:
    state = tuple(ID_state[ID])
    if state not in state_IDs:
        state_IDs[state] = []
    state_IDs[state].append(ID)

total = 0
H_Y_allX = 0.0
for IDs in state_IDs.values():
    if len(IDs) < 10:
        continue
    Y = [ID_score1[ID] for ID in IDs]
    total += len(IDs)
    H_Y_allX += len(IDs)*ee.entropyd(Y)
H_Y_allX = H_Y_allX / total

cdminfo_list = []
for i in range(len(names)):
    rstate_IDs = {}
    for ID in ID_state:
        state = ID_state[ID]
        rstate = tuple(state[:i] + state[i+1:])
        if rstate not in rstate_IDs:
            rstate_IDs[rstate] = []
        rstate_IDs[rstate].append(ID)
    total = 0
    H_Y_rX = 0.0
    for IDs in rstate_IDs.values():
        if len(IDs) < 10:
            continue
        Y = [ID_score1[ID] for ID in IDs]
        total += len(IDs)
        H_Y_rX += len(IDs)*ee.entropyd(Y)
    H_Y_rX = H_Y_rX / total
    cdminfo = H_Y_rX - H_Y_allX
    #assert cdminfo >= 0
    cdminfo_list.append(cdminfo)
    print names[i], cdminfo

fig = plt.figure()
plt.bar(range(len(cdminfo_list)), cdminfo_list, width=0.5, color='g')
plt.xticks(range(len(cdminfo_list)), names, rotation=90)
plt.ylabel("Bits")
plt.axhline(y=0, color='k', linestyle='--')
plt.title("Conditional mutual information with Condensability")
plt.savefig("bar_cdminfo.png",bbox_inches='tight')
plt.show()
plt.close()
"""
