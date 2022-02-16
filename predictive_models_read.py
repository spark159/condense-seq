import sys
import copy
import math
import random
import pickle

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sparse
from scipy.stats import pearsonr

from sklearn import linear_model
from sklearn.svm import SVR
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor


def read_anot_file(fname, target_names=None, jump=None, num_max=sys.maxsize):
    ID_chr, ID_pos = {}, {}
    name_ID_value = {}
    First = True
    count = 0
    counter = -1
    for line in open(fname):
        if count > num_max:
            break
        cols = line.strip().split()
        if First:
            names = cols[3:]
            First = False
            continue
        counter += 1
        if jump and counter % jump != 0:
            continue
        ID, chr, pos = int(cols[0]), cols[1], int(cols[2])
        ID_chr[ID] = chr
        ID_pos[ID] = pos
        cols = cols[3:]
        for i in range(len(cols)):
            name = names[i]
            if target_names and name not in target_names:
                continue
            if name not in name_ID_value:
                name_ID_value[name] = {}
            assert ID not in name_ID_value[name]
            try:
                value = float(cols[i])
            except:
                value = cols[i]
            if value == 'NA':
                value = np.NaN
            name_ID_value[name][ID] = value
        count += 1
    return ID_chr, ID_pos, name_ID_value


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

def get_density (ID_CNum, ID_meNum):
    ID_mefrac = {}
    for ID in ID_CNum:
        CNum = ID_CNum[ID]
        if CNum <= 0:
            ID_mefrac[ID] = 0.0
            continue
        meNum = ID_meNum[ID]
        mefrac = float(meNum) / (CNum)
        ID_mefrac[ID] = mefrac
    return ID_mefrac

# select feature sets
names = ['AT content', 'poly-G/C length', 'meCpG density', 'meCHG density', 'meCHH density', 'H2AFZ', 'H2AK5ac', 'H2BK120ac', 'H2BK12ac', 'H2BK15ac', 'H2BK20ac', 'H2BK5ac', 'H3K14ac', 'H3K18ac', 'H3K23ac', 'H3K23me2', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K56ac', 'H3K79me1', 'H3K79me2', 'H3K9ac', 'H3K9me3', 'H4K20me1', 'H4K5ac', 'H4K8ac', 'H4K91ac']


# load data
model_corrs = pickle.load(open('model_corrs.pickle', "rb"))
model_importances = pickle.load(open('model_importances.pickle', "rb"))


models = ["Linear reg", "SVR", "Boosting", "Random Forest", "Neural Network"]
for model in models:
    print (model, np.mean(model_corrs[model]), file=sys.stderr)

# plot prediction vs experiment
for model in models:
    X, Y = pickle.load(open("ExpVSPred_%s.pickle" % (model), "rb"))
    fig = plt.figure(figsize=(2.8, 2))
    plt.plot(X, Y, 'k.', ms=1)
    plt.xlabel("Experiment", fontsize=8)
    plt.ylabel("Prediction", fontsize=8)
    plt.gca().tick_params(axis='both', which='major', labelsize=6)
    plt.gca().tick_params(axis='both', which='minor', labelsize=6)
    minvalue, maxvalue = min(X+Y), max(X+Y)
    plt.plot([minvalue, 2.5], [minvalue, 2.5], 'r--', lw=1)
    #plt.xlim([minvalue, maxvalue])
    #plt.ylim([minvalue, maxvalue])
    plt.title("Model prediction (%s)" % (model), fontsize=10)
    plt.savefig('ExpVSPre_%s.png' % (model), dpi=1000, bbox_inches='tight')
    plt.close()


# plot correlations 
color_list = ['tab:red', 'tab:blue', 'tab:orange', 'tab:green', 'tab:purple']
fig = plt.figure(figsize=(2.8, 2))
for i in range(len(models)):
    model = models[i]
    corrs = model_corrs[model]
    bp = plt.boxplot(corrs, positions=[i], patch_artist=True, showfliers=False, widths=[0.4])
    for patch in bp['boxes']:
        patch.set_facecolor(color_list[i])
    for element in ['whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color='black')
#plt.boxplot(corrs_list, positions=range(len(mname_list)))
plt.xticks(range(len(models)), models, rotation=45, fontsize=8,  ha="right", rotation_mode="anchor")
plt.gca().tick_params('y', labelsize=6)
plt.ylabel("Pearson correlation", fontsize=8)
plt.title("10-fold cross validation", fontsize=10)
plt.xlim([-0.5, len(models)-0.5])
plt.ylim([0.55, 0.62])
plt.savefig('Pearson_model.svg', format='svg', bbox_inches='tight')
#plt.show()
plt.close()


# plot importances
models = ["Boosting", "Random Forest"]
for model in models:
    importances = np.asarray(model_importances[model])*100
    means = np.mean(importances, axis=0)
    stds = np.std(importances, axis=0)
    fig = plt.figure(figsize=(2.25, 6.4))
    ypos = [-i for i in range(len(names))]
    plt.barh(ypos, means, xerr=stds, align='center', color='tab:green', height=0.5, edgecolor='k')
    plt.xlabel("% of importance", fontsize=8)
    plt.yticks(ypos, names, fontsize=8)
    plt.title("Importance of variables (%s)" % (model) , fontsize=8)
    plt.gca().tick_params(axis='both', which='major', labelsize=8)
    plt.gca().tick_params(axis='both', which='minor', labelsize=8)
    #plt.savefig("Conditional_corr.png", bbox_inches='tight')
    plt.savefig(model+"_importance.svg", format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()
