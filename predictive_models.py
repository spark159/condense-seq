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


# load annotation data
print("Data reading start", file=sys.stderr)

path = ''
ID_chr, ID_pos, name_ID_value = read_anot_file(path+"H1_NCP_sp_chr1_extended_anot.cn")
ID_score = name_ID_value['work/2021_06_07_H1_sp_detail/H1-NCP-sp-8']
name_ID_value['AT content'] = name_ID_value['ATcontent']
name_ID_value['meCpG density'] = get_density(name_ID_value['CNumber(CpG)'],
                                             name_ID_value['meCNumber(CpG)'])
name_ID_value['meCHG density'] = get_density(name_ID_value['CNumber(CHG)'],
                                             name_ID_value['meCNumber(CHG)'])
name_ID_value['meCHH density'] = get_density(name_ID_value['CNumber(CHH)'],
                                             name_ID_value['meCNumber(CHH)'])

ID_seq = name_ID_value['Sequence']
ID_polyGC = {}
for ID in ID_seq:
    seq = ID_seq[ID]
    num_pos = poly_score(seq.upper(), nts='GC', pos=True)
    mean_len, count = 0.0, 0.0
    for num, pos in num_pos.items():
        mean_len += len(pos)*num
        count += len(pos)
    ID_polyGC[ID] = mean_len/count

name_ID_value['poly-G/C length'] = ID_polyGC

del ID_seq
del ID_polyGC


# select feature sets
names = ['AT content', 'poly-G/C length', 'meCpG density', 'meCHG density', 'meCHH density', 'H2AFZ', 'H2AK5ac', 'H2BK120ac', 'H2BK12ac', 'H2BK15ac', 'H2BK20ac', 'H2BK5ac', 'H3K14ac', 'H3K18ac', 'H3K23ac', 'H3K23me2', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K56ac', 'H3K79me1', 'H3K79me2', 'H3K9ac', 'H3K9me3', 'H4K20me1', 'H4K5ac', 'H4K8ac', 'H4K91ac']


# shuffle the IDs
IDs = list(ID_pos.keys())
random.seed(123)
random.shuffle(IDs)
IDs = IDs[:100000]


# get features and targets
features = [[] for i in range(len(IDs))]
for name in names:
    values = [name_ID_value[name][ID] for ID in IDs]
    min_value = min(values)
    max_value = max(values)
    for i in range(len(IDs)):
        value = values[i]
        re_value = float(value-min_value)/max_value
        features[i].append(re_value)
    del values
    del min_value
    del max_value
    del re_value

targets = [ID_score[ID] for ID in IDs]

del ID_pos
del ID_chr
del ID_score
del name_ID_value

print("Data reading is done", file=sys.stderr)

# k-fold cross validation
k = 10
test_num = int(len(IDs)/10)

model_corrs = {}
model_importances = {}
for i in range(k):
    print ("iteration %d/%d" % (i+1, k), file=sys.stderr)
    st, ed = test_num*i, test_num*(i+1)

    test_features, test_targets = features[st:ed], targets[st:ed]
    train_features = features[:st] + features[ed:] 
    train_targets = targets[:st] + targets[ed:]
    

    # Linear regression
    alpha = 0.5
    reg = linear_model.Ridge(alpha=alpha)
    reg.fit(train_features, train_targets)
    corr, _ = pearsonr(test_targets, reg.predict(test_features))
    if "Linear reg" not in model_corrs:
        model_corrs["Linear reg"] = []
    model_corrs["Linear reg"].append(corr)
    #print "Linear reg", corr


    # Supported Vector machine
    svm = SVR(gamma='scale')
    svm.fit(train_features, train_targets)
    corr, _ = pearsonr(test_targets, svm.predict(test_features))
    if "SVR" not in model_corrs:
        model_corrs["SVR"] = []
    model_corrs["SVR"].append(corr)
    #print "SVM reg", corr


    # Boosting
    #gb = GradientBoostingRegressor()
    gb = GradientBoostingRegressor(min_samples_leaf=3,
                                   n_estimators=198,
                                   max_depth=20,
                                   max_features=2,
                                   min_samples_split=1000)
    gb.fit(train_features, train_targets)
    corr, _ = pearsonr(test_targets, gb.predict(test_features))
    if "Boosting" not in model_corrs:
        model_corrs["Boosting"] = []
    model_corrs["Boosting"].append(corr)
    #print "Boosting reg", corr
    if "Boosting" not in model_importances:
        model_importances["Boosting"] = []
    model_importances["Boosting"].append(gb.feature_importances_)

    
    # Random Forest
    #rf = RandomForestRegressor()
    rf = RandomForestRegressor(min_samples_leaf=3,
                               n_estimators=100,
                               max_depth=20,
                               max_features=3)
    rf.fit(train_features, train_targets)
    corr, _ = pearsonr(test_targets, rf.predict(test_features))
    if "Random Forest" not in model_corrs:
        model_corrs["Random Forest"] = []
    model_corrs["Random Forest"].append(corr)
    #print "Random Forest reg", corr
    if "Random Forest" not in model_importances:
        model_importances["Random Forest"] = []
    model_importances["Random Forest"].append(rf.feature_importances_)


    # Neural Network (MLP)
    #mlp = MLPRegressor()
    mlp = MLPRegressor(hidden_layer_sizes=(32,64,16),
                       max_iter=1000,
                       activation='tanh')
    mlp.fit(train_features, train_targets)
    corr, _ = pearsonr(test_targets, mlp.predict(test_features))
    if "Neural Network" not in model_corrs:
        model_corrs["Neural Network"] = []
    model_corrs["Neural Network"].append(corr)
    #print "Neural Network(MLP) reg", corr


f5 = open(fnames[i] + "_pos_probs_" + str(i) + ".pickle", "rb")

# print average correlation 
models = ["Linear reg", "SVR", "Boosting", "Random Forest", "Neural Network"]
for model in models:
    print (model, np.mean(model_corrs[model]), file=sys.stderr)

# plot prediction vs experiment
for model, m in zip(models, [reg, svm, gb, rf, mlp]):
    X = list(test_targets)
    Y = list(m.predict(test_features))
    pickle.dump([X, Y], open("ExpVSPred_%s.pickle" % (model), "wb"))
    fig = plt.figure(figsize=(2.8, 2))
    plt.plot(X, Y, 'k.', ms=1)
    plt.xlabel("Experiment", fontsize=8)
    plt.ylabel("Prediction", fontsize=8)
    plt.gca().tick_params(axis='both', which='major', labelsize=6)
    plt.gca().tick_params(axis='both', which='minor', labelsize=6)
    #minvalue, maxvalue = min(X+Y), max(X+Y)
    #plt.xlim([minvalue, maxvalue])
    #plt.ylim([minvalue, maxvalue])
    plt.title("Model prediction (%s)" % (model), fontsize=10)
    plt.savefig('ExpVSPre_%s.svg' % (model), format='svg', bbox_inches='tight')
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
plt.ylim([0.3, 0.8])
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


# save result
pickle.dump(model_corrs, open("model_corrs.pickle", "wb"))
pickle.dump(model_importances, open("model_importances.pickle", "wb"))
