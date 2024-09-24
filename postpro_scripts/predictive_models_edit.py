import sys
import re
import copy
import glob
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
import load_file_edit as load_file
import statis_edit as statis

### set data information (fname/field) and key
path = '/Users/sangwoopark/jhu_rockfish/2024_01_05_GEO/processed_files/'
dinfo_dkey = {'H1_NCP_sp_1rep_deep_chr1_score_table.gtab.gz':None}

### load annotation data
#print("Data reading start", file=sys.stderr)
print "Data reading start"

dkey_ID_value = {}
for fkey in dinfo_dkey:
    field_dkey = dinfo_dkey[fkey]

    if field_dkey == None:
        field_choices = None
    else:
        field_choices = field_dkey.keys()
    
    for fname in glob.glob(path + '*'):
        if not re.match(fkey, fname.rsplit('/')[-1]):
            continue
        #print ("loading %s" % (fname.rsplit('/')[-1]), file=sys.stderr)
        print "loading %s" % (fname.rsplit('/')[-1])

        field_ID_value = load_file.read_gtab(fname,
                                             mode='col',
                                             field_choices=field_choices)

        if field_dkey == None:
            field_dkey = {field:field for field in field_ID_value.keys()}

        for field, dkey in field_dkey.items():
            ID_value = field_ID_value[field]
            if dkey not in dkey_ID_value:
                dkey_ID_value[dkey] = {}
            dkey_ID_value[dkey].update(ID_value)

### Compute other sequence features
# AT content
dkey_ID_value['AT content'] = copy.deepcopy(dkey_ID_value['ATcontent'])
del dkey_ID_value['ATcontent']

# methylation density
dkey_ID_value['meCpG density'] = statis.get_fract_dict(dkey_ID_value['CNumber(CpG)'],
                                                       dkey_ID_value['meCNumber(CpG)'])
dkey_ID_value['meCHG density'] = statis.get_fract_dict(dkey_ID_value['CNumber(CHG)'],
                                                       dkey_ID_value['meCNumber(CHG)'])
dkey_ID_value['meCHH density'] = statis.get_fract_dict(dkey_ID_value['CNumber(CHH)'],
                                                       dkey_ID_value['meCNumber(CHH)'])

# mean poly GC length
ID_polyGC = {}
for ID, seq in dkey_ID_value['Sequence'].items():
    num_pos = statis.polynt_count(seq.upper(), nts='GC', pos=True)
    mean_len, count = 0.0, 0.0
    for num, pos in num_pos.items():
        mean_len += len(pos)*num
        count += len(pos)
    ID_polyGC[ID] = mean_len/count

dkey_ID_value['poly-G/C length'] = ID_polyGC

del dkey_ID_value['Sequence']
del ID_polyGC


### select target and feature sets
feature_names = ['AT content', 'poly-G/C length', 'meCpG density', 'meCHG density', 'meCHH density', 'H2AFZ', 'H2AK5ac', 'H2BK120ac', 'H2BK12ac', 'H2BK15ac', 'H2BK20ac', 'H2BK5ac', 'H3K14ac', 'H3K18ac', 'H3K23ac', 'H3K23me2', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K56ac', 'H3K79me1', 'H3K79me2', 'H3K9ac', 'H3K9me3', 'H4K20me1', 'H4K5ac', 'H4K8ac', 'H4K91ac']

target_name = 'H1_NCP_sp_8_1rep_deep'
ID_score = dkey_ID_value[target_name]


### randomly select IDs for analysis
sample_size = 10**5
random.seed(123)
IDs = random.sample(ID_score.keys(),
                    sample_size)


### get features and targets
features = [[] for i in range(len(IDs))]
for feature_name in feature_names:
    values = [dkey_ID_value[feature_name][ID] for ID in IDs]
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

del ID_score
del dkey_ID_value

#print("Data reading is done", file=sys.stderr)
print "Data reading is done"

### k-fold cross validation
k = 10
test_num = int(len(IDs)/10)

model_corrs = {}
model_importances = {}
for i in range(k):
    #print ("iteration %d/%d" % (i+1, k), file=sys.stderr)
    print "iteration %d/%d" % (i+1, k)
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
    #print ("Linear reg", corr)


    # Supported Vector machine
    svm = SVR(gamma='scale')
    svm.fit(train_features, train_targets)
    corr, _ = pearsonr(test_targets, svm.predict(test_features))
    if "SVR" not in model_corrs:
        model_corrs["SVR"] = []
    model_corrs["SVR"].append(corr)
    #print ("SVM reg", corr)


    # Boosting
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
    #print ("Boosting reg", corr)
    if "Boosting" not in model_importances:
        model_importances["Boosting"] = []
    model_importances["Boosting"].append(gb.feature_importances_)

    
    # Random Forest
    rf = RandomForestRegressor(min_samples_leaf=3,
                               n_estimators=100,
                               max_depth=20,
                               max_features=3)
    rf.fit(train_features, train_targets)
    corr, _ = pearsonr(test_targets, rf.predict(test_features))
    if "Random Forest" not in model_corrs:
        model_corrs["Random Forest"] = []
    model_corrs["Random Forest"].append(corr)
    #print ("Random Forest reg", corr)
    if "Random Forest" not in model_importances:
        model_importances["Random Forest"] = []
    model_importances["Random Forest"].append(rf.feature_importances_)


    # Neural Network (MLP)
    mlp = MLPRegressor(hidden_layer_sizes=(32,64,16),
                       max_iter=1000,
                       activation='tanh')
    mlp.fit(train_features, train_targets)
    corr, _ = pearsonr(test_targets, mlp.predict(test_features))
    if "Neural Network" not in model_corrs:
        model_corrs["Neural Network"] = []
    model_corrs["Neural Network"].append(corr)
    #print ("Neural Network(MLP) reg", corr)

# save the result
pickle.dump(model_corrs, open("model_corrs.pickle", "wb"))
pickle.dump(model_importances, open("model_importances.pickle", "wb"))


### plot the machine learning result
# print average correlation 
models = ["Linear reg", "SVR", "Boosting", "Random Forest", "Neural Network"]
for model in models:
    #print (model, np.mean(model_corrs[model]), file=sys.stderr)
    print model, np.mean(model_corrs[model])

# plot prediction vs experiment scatter plot
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

# plot correlations bar plots
color_list = ['tab:red', 'tab:blue', 'tab:orange', 'tab:green', 'tab:purple']
fig = plt.figure(figsize=(2.8, 2))
for i in range(len(models)):
    model = models[i]
    corrs = model_corrs[model]
    bp = plt.boxplot(corrs,
                     positions=[i],
                     patch_artist=True,
                     showfliers=False,
                     widths=[0.4])
    for patch in bp['boxes']:
        patch.set_facecolor(color_list[i])
    for element in ['whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color='black')
plt.xticks(range(len(models)),
           models,
           rotation=45,
           fontsize=8,
           ha="right",
           rotation_mode="anchor")
plt.gca().tick_params('y', labelsize=6)
plt.ylabel("Pearson correlation", fontsize=8)
plt.title("10-fold cross validation", fontsize=10)
plt.xlim([-0.5, len(models)-0.5])
plt.ylim([0.3, 0.8])
plt.savefig('Pearson_model.svg', format='svg', bbox_inches='tight')
#plt.show()
plt.close()

# plot importances bar plots
models = ["Boosting", "Random Forest"]
for model in models:
    importances = np.asarray(model_importances[model])*100
    means = np.mean(importances, axis=0)
    stds = np.std(importances, axis=0)
    fig = plt.figure(figsize=(2.25, 6.4))
    ypos = [-i for i in range(len(names))]
    plt.barh(ypos,
             means,
             xerr=stds,
             align='center',
             color='tab:green',
             height=0.5,
             edgecolor='k')
    plt.xlabel("% of importance", fontsize=8)
    plt.yticks(ypos, names, fontsize=8)
    plt.title("Importance of variables (%s)" % (model) , fontsize=8)
    plt.gca().tick_params(axis='both', which='major', labelsize=8)
    plt.gca().tick_params(axis='both', which='minor', labelsize=8)
    #plt.savefig("Conditional_corr.png", bbox_inches='tight')
    plt.savefig(model+"_importance.svg", format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()
