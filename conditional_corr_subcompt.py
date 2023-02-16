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


### parameters
path = "/home/spark159/../../storage/"

#exp_list = [('H1', 'NCP', 'sp', 8),
#            ('H1', 'NCP', 'spd', 6),
#            ('H1', 'NCP', 'CoH', 5),
#            ('H1', 'NCP', 'PEG', 6),
#            ('H1', 'NCP', 'Ca', 5),
#            ('H1', 'NCP', 'Mg', 5),
#            ('H1', 'NCP', 'HP1a', 3),
#            ('H1', 'NCP', 'HP1bSUV', 4),
#            ('H1', 'NCP', 'LKH', 3),
#            ('H1', 'NCP', 'Ki67', 4),
#            ('H1', 'NCP', 'FUS', 5)]

exp_list = [('H1', 'NCP', 'sp', 8),
            ('H1', 'NCP', 'spd', 6),
            ('H1', 'NCP', 'CoH', 5),
            ('H1', 'NCP', 'PEG', 6),
            ('H1', 'NCP', 'Ca', 5),
            ('H1', 'NCP', 'HP1a', 3),
            ('H1', 'NCP', 'LKH', 3),
            ('H1', 'NCP', 'Ki67', 4)]

exp_list = [('H1', 'NCP', 'sp', 8),
            ('H1', 'NCP', 'HP1a', 3),
            ('H1', 'NCP', 'LKH', 3),
            ('H1', 'NCP', 'Ki67', 4)]

# set domain names
domain_names = {'NSpeckle':['SON'],
                 'Trnx':['POLR2A', 'POLR2AphosphoS5', 'H3K9ac', 'H3K4me3', 'H3K27ac'],
                 'Polycomb':['CBX8', 'EZH2', 'RNF2', 'SUZ12'],
                 'Hetero':['H3K9me3', 'CBX5'],
                 'Nucleolus':['Nucleolar'],
                 'Lamin':['LaminB1'],
                 'Compartment':['eigen'],
                 'other':['ATcontent']}

domains = ['NSpeckle', 'Trnx', 'Polycomb', 'Hetero', 'Nucleolus', 'Lamin']

names = []
for domain in domains:
    names += domain_names[domain]

# binsize
bin_size = 10000

# other parameters
#dtype = 'zscore'
dtype = 'score'


# get enrichment zscore for each data
exp_name_corr = {}
for exp in exp_list:
    cell, sample, agent, tnum = exp
    print exp
    
    # load data
    fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', 'anot']) + '.txt'
    score_field = "%s-%s-%s-%d" % (cell, sample, agent, tnum)
    
    name_ID_value = load_file.read_tabular_file(fname, mode='col')
    ID_chr = name_ID_value['Chromosome']
    ID_score = name_ID_value[score_field]

    IDs = sorted(ID_score.keys())
    
    # binning the features and get state
    corr_list = []
    ID_state = {}
    for name in names:
        values = [name_ID_value[name][ID] for ID in IDs] 
        scores = [ID_score[ID] for ID in IDs]
        corr = statis.get_spearman_corr(values, scores)
        print name, corr
        corr_list.append(corr)

        if exp not exp_name_corr:
            exp_name_corr[exp] = []
        if name not in exp_name_corr[exp]:
            exp_name_corr[exp][name] = corr

        min_value = min(values)
        max_value = max(values)
        for ID, value in zip(IDs, values):
            #if name.startswith('H'):
            #    re_value = round(float(value-min_value)/max_value,1)
            #else:
            #    re_value = round(float(value-min_value)/max_value,2)
            if ID not in ID_state:
                ID_state[ID] = []
            #ID_state[ID].append(re_value)
            ID_state[ID].append(value)
        del values
        del scores

    del ID_chr
    del name_ID_value

    print "Data reading is done"

    # plot spearman correlation
    yset1, dataset1 = [], []
    yset2, dataset2 = [], []
    for i in range(len(names)):
        corr = corr_list[i]
        if corr < 0:
            yset1.append(-i)
            dataset1.append(corr)
        else:
            yset2.append(-i)
            dataset2.append(corr)

    fig = plt.figure(figsize=(2.25, 6.4))
    plt.barh(yset1, dataset1, align='center', color='tab:red', height=0.5, edgecolor='k')
    plt.barh(yset2, dataset2, align='center', color='tab:blue', height=0.5, edgecolor='k')
    plt.axvline(x=0, linestyle='--', color='k', linewidth=1)
    plt.xlabel("Spearman correlation", fontsize=8)
    #plt.xlabel("Averaged correlation", fontsize=8, rotation=180) # for flip version
    plt.yticks([-i for i in range(len(names))], names, fontsize=8)
    #plt.title("Conditional Correlation", fontsize=8)
    plt.gca().tick_params(axis='both', which='major', labelsize=8)
    plt.gca().tick_params(axis='both', which='minor', labelsize=8)
    #plt.xticks(rotation=-90) # for flip version
    #plt.yticks(rotation=-20, ha="right", va='center', rotation_mode="anchor") # for flip version
    #plt.savefig("Conditional_corr.png", bbox_inches='tight')
    plt.savefig("Spearman_corr_%.svg" % (agent), format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

    fig = plt.figure()
    
    
    
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
            X = [ID_state[ID][i] for ID in IDs]
            Y = [ID_score[ID] for ID in IDs]

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
    color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'lime', 'salmon']*5

    #fig = plt.figure(figsize=(3,4))
    fig = plt.figure(figsize=(2.25, 6.4))
    for i in range(len(names)):
        corrs = np.asarray(corrs_list[i])
        weights = np.asarray(weights_list[i])
        frac_weights = weights / float(sum(weights))
        rgb_colors = np.zeros((len(corrs),4))
        rgb_colors[:,:3] = colors.to_rgba(color_list[i])[:3]
        rgb_colors[:,3] = 0.15 * weights / float(max(weights))
        order = np.argsort(weights)
        plt.scatter(corrs[order], [-i]*len(corrs), s=5000*frac_weights[order], color=rgb_colors[order])
        plt.annotate('x', (cdcorr_list[i], -i), ha='center', va='center')
    plt.axvline(x=0, linestyle='--', color='k')
    plt.xlim([-0.25, 0.25])
    #plt.xlabel("Spearman correlation", fontsize=8)
    plt.xlabel("Spearman correlation", fontsize=8, rotation=180) # for flip version
    plt.yticks([-i for i in range(len(names))], names, fontsize=8)
    #plt.title("Data Stratification", fontsize=8)
    plt.gca().tick_params(axis='both', which='major', labelsize=8)
    plt.gca().tick_params(axis='both', which='minor', labelsize=8)
    plt.xticks(rotation=-90) # for flip version
    plt.yticks(rotation=-20, ha="right", va='center', rotation_mode="anchor") # for flip version
    plt.savefig("Data_strat.png", bbox_inches='tight', dpi=500)
    #plt.savefig("Data_strat.svg", format='svg', bbox_inches='tight')
    #plt.show()
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

    #fig = plt.figure(figsize=(3,4))
    fig = plt.figure(figsize=(2.25, 6.4))
    plt.barh(yset1, dataset1, align='center', color='tab:red', height=0.5, edgecolor='k')
    plt.barh(yset2, dataset2, align='center', color='tab:blue', height=0.5, edgecolor='k')
    plt.axvline(x=0, linestyle='--', color='k', linewidth=1)
    #plt.xlabel("Averaged correlation", fontsize=8)
    plt.xlabel("Averaged correlation", fontsize=8, rotation=180) # for flip version
    plt.yticks([-i for i in range(len(names))], names, fontsize=8)
    #plt.title("Conditional Correlation", fontsize=8)
    plt.gca().tick_params(axis='both', which='major', labelsize=8)
    plt.gca().tick_params(axis='both', which='minor', labelsize=8)
    plt.xticks(rotation=-90) # for flip version
    plt.yticks(rotation=-20, ha="right", va='center', rotation_mode="anchor") # for flip version
    #plt.savefig("Conditional_corr.png", bbox_inches='tight')
    plt.savefig("Conditional_corr.svg", format='svg', bbox_inches='tight')
    #plt.show()
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
