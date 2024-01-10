import glob
import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import copy
import re
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy
import numpy as np
import sklearn.manifold
from sklearn.decomposition import PCA
from sklearn import preprocessing
from sklearn.neighbors import LocalOutlierFactor
#from pymol_graphics import Molecules
import random
from scipy.optimize import curve_fit
from sklearn import linear_model
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import pickle
from scipy.stats import norm

def rescale (value, old_st, old_ed, new_st, new_ed):
    new_value = new_st + (new_ed - new_st)*float(value-old_st)/(old_ed-old_st)
    return new_value


def GC_content(seq):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    return (num/float(len(seq)))*100


# Oncohistone libraray NGS informatino
def read_index (fname):
    index_name, name_index = {}, {}
    name_titr = {}
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        cols = line.strip().split('\t')
        index = cols[-1]
        name = cols[0]
        titr = int(cols[3])
        assert index not in index_name
        assert name not in name_index
        index_name[index] = name
        name_index[name] = index
        name_titr[name] = titr
    return index_name, name_index, name_titr
index_name, name_index, name_titr = read_index("Oncolib_NGS_information.csv")

# Oncohistone library barcode information
def read_table (fname):
    BC_BCnum, BCnum_BC = {}, {}
    BC_histone, histone_BC = {}, {}
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        cols = line.strip().split('\t')
        BC = cols[-2]
        histone = cols[0].strip()
        #BCnum = int(cols[-1].split('_')[1][2:])
        if BC == '-':
            continue
        assert BC not in BC_histone
        assert histone not in histone_BC
        #assert BCnum not in BCnum_BC
        #assert BC not in BC_BCnum
        BC_histone[BC] = histone
        histone_BC[histone] = BC
        #BC_BCnum[BC] = BCnum
        #BCnum_BC[BCnum] = BC
    return BC_histone, histone_BC
BC_histone, histone_BC = read_table("OncohistoneTable.csv")

# read amino acid information
def read_aa (fname):
    aa_info = {}
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            field = cols[3:]
            First = False
            continue
        name, triple, single = cols[:3]
        properties = cols[3:]

        assert single not in aa_info
        aa_info[single] = {"name":name, "triple":triple}

        for i in range(len(field)):
            key = field[i]
            try:
                value = float(properties[i])
            except:
                value = properties[i]
            aa_info[single][key] = value
    return aa_info
aa_info = read_aa("amino_acid.txt")

# read the published data
def read_old (fname):
    histone_mean = {}
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        if not line:
            continue
        
        #print (line)
        cols = line.strip().split(',')
        rawname, mean, std = cols
        mean = float(mean)
        #mean, std = float(mean), float(std)

        if rawname.startswith('H2A.V'):
            histone = 'H2A.V ' + rawname[5:]
        elif rawname.startswith('H2A.Z'):
            histone = 'H2A.Z ' + rawname[5:]
        elif rawname.startswith('H2A'):
            histone = 'H2A ' + rawname[3:]
        elif rawname.startswith('H2B'):
            histone = 'H2B ' + rawname[3:]
        elif rawname.startswith('H3.1'):
            histone = 'H3.1 ' + rawname[4:]
        elif rawname.startswith('H4'):
            histone = 'H4 ' + rawname[2:]
        elif rawname.startswith('H3.3'):
            histone = 'H3.3 ' + rawname[4:]
        else:
            histone = rawname
            
        histone = histone.strip()
        histone_mean[histone] = mean
    return histone_mean

histone_ACF = read_old("ACF.csv")
histone_Nap1 = read_old("Nap1.csv")
histone_freeDNA = read_old("Onco_freeDNA.csv")


# read sort file
sort_fname = "Sp-Spd-CoHex-PEG-HP1a-Oncolib_S1_L001_R1_001.sort"
validity_type_count = {}
name_count = {}
agent_num_histone_count = {}
for line in open(sort_fname):
    line = line.strip()
    if line.startswith('@'):
        validity, windows = line.split('::')[1].split(':')
        if validity not in validity_type_count:
            validity_type_count[validity] = {}
        if windows not in validity_type_count[validity]:
            validity_type_count[validity][windows] = 0
        validity_type_count[validity][windows] +=1
        if validity == 'invalid':
            continue
        BC, index = windows[1:-1].split('][')
        try:
            name = index_name[index]
            histone = BC_histone[BC]
        except:
            continue

        if name not in name_count:
            name_count[name] = 0
        name_count[name] +=1
        
        agent = name[:-1]
        num = int(name[-1])
        if agent not in agent_num_histone_count:
            agent_num_histone_count[agent] = {}
        if num not in agent_num_histone_count[agent]:
            agent_num_histone_count[agent][num] = {}
        if histone not in agent_num_histone_count[agent][num]:
            agent_num_histone_count[agent][num][histone] = 0
        agent_num_histone_count[agent][num][histone] +=1


# check data type
types, counts = [], []
for validity in validity_type_count:
    if validity == 'valid':
        type = 'valid'
        count = sum(validity_type_count[validity].values())
        types.append(type)
        counts.append(count)
    else:
        for type, count in validity_type_count[validity].items():
            types.append('invalid:' + type)
            counts.append(count)

fig = plt.figure()
plt.pie(counts, labels=types, shadow=True, startangle=90, autopct='%1.1f%%')
#plt.show()
plt.close()


# check count by sample name
X, Y = [], []
for name, count in name_count.items():
    X.append(name)
    Y.append(count)

fig = plt.figure()
plt.bar(X, Y)
plt.xticks(rotation=70)
#plt.show()
plt.close()


# check count by BC (input sample only)
all_histones = sorted(histone_BC.keys())
fig = plt.figure()
for agent in agent_num_histone_count.keys():
    histone_count = agent_num_histone_count[agent][0]
    X, Y = [], []
    for histone in all_histones:
        X.append(histone)
        Y.append(histone_count[histone])
    plt.plot(range(len(Y)), Y, '.-', alpha=0.8, label=agent)
plt.xlabel('Oncohistone IDs')
plt.ylabel('Read counts')
plt.legend()
#plt.show()
plt.close()


# parsing the mutation information
all_good_histones = list(set(all_histones) - set(['H3.1 K4M', 'H3.1 E97A', 'H4 G42V']))
cate_histones = {}
histone_minfo = {}
for histone in all_good_histones:
    if '-' in histone:
        tag, subunit = histone.split('-')
        if tag == 'Biotin':
            pos = 115-1
        else:
            assert tag in ['HA', 'FLAG']
            pos = 0
        histone_minfo[histone] = {}
        histone_minfo[histone][subunit] = {}
        histone_minfo[histone][subunit]["tag"] = {}
        histone_minfo[histone][subunit]["tag"][pos] = tag

        cate = 'tag'
        if cate not in cate_histones:
            cate_histones[cate] = []
        cate_histones[cate].append(histone)
        continue
    
    cols = histone.split(' ')
    if len(cols) != 2:
        assert len(cols) == 1
        histone_minfo[histone] = {}
        histone_minfo[histone][histone] = None

        cate = histone
        assert cate not in cate_histones
        cate_histones[cate] = [histone]
        continue
    
    subunit, mutations = cols

    if subunit.startswith('WT'):
        histone_minfo[histone] = None

        cate = 'WT'
        if cate not in cate_histones:
            cate_histones[cate] = []
        cate_histones[cate].append(histone)
        continue
    
    if mutations.startswith('DNA'):
        histone_minfo[histone] = 'freeDNA'

        cate = 'freeDNA'
        if cate not in cate_histones:
            cate_histones[cate] = []
        cate_histones[cate].append(histone)
        continue

    finds_list = [re.findall('\d+|\D+', mut) for mut in mutations.split('/')]
    for finds in finds_list: 
        if len(finds) == 3:
            pre_aa, pos, post_aa = finds
            pos = int(pos)
        elif len(finds) == 2:
            if finds[0].isalpha():
                pre_aa, pos = finds
                pos = int(pos)
                post_aa = finds_list[-1][-1]
            else:
                assert finds[1].isalpha()
                pos, post_aa = finds
                pos = int(pos)
                pre_aa = finds_list[0][0]
        else:
            assert len(finds) == 1
            pos = int(pos)
            pre_aa = finds_list[0][0]
            post_aa = finds_list[-1][-1]

        if post_aa not in aa_info.keys():
            type = 'tag'
        else:
            type = 'mut'

        if histone not in histone_minfo:
            histone_minfo[histone] = {}
        if subunit not in histone_minfo[histone]:
            histone_minfo[histone][subunit] = {}
        if type not in histone_minfo[histone][subunit]:
            histone_minfo[histone][subunit][type] = {}

        if type == 'tag':
            histone_minfo[histone][subunit][type][pos] = post_aa
        elif type == 'mut':
            histone_minfo[histone][subunit][type][pos] = (pre_aa, post_aa)

    if type == 'tag':
        cate = 'tag'
    else:
        assert type == 'mut'
        if subunit in ['H2A', 'H2B', 'H3.1', 'H4']:
            cate = 'WT+mut'
        else:
            cate = subunit + '+mut'

    if cate not in cate_histones:
        cate_histones[cate] = []
    cate_histones[cate].append(histone)

# categorize histones
histone_cate = {}
for cate, histones in cate_histones.items():
    for histone in histones:
        histone_cate[histone] = cate

# acidic patch amino acids
AP_info = {"H2A":["E56", "E61", "E64", "D90", "E91", "E92"], "H2B":["E105", "E113"]}
AP_sites = ["H2A E56", "H2A E61", "H2A E64", "H2A D90", "H2A E91", "H2A E92", "H2B E105", "H2B E113"]
APmtype_histones = {"Positive":["H2A E56K", "H2A E61K", "H2A E64K", "H2A D90K", "H2A E91K", "H2A E92K", "H2B E105K", "H2B E113K"]}

APmutants = []
for histone in cate_histones['WT+mut']:
    for subunit in histone_minfo[histone]:
        for pos in histone_minfo[histone][subunit]['mut']:
            pre_aa, post_aa = histone_minfo[histone][subunit]['mut'][pos]
            name = "%s %s%s" % (subunit, pre_aa, pos)
            #print (name)
            if name in AP_sites:
                APmutants.append(histone)
                break

# make absolute titration curve
def read_data (fname):
    conc_list = []
    mean_list = []
    std_list = []
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        if not line:
            continue
        cols = line.strip().split()
        conc, mean, std = float(cols[0]), float(cols[1]), float(cols[2])
        conc_list.append(conc)
        mean_list.append(mean)
        std_list.append(std)
    return conc_list, mean_list, std_list

agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
agent_fullname = {'sp':'Spermine(4+)', 'spd':'Spermidine(3+)', 'CoH':'Cobalt Hexammine(3+)', 'PEG':'PEG 8000', 'HP1a':'HP1 $\\alpha$'}
agent_filename = {'sp':"Oncolib_spermine.csv", 'spd':"Oncolib_spermidine.csv", 'CoH':"Oncolib_CoHex.csv", 'PEG':"Oncolib_PEG.csv", 'HP1a':"Oncolib_HP1a.csv"}

agent_num_histone_fraction = {}
for agent in agent_list:
    for num in agent_num_histone_count[agent]:
        total = sum(agent_num_histone_count[agent][num].values())
        histone_fraction = {}
        for histone in agent_num_histone_count[agent][num]:
            count = agent_num_histone_count[agent][num][histone]
            fraction = float(count) / total
            histone_fraction[histone] = fraction
        if agent not in agent_num_histone_fraction:
            agent_num_histone_fraction[agent] = {}
        agent_num_histone_fraction[agent][num] = copy.deepcopy(histone_fraction)

agent_histone_titration = {}
for agent in agent_list:
    full_conc_list, full_mean_list, _ = read_data(agent_filename[agent])
    #print (full_conc_list)
    num_list = sorted(agent_num_histone_fraction[agent].keys())
    for num in num_list:
        if num == 0:
            conc = 0.0
            mean = 1.0
        else:
            titr = name_titr[agent + str(num)] - 2
            if agent in ['PEG', 'HP1a']:
                assert full_conc_list[0] == 0
                titr +=1
            conc = full_conc_list[titr]
            mean = full_mean_list[titr]
            #print (titr)
        for histone, fraction in agent_num_histone_fraction[agent][num].items():
            survival = mean * fraction
            #print (survival)
            if agent not in agent_histone_titration:
                agent_histone_titration[agent] = {}
            if histone not in agent_histone_titration[agent]:
                agent_histone_titration[agent][histone] = {}
            if 'conc' not in agent_histone_titration[agent][histone]:
                agent_histone_titration[agent][histone]['conc'] = []
            if 'survival' not in agent_histone_titration[agent][histone]:
                agent_histone_titration[agent][histone]['survival'] = []
            agent_histone_titration[agent][histone]['conc'].append(conc)
            agent_histone_titration[agent][histone]['survival'].append(survival)

for agent in agent_list:
    for histone in list(agent_histone_titration[agent].keys()):
        input = float(agent_histone_titration[agent][histone]['survival'][0])
        for i in range(len(agent_histone_titration[agent][histone]['survival'])):
            agent_histone_titration[agent][histone]['survival'][i] /= input

# define the metrics
# fitting with sigmoidal curve and get "C-half"
# get condensabiltiy "Score" by calculating <-log(survival)>
def sigmoid(x, L ,x0, k):
    y = L / (1 + np.exp(k*(x-x0)))
    return (y)

agent_histone_Chalf = {}
agent_histone_fitting = {}
agent_histone_score = {}
for agent in agent_list:
    #fig = plt.figure()
    fig = plt.figure(figsize=(2, 1.4))
    for histone in all_good_histones:
        #if histone not in ['WT #1', 'WT #2', 'Cuttable DNA', 'Uncuttable DNA']:
        #    continue
        #if not histone.endswith('ub'):
        #    continue

        #fig = plt.figure()
        
        X = agent_histone_titration[agent][histone]['conc']
        Y = agent_histone_titration[agent][histone]['survival']
        if agent == 'sp':
            X = X[:6] + X[7:]
            Y = Y[:6] + Y[7:]
            
        #X, Y = X[1:], Y[1:]

        p0 = [max(Y), np.median(X), 1]
        bounds = ([0.0, 0.0, 0.0], [max(Y)+max(Y)*0.1, np.inf, np.inf])

        popt, pcov = curve_fit(sigmoid, X, Y, p0, bounds = bounds,  method='dogbox')
        residuals = np.asarray(Y)- sigmoid(X, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((np.asarray(Y)-np.mean(Y))**2)
        r_squared = 1 - (ss_res / ss_tot)
        pred_X = np.linspace(min(X[1:]), max(X[1:]), 1000)
        pred_Y = sigmoid(pred_X, *popt)

        c = plt.plot(X[1:], Y[1:], '.', alpha=0.3)
        plt.plot(pred_X, pred_Y, '-', color=c[0].get_color(), alpha=0.3, label=histone)
        #plt.axvline(x=popt[1], linestyle='--', color=c[0].get_color(), alpha=0.5)

        #plt.title("%s %s" % (agent, histone))
        #plt.xlabel("Concentration")
        #plt.ylabel("Soluble fraction")
        #if agent in ['HP1a']:
        #    plt.xscale('log', base=2)
        #elif agent in ['sp', 'spd', 'CoH']:
        #    plt.xscale('log', base=10)
        #plt.show()
        #plt.close()

        if agent not in agent_histone_Chalf:
            agent_histone_Chalf[agent] = {}
        if histone not in agent_histone_Chalf[agent]:
            agent_histone_Chalf[agent][histone] = popt[1]

        if agent not in agent_histone_fitting:
            agent_histone_fitting[agent] = {}
        if histone not in agent_histone_fitting[agent]:
            agent_histone_fitting[agent][histone] = {}
        agent_histone_fitting[agent][histone]['para'] = popt
        agent_histone_fitting[agent][histone]['r-square'] = r_squared

        score = np.mean(-np.log2(np.asarray(Y[1:])))

        if agent not in agent_histone_score:
            agent_histone_score[agent] = {}
        agent_histone_score[agent][histone] = score

        #print (popt)


    plt.gca().tick_params(axis='both', which='major', labelsize=5)
    plt.gca().tick_params(axis='both', which='minor', labelsize=5)

    #plt.xlabel("Spermine concentration (mM)", fontsize=8)
    #plt.title("Oncohistone library", fontsize=8)

        
    #plt.title(agent_fullname[agent])
    plt.xlabel("Concentration")
    plt.ylabel("Soluble fraction", fontsize=8)

    #if agent in ['HP1a']:
    #    plt.xscale('log', base=2)
    #elif agent in ['sp', 'spd', 'CoH']:
    #    plt.xscale('log', base=10)

    if agent in ['HP1a']:
        plt.xscale('log', basex=2)
    elif agent in ['sp', 'spd', 'CoH']:
        plt.xscale('log', basex=10)
    #plt.legend()
    #plt.show()
    #plt.savefig(agent+'.png')
    #plt.savefig(agent+'.svg', format='svg', bbox_inches='tight')
    plt.close()    

#sys.exit(1)
    
# check r-square of fitting
for agent in agent_list:
    data = []
    for histone in agent_histone_fitting[agent]:
        r_squared = agent_histone_fitting[agent][histone]['r-square']
        data.append(r_squared)

    fig = plt.figure()
    p = plt.hist(data, bins=50, range=(0.5,1))
    plt.xlim([0.5, 1])
    #plt.title(agent)
    plt.text(np.mean(plt.gca().get_xlim()), np.mean(plt.gca().get_ylim()), '$R^{2}=$%.2f' % (np.mean(data)), fontsize=20, horizontalalignment='left',  verticalalignment='center')
    plt.xlabel("R-squared")
    plt.ylabel("Counts")
    #plt.savefig(agent+'_r.png')
    #plt.show()
    plt.close()

# save the scores
#f = open("Oncolib_scores.txt", 'w')
#s = ['name'] + [agent for agent in agent_list]
#print >> f, '\t'.join(s)
#for histone in all_good_histones:
#    s = [histone]
#    for agent in agent_list:
#        score = agent_histone_score[agent][histone]
#        s.append(str(score))
#    print >> f, '\t'.join(s)
#f.close()

#sys.exit(1)

# get difference w.r.t. wild type
agent_histone_dChalf = {}
agent_histone_dscore = {}
all_good_histones_exceptWT = list(set(all_good_histones) - set(cate_histones['WT']))
for agent in agent_list:
    WT_Chalf, WT_score = [], []
    for histone in cate_histones['WT']:
        WT_Chalf.append(agent_histone_Chalf[agent][histone])
        WT_score.append(agent_histone_score[agent][histone])
    WT_Chalf = np.mean(WT_Chalf)
    WT_score = np.mean(WT_score)

    for histone in all_good_histones_exceptWT:
        dChalf = agent_histone_Chalf[agent][histone] - WT_Chalf
        dscore = agent_histone_score[agent][histone] - WT_score

        if agent not in agent_histone_dChalf:
            agent_histone_dChalf[agent] = {}
        agent_histone_dChalf[agent][histone] = dChalf

        if agent not in agent_histone_dscore:
            agent_histone_dscore[agent] = {}
        agent_histone_dscore[agent][histone] = dscore


# find out outliers based on two metircs
histone_list_list = [all_good_histones_exceptWT, cate_histones['WT+mut']]
agent_outliers_list = []
for k in range(len(histone_list_list)):
    histone_list = histone_list_list[k]
    agent_outliers = {}
    for agent in agent_list:
        data_list = []
        for histone in histone_list:
            dChalf = agent_histone_dChalf[agent][histone]
            dscore = agent_histone_dscore[agent][histone]
            data_list.append([dChalf, dscore])

        clf = LocalOutlierFactor()
        outcheck = clf.fit_predict(data_list)

        for i in range(len(outcheck)):
            if outcheck[i] < 0:
                outlier = histone_list[i]
                if agent not in agent_outliers:
                    agent_outliers[agent] = []
                agent_outliers[agent].append(outlier)


        #fig = plt.figure()
        fig = plt.figure(figsize=(2, 1.4))
        for histone in histone_list:
            dChalf = agent_histone_dChalf[agent][histone]
            dscore = agent_histone_dscore[agent][histone]
            if histone not in agent_outliers[agent]:
                plt.plot(dChalf, dscore, 'k.', markersize=2.5)
            else:
                plt.plot(dChalf, dscore, 'k.', markersize=2.5)
                #plt.annotate(histone, (dChalf, dscore))

        plt.axvline(x=0, linestyle='--', color='k', alpha=0.5)
        plt.axhline(y=0, linestyle='--', color='k', alpha=0.5)
        plt.gca().tick_params(axis='both', which='major', labelsize=5)
        plt.gca().tick_params(axis='both', which='minor', labelsize=5)
        plt.title(agent_fullname[agent], fontsize=8)
        plt.xlabel("$\Delta$ C-half", fontsize=8)
        plt.ylabel("$\Delta$ Score", fontsize=8)
        #plt.savefig(agent + "_ChalfVSScore.svg", format='svg', bbox_inches='tight')
        #plt.title(agent_fullname[agent])
        #plt.xlabel("$\Delta$ C-half")
        #plt.ylabel("$\Delta$ Score")
        #plt.savefig(str(k) + '_' + agent + "_ChalfVSScore.png")
        #plt.show()
        plt.close()


    agent_outliers_list.append(copy.deepcopy(agent_outliers))
all_agent_outliers, WTmut_agent_outliers = agent_outliers_list

#sys.exit(1)

# get pseudo p-value based on WT distributions
agent_histone_pvalue = {}
for agent in agent_list:
    WTscores = [agent_histone_score[agent][histone] for histone in cate_histones['WT']]
    WTmean = np.mean(WTscores)
    WTstd = np.std(WTscores)
    for histone in all_good_histones_exceptWT:
        score = agent_histone_score[agent][histone]
        re_score = float(score-WTmean)/WTstd
        if re_score >=0:
            pvalue = 1 - norm.cdf(re_score)
        else:
            pvalue = norm.cdf(re_score)
        if agent not in agent_histone_pvalue:
            agent_histone_pvalue[agent] = {}
        agent_histone_pvalue[agent][histone] = pvalue

# dscore VS p-value scatter plot (WT single mutant only)
# find outliers with p-value < 0.01
agent_outliers = {}
for agent in agent_list:
    if agent not in agent_outliers:
        agent_outliers[agent] = []
    
    fig = plt.figure()
    for histone in cate_histones['WT+mut']:
        if '/' in histone:
            continue
        
        dscore = agent_histone_dscore[agent][histone]
        pvalue = agent_histone_pvalue[agent][histone]
        if pvalue < 0.01:
            plt.plot(dscore, -np.log(pvalue), 'r.')
            plt.annotate(histone, (dscore, -np.log(pvalue)))
            agent_outliers[agent].append(histone)
        else:
            plt.plot(dscore, -np.log(pvalue), 'k.')
            
    plt.axvline(x=0, linestyle='--', color='k', alpha=0.5)
    plt.axhline(y=0, linestyle='--', color='k', alpha=0.5)
    plt.xlabel("$\Delta$ Score")
    plt.ylabel("-log(p-value)")
    plt.title(agent)
    #plt.show()
    plt.close()

#sys.exit(1)

# compare between condensing agents
histone_list = all_good_histones
pair_corr = {}
corr_matrix = [ [np.nan for i in range(len(agent_list))] for i in range(len(agent_list)) ]
for i in range(len(agent_list)-1):
    for j in range(i+1, len(agent_list)):
        agent1 = agent_list[i]
        agent2 = agent_list[j]

        X, Y = [], []
        for histone in histone_list:
            X.append(agent_histone_score[agent1][histone])
            Y.append(agent_histone_score[agent2][histone])

        corr = scipy.stats.spearmanr(X, Y)[0]
        #print ("%s VS %s: %1.2f" % (agent1, agent2, corr))
            
        fig = plt.figure()
        plt.plot(X, Y, '.')
        plt.annotate("Spearman %1.2f" % (corr), xy=(0.2, 0.75), fontsize=12, xycoords='axes fraction')
        plt.title("%s VS %s" % (agent1, agent2))
        plt.xlabel("fold change (%s)" % (agent1))
        plt.ylabel("fold change (%s)" % (agent2))
        #plt.xscale('log', base=2)
        #plt.yscale('log', base=2)
        #plt.show()
        plt.close()

        pair_corr[(agent1, agent2)] = corr
        corr_matrix[i][j] = corr
        corr_matrix[j][i] = corr


fig, axes = plt.subplots(figsize=(3,2.4), nrows=len(agent_list), ncols=len(agent_list))

for i in range(len(agent_list)):
    for j in range(len(agent_list)):
        idx = len(agent_list)*i + j
        agent1, agent2 = agent_list[i], agent_list[j]
        if i > j:
            X = [ agent_histone_score[agent1][histone] for histone in histone_list ]
            Y = [ agent_histone_score[agent2][histone] for histone in histone_list ]
            axes[i,j].plot(X, Y, 'k.', markersize=0.5)
            wspace = 0.1*(max(X) - min(X))
            hspace = 0.1*(max(Y) - min(Y))
            axes[i,j].set_xticks([min(X), max(X)])
            axes[i,j].set_xticklabels([str(round(min(X),1)), str(round(max(X),1))], rotation=45)
            axes[i,j].set_yticks([min(Y), max(Y)])
            axes[i,j].set_yticklabels([str(round(min(Y),1)), str(round(max(Y),1))])
            axes[i,j].set_xlim(min(X)-wspace, max(X)+wspace)
            axes[i,j].set_ylim(min(Y)-hspace, max(Y)+hspace)
            axes[i,j].tick_params(axis='both', which='major', labelsize=5)
            axes[i,j].tick_params(axis='both', which='minor', labelsize=5)
            #axes[i,j].set_xscale('log', base=2)
            #axes[i,j].set_yscale('log', base=2)
            if j > 0 and i < len(agent_list) -1:
                axes[i,j].tick_params(axis='both', which='both', labelbottom=False, labelleft=False)
            if j == 0 and i < len(agent_list) -1:
                axes[i,j].tick_params(axis='x', which='both', labelbottom=False)
            if j > 0 and i == len(agent_list) - 1:
                axes[i,j].tick_params(axis='y', which='both', labelleft=False)
                
        elif i == j:
            matrix = np.zeros((len(agent_list), len(agent_list)))
            matrix[:] = np.nan
            axes[i,j].imshow(matrix, origin='lower')
            s = agent1
            axes[i,j].text(len(agent_list)/2, len(agent_list)/2, s, ha="center", va="center", fontsize=7, weight='bold')
            axes[i,j].set_xlim([0, len(agent_list)-1])
            axes[i,j].set_ylim([0, len(agent_list)-1])
            axes[i,j].set_axis_off()
        else:
            assert i < j
            value = pair_corr[(agent1, agent2)]
            matrix = np.zeros((len(agent_list), len(agent_list)))
            matrix[:] = value
            img = axes[i,j].imshow(matrix, cmap="jet", vmin=0.1, vmax=0.9, origin='lower')
            if value < 0.4:
                color = "white"
            else:
                color = "black"
            axes[i,j].text(len(agent_list)/2, len(agent_list)/2, str(round(value,2)), ha="center", va="center", fontsize=6, color=color, weight='bold')
            axes[i,j].set_xlim([0, len(agent_list)-1])
            axes[i,j].set_ylim([0, len(agent_list)-1])
            axes[i,j].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labeltop=False, labelleft=False, labelright=False)

plt.subplots_adjust(wspace=0.2, hspace=0.2)
cbar=fig.colorbar(img, ax=axes, location='right', shrink=0.6, aspect=30, ticks=[0.1, 0.9])
cbar.ax.set_yticklabels([str(0.1), str(0.9)], fontsize=5)
cbar.ax.set_ylabel('Spearman correlation', rotation=-90, va="bottom", fontsize=5, labelpad=-5)
#plt.suptitle("Corrrelation betwen condensing agents")
#plt.savefig("onco_corr_btw_agent.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()

#sys.exit(1)

# check pvalue for histone body VS tail
histone_list = cate_histones['WT+mut']
subunit_foldrange = {'H2A':[17, 96], 'H2B':[38, 122], 'H3.1':[45, 131], 'H4':[31, 93]}
subunit_domain_histones = {}
for histone in histone_list:
    for subunit in histone_minfo[histone]:
        for pos in histone_minfo[histone][subunit]['mut']:
            if pos >= subunit_foldrange[subunit][0] and pos <= subunit_foldrange[subunit][1]:
                domain = 'fold'
            else:
                domain = 'tail'
            if subunit not in subunit_domain_histones:
                subunit_domain_histones[subunit] = {}
            if domain not in subunit_domain_histones[subunit]:
                subunit_domain_histones[subunit][domain] = []
            subunit_domain_histones[subunit][domain].append(histone)

agent_subunit_domain_mpvalue = {}
mpvalue_list = []
for agent in agent_list:
    for subunit in subunit_domain_histones:
        for domain in subunit_domain_histones[subunit]:
            pvalues = []
            for histone in subunit_domain_histones[subunit][domain]:
                pvalue = agent_histone_pvalue[agent][histone]
                #pvalues.append(pvalue)
                pvalues.append(1 - pvalue) # actually collect 1-p, not p
            mpvalue = np.mean(pvalues)
            if agent not in agent_subunit_domain_mpvalue:
                agent_subunit_domain_mpvalue[agent] = {}
            if subunit not in agent_subunit_domain_mpvalue[agent]:
                agent_subunit_domain_mpvalue[agent][subunit] = {}
            agent_subunit_domain_mpvalue[agent][subunit][domain] = mpvalue
            mpvalue_list.append(mpvalue)
pmin, pmax = min(mpvalue_list), max(mpvalue_list)

def uniform_circle (R, repeat):
    theta = np.random.uniform(0,2*np.pi, repeat)
    radius = np.random.uniform(0,R, repeat) ** 0.5
    X = radius * np.cos(theta)
    Y = radius * np.sin(theta)
    return X, Y

tuple_subunit = {(0,0):'H2A', (0,1):'H2B', (1,0):'H3.1', (1,1):'H4'}
tuple_sign = {(0,0):(1,-1), (0,1):(-1,-1), (1,0):(1,1), (1,1):(-1,1)}
agent_color = {'sp':'tab:blue', 'spd':'tab:orange', 'CoH':'tab:green', 'PEG':'tab:red', 'HP1a':'tab:purple'}
size = 20

fig, axes = plt.subplots(figsize=(3.8,3.4), nrows=2, ncols=2)
for i in range(2):
    for j in range(2):
        subunit = tuple_subunit[(i,j)]
        sign = tuple_sign[(i,j)]
        for domain in ['fold', 'tail']:
            if domain == 'fold':
                xc, yc = sign[0]*0.25*size, sign[1]*0.25*size
            else:
                xc, yc = -sign[0]*0.25*size, -sign[1]*0.25*size

            for agent in agent_list:
                mpvalue = agent_subunit_domain_mpvalue[agent][subunit][domain]

                R = (0.5)*size
                dX, dY = uniform_circle(R, 1)
                x = xc + dX[0]
                y = yc + dY[0]

                s = rescale(mpvalue, pmin, pmax, 0.5, 20*size)
                alpha = rescale(mpvalue, pmin, pmax, 0, 1)                
                color = agent_color[agent]
                axes[i][j].scatter([x], [y], s=[s], c=color, alpha=alpha)
                
        axes[i][j].set_xlim([-0.5*size, 0.5*size])
        axes[i][j].set_ylim([-0.5*size, 0.5*size])
        axes[i][j].spines['left'].set_visible(False)
        axes[i][j].spines['right'].set_visible(False)
        axes[i][j].spines['top'].set_visible(False)
        axes[i][j].spines['bottom'].set_visible(False)
        axes[i][j].tick_params(left='off', labelleft='off', right='off', labelright='off',
                               top='off', labeltop='off', bottom='off', labelbottom='off')

plt.subplots_adjust(wspace=0.02, hspace=0.02)
#plt.savefig("foldvstail_bubble.svg", format="svg", bbox_inches="tight")
#plt.show()
plt.close()
#sys.exit(1)


### Accumulation plot (all single WT mutations)
histone_list = cate_histones['WT+mut']
subunit_loc_minfos = {}
for histone in histone_list:
    subunit, mutation = histone.split(' ')
    if '/' in mutation:
        continue
    before, pos, after = mutation[0], mutation[1:-1], mutation[-1]
    pos = int(pos)
    minfo = (before, after)
    if subunit not in subunit_loc_minfos:
        subunit_loc_minfos[subunit] = {}
    if pos not in subunit_loc_minfos[subunit]:
        subunit_loc_minfos[subunit][pos] = []
    subunit_loc_minfos[subunit][pos].append(minfo)
    
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
agent_color = {'sp':'tab:blue', 'spd':'tab:orange', 'CoH':'tab:green', 'PEG':'tab:red', 'HP1a':'tab:purple'}

aa_type = {'C':'ap', 'D':'neg', 'S':'p', 'Q':'p', 'K':'pos',
           'I':'ap', 'P':'ap', 'T':'p', 'F':'ap', 'N':'p', 
           'G':'ap', 'H':'pos', 'L':'ap', 'R':'pos', 'W':'ap', 
           'A':'ap', 'V':'ap', 'E':'neg', 'Y':'p', 'M':'ap'}

aatype_supscript = {'ap':'', 'p':'0', 'pos':'+', 'neg':'-'}
aatype_color = {'ap':'black', 'p':'tab:green', 'pos':'tab:blue', 'neg':'tab:red'}

# group by location
for subunit in subunit_loc_minfos:
    loc_minfos = subunit_loc_minfos[subunit]
    location_list = sorted(loc_minfos.keys())
    width_ratios = [len(loc_minfos[loc]) for loc in location_list]
    width = int(round(sum(width_ratios)*0.195))
    fig, axes = plt.subplots(figsize=(width, 1.2), nrows=1, ncols=len(location_list), sharey=True,
    gridspec_kw={'width_ratios':width_ratios})
    #fig, axes = plt.subplots( nrows=1, ncols=len(location_list), sharey=True, gridspec_kw={'width_ratios':width_ratios})
    for i in range(len(location_list)):
        pos = location_list[i]
        xticks = []
        for j in range(len(loc_minfos[pos])):
            minfo = sorted(loc_minfos[pos])[j]
            before, after = minfo
            #text = '${\\mathrm{%s}}^{%s}\\rightarrow{\\mathrm{%s}}^{%s}$' % (before,
                                                                             #aatype_supscript[aa_type[before]],
                                                                             #after,
                                                                             #aatype_supscript[aa_type[after]])
            #text = '${%s}^{%s}\\rightarrow{%s}^{%s}$' % (before, aatype_supscript[aa_type[before]],
                                                         #after, aatype_supscript[aa_type[after]])

            text = '${\\mathrm{%s}}^{%s}$' % (after, aatype_supscript[aa_type[after]])
            xticks.append(text)
            mutation = ''.join([before, str(pos), after])
            histone = ' '.join([subunit, mutation])
            positive, negative = 0.0, 0.0
            for k in range(len(agent_list)):
                agent = agent_list[k]
                dscore = agent_histone_dscore[agent][histone]
                if dscore > 0:
                    axes[i].bar(j, dscore, bottom=positive, color=agent_color[agent], label=agent, width=0.5)
                    positive += dscore
                elif dscore < 0:
                    axes[i].bar(j, dscore, bottom=negative, color=agent_color[agent], label=agent, width=0.5)
                    negative += dscore

        if i > 0:
            axes[i].spines['left'].set_visible(False)
            axes[i].tick_params(left=False)
            
        axes[i].spines['top'].set_visible(False)
        axes[i].spines['right'].set_visible(False)
        axes[i].axhline(y=0, linestyle='-', color='k', alpha=1.0, lw=0.7)

        #if histone in APmutants:
        #    axes[i].set_title(str(pos), fontsize=6, weight='bold', color='red')
        #else:
        #    axes[i].set_title(str(pos), fontsize=6, color=aatype_color[aa_type[before]])

        text = str(pos) + '\n' + '${\\mathrm{%s}}^{%s}$' % (before, aatype_supscript[aa_type[before]])
        color = aatype_color[aa_type[before]]
        axes[i].set_title(text, fontsize=6, weight='bold', color=color, ha='center')

        axes[i].set_xticks(range(len(xticks)))
        #axes[i].set_xticklabels(xticks, fontsize=6, rotation=50, ha="right", rotation_mode="anchor")
        axes[i].set_xticklabels(xticks, fontsize=6)
        axes[i].set_xlim([-0.5, j+0.5])
        axes[i].set_ylim([-4.8, 3])
        #axes[i].set_ylim([-5, 3])
        if i == 0:
            axes[i].set_ylabel('Accumulated Score', fontsize=5)
        axes[i].tick_params(axis='y', which='major', labelsize=5)
        axes[i].tick_params(axis='y', which='minor', labelsize=5)

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    #ax.legend(handles[:len(agent_list)], labels[:len(agent_list)], loc='lower right', frameon=False)
    #plt.legend()
    #plt.tight_layout()
    #plt.savefig("single_WTmutant_%s.svg" % (subunit), format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()
    
#sys.exit(1)

aa_type = {'C':'ap', 'D':'neg', 'S':'p', 'Q':'p', 'K':'pos',
           'I':'ap', 'P':'ap', 'T':'p', 'F':'ap', 'N':'p', 
           'G':'ap', 'H':'pos', 'L':'ap', 'R':'pos', 'W':'ap', 
           'A':'ap', 'V':'ap', 'E':'neg', 'Y':'p', 'M':'ap'}
aatype_supscript = {'ap':'', 'p':'0', 'pos':'+', 'neg':'-'}

aatype_vec = {'ap':(0, 0), 'p':(0, 1), 'pos':(1, 1), 'neg':(-1, 1)}
aatype_color = {'ap':'black', 'p':'tab:green', 'pos':'tab:blue', 'neg':'tab:red'}


# group by transition type and plot transition matrix
histone_list = cate_histones['WT+mut']
before_after_locations = {}
for histone in histone_list:
    subunit, mutation = histone.split(' ')
    if '/' in mutation:
        continue
    before, pos, after = mutation[0], mutation[1:-1], mutation[-1]
    location = (subunit, int(pos))
    if before not in before_after_locations:
        before_after_locations[before] = {}
    if after not in before_after_locations[before]:
        before_after_locations[before][after] = []
    before_after_locations[before][after].append(location)

before_after_agent_mpvalue = {}
for before in before_after_locations:
    for after in before_after_locations[before]:
        for agent in agent_list:
            locations = before_after_locations[before][after]
            pvalues = []
            for subunit, pos in locations:
                mutation = ''.join([before, str(pos), after])
                histone = ' '.join([subunit, mutation])
                pvalue = agent_histone_pvalue[agent][histone]
                pvalues.append(pvalue)
            mpvalue = np.mean(pvalues)
            if before not in before_after_agent_mpvalue:
                before_after_agent_mpvalue[before] = {}
            if after not in before_after_agent_mpvalue[before]:
                before_after_agent_mpvalue[before][after] = {}
            before_after_agent_mpvalue[before][after][agent] = mpvalue
            

def aa_text (aa):
    aatype = aa_type[aa]
    supscript = aatype_supscript[aatype]
    #return aa
    return '${\\mathrm{%s}}^{\\mathrm{%s}}$' % (aa, supscript)

def aa_color (aa):
    return aatype_color[aa_type[aa]]

def aa_cmp (aa1, aa2):
    aatype1 = aa_type[aa1]
    aatype2 = aa_type[aa2]
    if aatype1 != 'ap' and aatype2 == 'ap':
        return -1
    elif aatype1 == 'ap' and aatype2 != 'ap':
        return 1
    elif aatype1 != 'ap' and aatype2 != 'ap':
        vec1 = aatype_vec[aatype1]
        vec2 = aatype_vec[aatype2]
        if vec1 < vec2:
            return -1
        elif vec1 > vec2:
            return 1
        return 0
    else:
        if aa1 < aa2:
            return -1
        elif aa1 > aa2:
            return 1
        return 0

"""
idx_before = sorted(before_after_locations.keys(), cmp=aa_cmp)
idx_after = set([])
for before in idx_before:
    idx_after |= set(before_after_locations[before].keys())
idx_after = sorted(list(idx_after), cmp=aa_cmp)

fig, axes = plt.subplots(figsize=(7.5, 10), nrows=len(idx_after), ncols=len(idx_before))
for i in range(len(idx_after)):
    after = idx_after[i]
    for j in range(len(idx_before)):
        before = idx_before[j]
        try:
            agent_mpvalue = before_after_agent_mpvalue[before][after]
            r = 1
            axes[i][j].plot(0, 0, '.', ms=r**2, color='white', zorder=100)
            for k in range(len(agent_list)):
                agent = agent_list[k]
                mpvalue = agent_mpvalue[agent]
                #r += 4*mpvalue
                r += (1-mpvalue)
                #r += -np.log2(mpvalue)*0.1
                axes[i][j].plot(0, 0, '.', ms=r**2, color=agent_color[agent], zorder=-k)
        except:
            pass

        axes[i][j].spines['left'].set_visible(False)
        axes[i][j].spines['right'].set_visible(False)
        axes[i][j].spines['top'].set_visible(False)
        axes[i][j].spines['bottom'].set_visible(False)
        axes[i][j].tick_params(left='off', labelleft='off', right='off', labelright='off',
                               top='off', labeltop='off', bottom='off', labelbottom='off')
        axes[i][j].set_xlim([-50, 50])
        axes[i][j].set_ylim([-50, 50])

        if i == 0:
            axes[i][j].spines['top'].set_visible(True)
            axes[i][j].tick_params(labeltop='on')
            axes[i][j].set_xticks([0])
            axes[i][j].set_xticklabels([aa_text(before)], color=aa_color(before))

            for ticklabel in axes[i][j].get_xticklabels():
                ticklabel.set_fontweight('bold')

        if j == 0:
            axes[i][j].spines['left'].set_visible(True)
            axes[i][j].tick_params(labelleft='on')
            axes[i][j].set_yticks([0])
            axes[i][j].set_yticklabels([aa_text(after)], color=aa_color(after))

            for ticklabel in axes[i][j].get_yticklabels():
                ticklabel.set_fontweight('bold')


#plt.show()
#plt.savefig("WTtransition_aa.svg", format='svg', bbox_inches='tight')
plt.close()

#sys.exit(1)
    


    
# group by mutation type
histone_list = cate_histones['WT+mut']
mtype_mut_locations = {}                             
for histone in histone_list:
    subunit, mutation = histone.split(' ')
    if '/' in mutation:
        continue
    before, pos, after = mutation[0], mutation[1:-1], mutation[-1]
    mtype = (aa_type[before], aa_type[after])
    location = (subunit, int(pos))
    mut = (before, after)
    if mtype not in mtype_mut_locations:
        mtype_mut_locations[mtype] = {}
    if mut not in mtype_mut_locations[mtype]:
        mtype_mut_locations[mtype][mut] = []
    mtype_mut_locations[mtype][mut].append(location)

for mtype in mtype_mut_locations:
    mut_locations = mtype_mut_locations[mtype]
    mutations = mut_locations.keys()

    width_ratios = [len(mut_locations[mut]) for mut in mutations]
    width = sum(width_ratios)*0.15

    fig, axes = plt.subplots(figsize=(width, 0.6), nrows=1, ncols=len(mutations), sharey=True, gridspec_kw={'width_ratios':width_ratios})

    for i in range(len(mutations)):

        if len(mutations) == 1:
            ax = axes
        else:
            ax = axes[i]
        
        mut = mutations[i]
        before, after = mut
        locations = sorted(mut_locations[mut])
        xticks = []
        for j in range(len(locations)):
            subunit, pos = locations[j]
            mutation = ''.join([before, str(pos), after])
            histone = ' '.join([subunit, mutation])
            xticks.append(' '.join([subunit, str(pos)]))
            positive, negative = 0.0, 0.0
            for k in range(len(agent_list)):
                agent = agent_list[k]
                dscore = agent_histone_dscore[agent][histone]
                if dscore > 0:
                    ax.bar(j, dscore, bottom=positive, color=agent_color[agent], label=agent, width=0.5)
                    positive += dscore
                elif dscore < 0:
                    ax.bar(j, dscore, bottom=negative, color=agent_color[agent], label=agent, width=0.5)
                    negative += dscore

        if i > 0:
            ax.spines['left'].set_visible(False)
            ax.tick_params(left=False)
            
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.axhline(y=0, linestyle='-', color='k', alpha=1.0, lw=0.7)

        before_sup = aatype_supscript[aa_type[before]]
        after_sup = aatype_supscript[aa_type[after]]
        text = '${\\mathrm{%s}}^{%s}\\rightarrow{\\mathrm{%s}}^{%s}$' % (before,before_sup,after,after_sup)
        #text = str(pos) + '\n' + '${\\mathrm{%s}}^{%s}$' % (before, aatype_supscript[aa_type[before]])
        #color = aatype_color[aa_type[before]]
        ax.set_title(text, fontsize=5, weight='bold', rotation=45, ha="left", rotation_mode="anchor")

        ax.set_xticks(range(len(xticks)))
        ax.set_xticklabels(xticks, fontsize=5, rotation=50, ha="right", rotation_mode="anchor")
        #ax.set_xticklabels(xticks, fontsize=6)
        ax.set_xlim([-0.5, j+0.5])
        ax.set_ylim([-4.8, 3])
        #ax.set_ylim([-5, 3])
        if i == 0:
            ax.set_ylabel('Accumulated Score', fontsize=5)
        ax.tick_params(axis='y', which='major', labelsize=5)
        ax.tick_params(axis='y', which='minor', labelsize=5)

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=(width-5)*0.1, hspace=None)
    #ax = plt.gca()
    #handles, labels = ax.get_legend_handles_labels()
    #ax.legend(handles[:len(agent_list)], labels[:len(agent_list)], loc='lower right', frameon=False)
    #plt.legend()
    #plt.tight_layout()
    #plt.savefig("single_WTmutant_%s%s.svg" % (mtype[0], mtype[1]), format='svg', bbox_inches='tight)'
    #plt.show()
    plt.close()


#sys.exit(1)
"""

# plot ranking bar (pretty version)
aa_type = {'C':'ap', 'D':'neg', 'S':'p', 'Q':'p', 'K':'pos',
           'I':'ap', 'P':'ap', 'T':'p', 'F':'ap', 'N':'p', 
           'G':'ap', 'H':'pos', 'L':'ap', 'R':'pos', 'W':'ap', 
           'A':'ap', 'V':'ap', 'E':'neg', 'Y':'p', 'M':'ap'}
aatype_supscript = {'ap':'', 'p':'0', 'pos':'+', 'neg':'-'}

aatype_vec = {'ap':(0, 0), 'p':(0, 1), 'pos':(1, 1), 'neg':(-1, 1)}
dcomp1_color = {-2:'purple', -1:'tab:red', 0:'black', 1:'tab:blue', 2:'blue'}
dcomp2_shapealphasize = {-1:('D',1, None), 0:('o',0, 7), 1:('o',1, 7)}
#aatype_color = {'ap':'black', 'p':'purple', 'pos':'blue', 'neg':'red'}
#aatype_color = {'ap':'black', 'p':'tab:green', 'pos':'tab:blue', 'neg':'tab:red'}

def aa_text (aa):
    aatype = aa_type[aa]
    supscript = aatype_supscript[aatype]
    return aa
    #return '${\\mathrm{%s}}^{\\mathrm{%s}}$' % (aa, supscript)

def mutation_text (mutation):
    pre_aa, post_aa = mutation
    aatype1 = aa_type[pre_aa]
    aatype2 = aa_type[post_aa]
    vec1 = aatype_vec[aatype1]
    vec2 = aatype_vec[aatype2]
    dcomp1 = vec2[0] - vec1[0]
    dcomp2 = vec2[1] - vec1[1]
    color = dcomp1_color[dcomp1]
    shape, alpha, size = dcomp2_shapealphasize[dcomp2]
    #if pre_aa in ['Y', 'R'] or post_aa in ['Y', 'R']:
    #    return aa_text(post_aa), 'red', 'D', 1, None
    #return aa_text(post_aa), 'k', 'o', 0, None
    return aa_text(post_aa), color, shape, alpha, size

histone_list = cate_histones['WT+mut']
additional_list = cate_histones['H2A.Z'] + cate_histones['H2A.V'] + cate_histones['H3.3']
additional_list += cate_histones['freeDNA']

subunits = ['H2A', 'H2B', 'H3.1', 'H4']
AP_info = {"H2A":["E56", "E61", "E64", "D90", "E91", "E92"], "H2B":["E105", "E113"]}
subunit_foldrange = {'H2A':[17, 96], 'H2B':[38, 122], 'H3.1':[45, 131], 'H4':[31, 93]}
subunit_idx_posaa = {}
AP_histones = set([])
for histone in histone_list:
    if histone.startswith('WT'):
        continue
    minfo = histone_minfo[histone]
    for subunit in minfo:
        for pos in minfo[subunit]['mut']:
            pre_aa, post_aa = minfo[subunit]['mut'][pos]

            if subunit not in subunit_idx_posaa:
                subunit_idx_posaa[subunit] = set([])
            subunit_idx_posaa[subunit].add((pos, pre_aa))

            try:
                assert pre_aa + str(pos) in AP_info[subunit]
                AP_histones.add(histone)
            except:
                pass
        
subunit_pos_idx = {}
for subunit in subunit_idx_posaa:
    idx_posaa = sorted(list(subunit_idx_posaa[subunit]))
    subunit_idx_posaa[subunit] = idx_posaa
    pos_idx = {idx_posaa[k][0]:k for k in range(len(idx_posaa))}
    subunit_pos_idx[subunit] = pos_idx

AP_histones = list(AP_histones)

for agent in agent_list:
    dscore_histone = []
    for histone in histone_list + additional_list:
        dscore = agent_histone_dscore[agent][histone]
        dscore_histone.append([dscore, histone])
    dscore_histone = sorted(dscore_histone)

    subunit_table = {}
    sorted_dscore_list = []
    sorted_histone_list = []
    for dscore, histone in dscore_histone:
        if histone not in additional_list:
            minfo = histone_minfo[histone]
            for subunit in subunits:
                pos_idx = subunit_pos_idx[subunit]
                row = [[] for k in range(len(pos_idx))]
                if subunit in minfo:
                    for pos in minfo[subunit]['mut']:
                        pre_aa, post_aa = minfo[subunit]['mut'][pos]
                        idx = pos_idx[pos]
                        row[idx].append((pre_aa, post_aa))    
                if subunit not in subunit_table:
                    subunit_table[subunit] = []
                subunit_table[subunit].append(row)
        else:
            for subunit in subunits:
                pos_idx = subunit_pos_idx[subunit]
                row = [[] for k in range(len(pos_idx))]
                if subunit not in subunit_table:
                    subunit_table[subunit] = []
                subunit_table[subunit].append(row)

        sorted_histone_list.append(histone)
        sorted_dscore_list.append(dscore)

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(7,11), sharey=True, gridspec_kw={'width_ratios': [8, 0.8]})

    # plot mutation table
    offset = 0
    space = 1
    xticks = []
    xticklabels1, xticklabels2 = [], []
    vlines = []
    for i in range(len(subunits)):
        subunit = subunits[i]
        table = subunit_table[subunit]
        
        for j in range(len(table)):
            row = table[j]
            for k in range(len(row)):
                if len(row[k]) > 0:
                    mutation = row[k][0]
                    text, color, shape, alpha, size = mutation_text(mutation)
                    axes[0].plot([offset+k], [j], shape, alpha=alpha,
                                 color='k', mfc='white', mew= 0.5, ms=size)
                    axes[0].annotate(text, xy=(offset+k, j), ha='center', va='center',
                                     fontsize=5, color=color, weight='bold')

                    pos, pre_aa = subunit_idx_posaa[subunit][k]
                    try:
                        assert pre_aa + str(pos) in AP_info[subunit]
                        vlines.append(offset+k)
                    except:
                        pass

        idx_posaa = subunit_idx_posaa[subunit]
        xticks += [offset + idx for idx in range(len(idx_posaa))]
        xticklabels1 += [str(pos) for pos, aa in idx_posaa]
        xticklabels2 += [aa_text(aa) for pos, aa in idx_posaa]

        offset += len(idx_posaa)
        if i < len(subunits)-1:
            offset += space
    
    axes[0].set_xlim([-1, offset])
    axes[0].invert_yaxis()

    ax = axes[0].twiny()
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels2, fontsize=5)
        
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    ax.tick_params(left='off', labelleft='off', right='off', labelright='off', top='off', labeltop='off', bottom='off', labelbottom='on')
    ax.set_xlim([-1, offset])
  
    axes[0].set_xticks(xticks)

    axes[0].set_xticklabels(xticklabels1, fontsize=5, rotation=90)
    axes[0].set_yticks(range(len(sorted_histone_list)))

    yticklabels = []
    for histone in sorted_histone_list:
        if agent_histone_pvalue[agent][histone] < 0.01:
            yticklabels.append('$\\ast$'+histone)
        else:
            yticklabels.append(histone)
    
    #axes[0].set_yticklabels(sorted_histone_list, fontsize=5)
    axes[0].set_yticklabels(yticklabels, fontsize=5)
    
    axes[0].spines['right'].set_visible(False)
    axes[0].tick_params(right='off', labelright='off', top='on', labeltop='on', bottom='on', labelbottom='off')

    #for ticklabel in axes[0].get_xticklabels():
    #    pos = int(ticklabel.get_text())
    #    #color = aatype_color[aa_type[aa]]
    #    #ticklabel.set_color(color)
    #    ticklabel.set_fontweight('bold')

    for ticklabel in axes[0].get_yticklabels():
        histone = ticklabel.get_text().strip('$\\ast$')
        if histone in AP_histones:
            ticklabel.set_fontweight('bold')
        if histone in cate_histones['H2A.Z']:
            ticklabel.set_color('m')
            ticklabel.set_fontweight('bold')
        if histone in cate_histones['H2A.V']:
            ticklabel.set_color('tab:green')
            ticklabel.set_fontweight('bold')
        if histone in cate_histones['H3.3']:
            ticklabel.set_color('tab:blue')
            ticklabel.set_fontweight('bold')
        if histone in cate_histones['freeDNA']:
            ticklabel.set_color('tab:red')
            ticklabel.set_fontweight('bold')
        

    for ticklabel in ax.get_xticklabels():
        aa = ticklabel.get_text()
        #color = aatype_color[aa_type[aa]]
        #ticklabel.set_color(color)
        ticklabel.set_fontweight('bold')

    for vline in vlines:
        axes[0].axvline(x=vline, color='k', linestyle='dotted', zorder=-100, alpha=0.2, linewidth=0.8)

    for u in range(len(sorted_histone_list)):
        histone = sorted_histone_list[u]
        if histone in AP_histones:
            axes[0].axhline(y=u, color='k', linestyle='dotted', zorder=-100, alpha=0.4, linewidth=1)

    # plot score bar graph
    barlist = axes[1].barh(range(len(sorted_dscore_list)), sorted_dscore_list, height=0.7)
    for bar, dscore in zip(barlist, sorted_dscore_list):
        if dscore >=0:
            color = 'tab:blue'
        else:
            color = 'tab:red'
        bar.set_color(color)
        #bar.set_edgecolor('k')
        #bar.set_linewidth(0.7)

    space = abs(max(sorted_dscore_list) - min(sorted_dscore_list))*0.1
    for u in range(len(sorted_histone_list)):
        histone = sorted_histone_list[u]
        dscore = sorted_dscore_list[u]
        if agent_histone_pvalue[agent][histone] < 0.01:
            if dscore > 0 :
                offset = space
            else:
                offset = -space
            axes[1].annotate('$\\ast$', (dscore+offset, u), ha='center', va='center', weight='bold', size=4)
        
    axes[1].tick_params(top='on', left='on', bottom='off', right='off', labeltop='on', labelleft='off', labelbottom='off')
    axes[1].tick_params(axis='both', which='major', labelsize=5)
    axes[1].tick_params(axis='both', which='minor', labelsize=5)
    axes[1].invert_yaxis()
    axes[1].axvline(x=0, color='k', linestyle='-', alpha=1, linewidth=0.8)

    axes[1].set_xticks([min(sorted_dscore_list)-2*space, 0, max(sorted_dscore_list)+2*space])
    axes[1].set_xticklabels([str(round(min(sorted_dscore_list),1)), 0, str(round(max(sorted_dscore_list),1))], fontsize=5, rotation=45, ha="left", rotation_mode="anchor")

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=None)
    plt.ylim([-1, len(sorted_dscore_list)])

    #plt.show()
    plt.savefig("%s_table_bar.svg" % (agent), format='svg', bbox_inches='tight')
    plt.close()

sys.exit(1)



# check cation-pi interactions
cat_aas = ['R', 'K']
pi_aas = ['F', 'Y']
histone_list = cate_histones['WT+mut']

histone_dnum = {}
for histone in histone_list:
    minfo = histone_minfo[histone]
    dnum = [0, 0]
    for subunit in minfo:
        for pos in minfo[subunit]['mut']:
            pre_aa, post_aa = minfo[subunit]['mut'][pos]

            # check cat aa change
            if pre_aa in cat_aas and post_aa not in cat_aas:
                dnum[0] -=1
            if pre_aa not in cat_aas and post_aa in cat_aas:
                dnum[0] +=1

            # check pi aa change
            if pre_aa in pi_aas and post_aa not in pi_aas:
                dnum[1] -=1
            if pre_aa not in pi_aas and post_aa in pi_aas:
                dnum[1] +=1

    histone_dnum[histone] = dnum

for agent in agent_list:
    dscore_histone = []
    for histone in histone_list:
        dscore = agent_histone_dscore[agent][histone]
        dscore_histone.append([dscore, histone])
    dscore_histone = sorted(dscore_histone)

    X, Y = [], []
    C = []
    for dscore, histone in dscore_histone:
        dnum = histone_dnum[histone]
        X.append(dnum[0]+random.gauss(0, 0.2))
        Y.append(dnum[1]+random.gauss(0, 0.2))
        C.append(dscore)

    fig = plt.figure()
    plt.scatter(X, Y, c=C, cmap='jet', s=5, alpha=1)
    plt.axvline(x=0, linestyle='--', color='k')
    plt.axhline(y=0, linestyle='--', color='k')
    plt.xlabel("Cation amino acid number change")
    plt.ylabel("Pi amino acid number change")
    plt.title(agent)
    #plt.show()
    plt.close()
            

#sys.exit(1)

# categorize histones by combining all condensing agents scores
cate_color = {'freeDNA':"tab:red",
              'tag':"saddlebrown",
              'WT':"black",
              'WT+mut':"tab:gray",
              'H2A.Z':"magenta",
              'H2A.Z+mut':"tab:pink",
              'H2A.V':"green",
              'H2A.V+mut':"tab:green",
              'H3.3':"blue",
              'H3.3+mut':"tab:blue"}
cate_marker = {'freeDNA':"X",
              'tag':"D",
              'WT':"o",
              'WT+mut':"o",
              'H2A.Z':"^",
              'H2A.Z+mut':"^",
              'H2A.V':"s",
              'H2A.V+mut':"s",
              'H3.3':"P",
              'H3.3+mut':"P"}

#cate_list = ['WT', 'WT+mut', 'H2A.Z', 'H2A.Z+mut', 'H2A.V', 'H2A.V+mut', 'H3.3', 'H3.3+mut', 'tag', 'freeDNA']
cate_list = ['WT', 'WT+mut', 'H2A.Z', 'H2A.Z+mut', 'H2A.V', 'H2A.V+mut', 'H3.3', 'H3.3+mut', 'freeDNA']
#histone_list = all_good_histones
histone_list = list(set(all_good_histones) - set(cate_histones['tag']))

histone_state = {}
for agent in agent_list:
    for histone in histone_list:
        score = agent_histone_score[agent][histone]
        if histone not in histone_state:
            histone_state[histone] = []
        histone_state[histone].append(score)


X = [histone_state[histone] for histone in histone_list]
        
#scaler = preprocessing.StandardScaler(with_std=False).fit(X)
scaler = preprocessing.StandardScaler().fit(X)
X_scaled = scaler.transform(X)
pca = PCA().fit(X_scaled)
Xr = pca.transform(X_scaled)

#pca = PCA().fit(X)
#Xr = pca.transform(X)

variance_ratio_list = 100* pca.explained_variance_ratio_

fig = plt.figure()
plt.plot(range(1, len(variance_ratio_list)+1), variance_ratio_list, '.-')
plt.xlabel("PCA component")
plt.ylabel("Variance (%)")
plt.xticks(range(1, len(variance_ratio_list)+1))
#plt.show()
plt.close()

clf = LocalOutlierFactor()
#outcheck = clf.fit_predict(X)
outcheck = clf.fit_predict( [row[:2] for row in Xr] )

outliers = []
histone_PCA = {}
x_list, y_list = [], []
for i in range(len(Xr)):
    histone = histone_list[i]
    x, y = Xr[i][0], Xr[i][1]
    histone_PCA[histone] = (x, y)
    x_list.append(x)
    y_list.append(y)
    if outcheck[i] < 0:
        outliers.append(histone)

fig = plt.figure(figsize=(8,7))
#plt.scatter(x_list, y_list, s=5, c=C, cmap='jet')
for i in range(len(Xr)):
    plt.plot(Xr[i][0], Xr[i][1], 'k.')
    if outcheck[i] < 0:
        plt.plot(Xr[i][0], Xr[i][1], 'r.')
        plt.annotate(histone_list[i], (Xr[i][0], Xr[i][1]))
plt.title("PCA plot")
#cbar = plt.colorbar()
#cbar.set_label("AT content (%)")
plt.xlabel('PC1')
plt.ylabel('PC2')
#plt.tight_layout()
#plt.show()
plt.close()

#fig = plt.figure(figsize=(8,7))
fig = plt.figure(figsize=(6,5))
for cate in cate_list:
    histones = cate_histones[cate]
    color = cate_color[cate]
    marker = cate_marker[cate]
    if cate == 'WT+mut':
        alpha=0.5
    else:
        alpha=1.0
    for i in range(len(histones)):
        histone = histones[i]
        x, y = histone_PCA[histone]
        if i == 0:
            if cate == 'WT':
                plt.plot(x, y, '.', color=color, marker=marker, alpha=alpha, label=cate, zorder=100)
            else:
                plt.plot(x, y, '.', color=color, marker=marker, alpha=alpha, label=cate)
        else:
            if cate == 'WT':
                plt.plot(x, y, '.', color=color, marker=marker, alpha=alpha, zorder=100)
            else:
                plt.plot(x, y, '.', color=color, marker=marker, alpha=alpha)
        if histone in outliers:
            plt.annotate(histone, (x,y), fontsize=6)

plt.title("PCA plot", fontsize=10)
plt.xlabel('PC1', fontsize=8)
plt.ylabel('PC2', fontsize=8)
plt.gca().tick_params(axis='both', which='major', labelsize=5)
plt.gca().tick_params(axis='both', which='minor', labelsize=5)
plt.legend(fontsize=6)
#plt.savefig("PCA_allagent.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()

#sys.exit(1)

"""
# check the WT outliers on the structure
subunit_chains = {'H2A':['C'], 'H2B':['D'], 'H3.1':['A'], 'H4':['B']}
subunit_foldrange = {'H2A':[17, 96], 'H2B':[38, 122], 'H3.1':[45, 131], 'H4':[31, 93]}
ltype_locs = {'D-T':[('H2B', 68), ('H2B', 70), ('H2B', 71), ('H2B', 76), ('H4', 91), ('H4', 92)], 
              'H3-H4':[('H3.1', 50), ('H3.1', 97), ('H4', 42)],
              'H-DNA':[('H2A', 29), ('H2A', 74), ('H2A', 77), ('H3.1', 83), ('H3.1', 118)]}
loc_ltype = {}
for ltype, locs in ltype_locs.items():
    for loc in locs:
        loc_ltype[loc] = ltype
        
AP_info = {"H2A":["E56", "E61", "E64", "D90", "E91", "E92"], "H2B":["E105", "E113"]}

for subunit in AP_info:
    for residue in AP_info[subunit]:
        loc = (subunit, int(residue[1:]))
        assert loc not in loc_ltype
        loc_ltype[loc] = 'AP'

ltype_color = {'tail':'green', 'AP':'red', 'D-T':'magenta', 'H3-H4':'pink', 'H-DNA':'blue', 'fold':'brown'}

# color by location typex
value_type = 'up'
note='_side'
for agent in agent_list:
    NCP = Molecules("1kx5")
    ltype_chain_resi_values = {}
    #chain_resi_values_up = {}
    #chain_resi_values_down = {}
    for histone in agent_outliers[agent]:
        value = agent_histone_dscore[agent][histone]
        if value_type == 'down' and value >0:
            continue
        if value_type == 'up' and value <0:
            continue
        for subunit in histone_minfo[histone]:
            for pos in histone_minfo[histone][subunit]['mut']:
                pre_aa, post_aa = histone_minfo[histone][subunit]['mut'][pos]

                ltype = 'fold'
                st, ed = subunit_foldrange[subunit]
                try:
                    ltype = loc_ltype[(subunit, pos)]
                    assert pos >= st or pos <= ed
                except:
                    if pos < st or pos > ed:
                        ltype = 'tail'

                #print (ltype)
                chains = subunit_chains[subunit]
                for chain in chains:

                    if ltype not in ltype_chain_resi_values:
                        ltype_chain_resi_values[ltype] = {}
                    if chain not in ltype_chain_resi_values[ltype]:
                        ltype_chain_resi_values[ltype][chain] = {}
                    if pos not in ltype_chain_resi_values[ltype][chain]:
                        ltype_chain_resi_values[ltype][chain][pos] = []
                    ltype_chain_resi_values[ltype][chain][pos].append(value)


    NCP.remove_ions()
    NCP.stylish()
    NCP.coloring(NCP.chain_resi_resn, 'white')
    for ltype in ltype_chain_resi_values:
        chain_resi_values = ltype_chain_resi_values[ltype]
        color = ltype_color[ltype]
        NCP.make_sphere(chain_resi_values)
        NCP.coloring(chain_resi_values, color)
    NCP.save_session(agent + '_' + value_type + note)
    NCP.clear_up()

sys.exit(1)


# color by up/down
for agent in agent_list:
    NCP = Molecules("1kx5")
    chain_resi_values_up = {}
    chain_resi_values_down = {}
    for histone in agent_outliers[agent]:
        value = agent_histone_dscore[agent][histone]
        for subunit in histone_minfo[histone]:
            for pos in histone_minfo[histone][subunit]['mut']:
                pre_aa, post_aa = histone_minfo[histone][subunit]['mut'][pos]
                
                chains = subunit_chains[subunit]
                for chain in chains:

                    if value > 0:        
                        if chain not in chain_resi_values_up:
                            chain_resi_values_up[chain] = {}
                        if pos not in chain_resi_values_up[chain]:
                            chain_resi_values_up[chain][pos] = []
                        chain_resi_values_up[chain][pos].append(value)

                    else:
                        if chain not in chain_resi_values_down:
                            chain_resi_values_down[chain] = {}
                        if pos not in chain_resi_values_down[chain]:
                            chain_resi_values_down[chain][pos] = []
                        assert value < 0
                        chain_resi_values_down[chain][pos].append(value)


    NCP.remove_ions()
    NCP.stylish()
    NCP.coloring(NCP.chain_resi_resn, 'white')
    NCP.make_sphere(chain_resi_values_up)
    NCP.coloring(chain_resi_values_up, 'blue')
    NCP.make_sphere(chain_resi_values_down)
    NCP.coloring(chain_resi_values_down, 'red')
    NCP.save_session(agent)
    NCP.clear_up()

#sys.exit(1)
"""

    



# plot ranking bar graph for each metrics
histone_list = all_good_histones_exceptWT
agent_outliers = all_agent_outliers
#histone_list = cate_histones['WT+mut']
#agent_outliers = WTmut_agent_outliers

cate_color = {'freeDNA':"tab:red",
              'tag':"saddlebrown",
              'WT':"black",
              'WT+mut':"black",
              'H2A.Z':"magenta",
              'H2A.Z+mut':"magenta",
              'H2A.V':"green",
              'H2A.V+mut':"green",
              'H3.3':"blue",
              'H3.3+mut':"blue"}


legend_elements = [Line2D([0], [0], color='b', lw=4, label='Line'),
                   Line2D([0], [0], marker='o', color='w', label='Scatter',
                          markerfacecolor='g', markersize=15),
                   Patch(facecolor='orange', edgecolor='r',
                         label='Color Patch')]

legend_elements = [Line2D([0], [0], marker='o', color='w', label='freeDNA', markerfacecolor='w'),
                   Line2D([0], [0], marker='o', color='w', label='tag', markerfacecolor='w'),
                   Line2D([0], [0], marker='o', color='w', label='WT+mut', markerfacecolor='w'),
                   Line2D([0], [0], marker='o', color='w', label='AP mut', markerfacecolor='w'),
                   Line2D([0], [0], marker='o', color='w', label='H2A.Z+mut', markerfacecolor='w'),
                   Line2D([0], [0], marker='o', color='w', label='H2A.V+mut', markerfacecolor='w'),
                   Line2D([0], [0], marker='o', color='w', label='H3.3+mut', markerfacecolor='w'),
                   Line2D([0], [0], marker='o', color='w', label='outliers', markerfacecolor='w')]


for agent in agent_list:
    dChalf_histone, dscore_histone = [], []
    for histone in histone_list:
        dChalf = agent_histone_dChalf[agent][histone]
        dscore = agent_histone_dscore[agent][histone]
        dChalf_histone.append([dChalf, histone])
        dscore_histone.append([dscore, histone])
    dChalf_histone = sorted(dChalf_histone)
    dscore_histone = sorted(dscore_histone)

    value_histone_list = [dChalf_histone, dscore_histone]
    ylabel_list = ["$\Delta$ C-half", "$\Delta$ Score"]
    #value_histone_list = [dscore_histone]
    #ylabel_list = ["$\Delta$ Score"]
    
    for value_histone, ylabel in list(zip(value_histone_list, ylabel_list)):        
        label, X, Y = [], [], []
        outlabel, outX, outY = [], [], []
        for i in range(len(value_histone)):
            value, histone = value_histone[i]
            if histone not in agent_outliers[agent]:
                label.append(histone)
                X.append(i)
                Y.append(value)
            else:
                outlabel.append(histone)
                outX.append(i)
                outY.append(value)

        fig = plt.figure(figsize=(18,5))
        #fig = plt.figure(figsize=(5,18))
        #plt.plot(range(len(Y)), Y, 'o-')
        plt.bar(X, Y)
        #plt.barh(X, Y)
        #plt.xticks(X, label, fontsize=6, rotation=80)
        plt.bar(outX, outY, color='r')
        #plt.barh(outX, outY, color='r')
        #plt.yticks(X+outX, label+outlabel, fontsize=6)
        plt.xticks(X+outX, label+outlabel, fontsize=6, rotation=80, ha="right", rotation_mode="anchor")

        for ticklabel in plt.gca().get_xticklabels():
            histone = ticklabel.get_text()
            cate = histone_cate[histone]
            color = cate_color[cate]
            ticklabel.set_color(color)
            if histone in APmutants:
                ticklabel.set_weight(1000)
            #if histone in APmtype_histones["Positive"]:
            #    ticklabel.set_color('b')

        #for ticklabel in plt.gca().get_xticklabels():
        #    if ticklabel.get_text() in APmutants:
        #        ticklabel.set_weight(1000)
        #    if ticklabel.get_text() in agent_outliers[agent]:
        #        ticklabel.set_color('r')
        #    if ticklabel.get_text() in APmtype_histones["Positive"]:
        #        ticklabel.set_color('b')

        ax = plt.gca()

        if ylabel == "$\Delta$ C-half":
            loc = 'upper left'
        elif ylabel == "$\Delta$ Score": 
            loc = 'lower right'
        
        leg = ax.legend(handles=legend_elements, loc=loc, frameon=False)
        for text in leg.get_texts():
            cate = text.get_text()
            if cate in cate_color:
                text.set_color(cate_color[cate])
            if cate == 'AP mut':
                text.set_weight(1000)
            if cate == 'outliers':
                text.set_bbox(dict(facecolor='red', edgecolor='w', alpha=0.7))
        
        plt.title(agent_fullname[agent])
        #plt.yscale('log', base=2)
        plt.ylabel(ylabel)
        plt.ylabel(ylabel)
        #plt.ylabel('C$_{1/2}$')
        plt.tight_layout()
        #plt.savefig("%s_%s_bar.png" % (agent, ylabel))
        #plt.show()
        plt.close()

#sys.exit(1)


# check the correlation with freeDNA contamination
histone_list = all_good_histones_exceptWT
agent_outliers = all_agent_outliers
#histone_list = cate_histones['WT+mut']
#agent_outliers = WTmut_agent_outliers
agent_histone_value_list = [agent_histone_Chalf, agent_histone_score]
ylabel_list = ['C-half', 'Score']
for agent in agent_list:
    for agent_histone_value, ylabel in list(zip(agent_histone_value_list, ylabel_list)):

        X, Y = [], []
        for histone in histone_list:
            freeDNA = histone_freeDNA[histone]
            X.append(freeDNA*100)
            value = agent_histone_value[agent][histone]
            Y.append(value)

        corr = scipy.stats.pearsonr(X, Y)[0]
        #print ("%s VS %s: %1.2f" % ('FreeDNA', agent, corr))

        feature_list = [[x] for x in X]
        test_list = [[y] for y in Y]

        #reg = linear_model.Ridge(alpha=0.5)
        #reg.fit (feature_list, test_list)
        #Ypred = reg.predict(feature_list)
        #X_Ypred = sorted([[X[i], Ypred[i][0]] for i in range(len(Ypred))])

        reg = linear_model.HuberRegressor()
        reg.fit (feature_list, Y)
        Ypred = reg.predict(feature_list)
        X_Ypred = sorted([[X[i], Ypred[i]] for i in range(len(Ypred))])        

        fig = plt.figure()
        #plt.plot(X, Y, '.')
        for histone in histone_list:
        #for histone in cate_histones['WT'] + cate_histones['WT+mut']:
            x = histone_freeDNA[histone]*100
            y = agent_histone_value[agent][histone]
            if histone not in agent_outliers[agent]:
                plt.plot(x, y, 'k.')
            else:
                plt.plot(x, y, 'r.')
                plt.annotate(histone, (x, y))
        plt.plot([X_Ypred[0][0], X_Ypred[-1][0]], [X_Ypred[0][1], X_Ypred[-1][1]], 'b--')
        #xloc, yloc = np.mean([np.median(X), max(X)]), np.mean([np.median(Y), max(Y)])
        #plt.text(xloc, yloc, str(round(corr,3)), fontsize=20, va='center', ha='center')
        plt.title(agent_fullname[agent])
        plt.xlabel("Free DNA (%)")
        #plt.ylabel("C-half")
        plt.ylabel(ylabel)
        #plt.xlim([-10, 30])
        #plt.ylim([0, 2.5])
        if ylabel == 'C-half':
            plt.yscale('log', basey=2)
        #plt.xscale('log', base=2)
        #plt.savefig("%s_%s_freeDNA.png" % (agent, ylabel))
        #plt.show()
        plt.close()

# check the correlation with input count
histone_list = all_good_histones_exceptWT
agent_outliers = all_agent_outliers
#histone_list = cate_histones['WT+mut']
#agent_outliers = WTmut_agent_outliers
for agent in agent_list:
    histone_count = agent_num_histone_count[agent][0]
    X, Y = [], []
    for histone in histone_list:
        count = histone_count[histone]
        score = agent_histone_score[agent][histone]
        X.append(count)
        Y.append(score)


    corr = scipy.stats.pearsonr(X, Y)[0]
    #print ("%s VS %s: %1.2f" % ('Input', agent, corr))

    feature_list = [[x] for x in X]
    test_list = [[y] for y in Y]

    reg = linear_model.Ridge(alpha=0.5)
    reg.fit (feature_list, test_list)
    Ypred = reg.predict(feature_list)
    X_Ypred = sorted([[X[i], Ypred[i][0]] for i in range(len(Ypred))])
    

    fig = plt.figure()
    #plt.plot(X, Y, '.')
    for histone in histone_list:
        x = histone_count[histone]
        y = agent_histone_score[agent][histone]
        if histone not in agent_outliers[agent]:
            plt.plot(x, y, 'k.')
        else:
            plt.plot(x, y, 'r.')
            plt.annotate(histone, (x, y))
    plt.plot([X_Ypred[0][0], X_Ypred[-1][0]], [X_Ypred[0][1], X_Ypred[-1][1]], 'b--')
    plt.title(agent_fullname[agent])
    plt.xlabel("Input count")
    plt.ylabel("Score")
    #plt.savefig("Input_%s.png" % (agent))
    #plt.show()
    plt.close()


# check the correlation with AT content
histone_list = all_good_histones_exceptWT
agent_outliers = all_agent_outliers
#histone_list = cate_histones['WT+mut']
#agent_outliers = WTmut_agent_outliers
for agent in agent_list:
    X, Y = [], []
    for histone in histone_list:
        AT = 100.0 - GC_content(histone_BC[histone])
        score = agent_histone_score[agent][histone]
        X.append(AT)
        Y.append(score)


    corr = scipy.stats.pearsonr(X, Y)[0]
    #print ("%s VS %s: %1.2f" % ('Input', agent, corr))

    feature_list = [[x] for x in X]
    test_list = [[y] for y in Y]

    reg = linear_model.Ridge(alpha=0.5)
    reg.fit (feature_list, test_list)
    Ypred = reg.predict(feature_list)
    X_Ypred = sorted([[X[i], Ypred[i][0]] for i in range(len(Ypred))])
    

    fig = plt.figure()
    #plt.plot(X, Y, '.')
    for histone in histone_list:
        x = 100.0 - GC_content(histone_BC[histone])
        y = agent_histone_score[agent][histone]
        if histone not in agent_outliers[agent]:
            plt.plot(x, y, 'k.')
        else:
            plt.plot(x, y, 'r.')
            plt.annotate(histone, (x, y))
    plt.plot([X_Ypred[0][0], X_Ypred[-1][0]], [X_Ypred[0][1], X_Ypred[-1][1]], 'b--')
    plt.title(agent_fullname[agent])
    plt.xlabel("AT content(%)")
    plt.ylabel("Score")
    #plt.show()
    plt.close()


#sys.exit(1)


# compare with chromatin remodeling data
histone_list = all_good_histones
fig1, axes1 = plt.subplots(nrows=len(agent_list), ncols=2)
fig2, axes2 = plt.subplots(nrows=len(agent_list), ncols=2)
for i in range(len(agent_list)):
    agent = agent_list[i]
    X1 = [histone_ACF[histone] for histone in histone_list]
    X2 = [histone_Nap1[histone] for histone in histone_list]
    Y = [agent_histone_score[agent][histone] for histone in histone_list]

    corr1 = scipy.stats.spearmanr(X1, Y)[0]
    #print ("%s VS %s: %1.2f" % ('ACF', agent, corr1))

    corr2 = scipy.stats.spearmanr(X2, Y)[0]
    #print ("%s VS %s: %1.2f" % ('Nap1', agent, corr2))

    axes1[i, 0].plot(X1, Y, '.')
    #axes1[i, 0].set_yscale('log', base=2)

    axes1[i, 1].plot(X2, Y, '.')
    #axes1[i, 1].set_yscale('log', base=2)

    axes1[i, 0].set_ylabel(agent, rotation=0, labelpad=20)
    axes1[i, 1].tick_params(axis='y', which='both', labelleft=False)

    if i >= len(agent_list) - 1:
        axes1[i, 0].set_xlabel('ACF data')
        axes1[i, 1].set_xlabel('Nap1 data')
    else:
        axes1[i,0].tick_params(axis='x', which='both', labelbottom=False)
        axes1[i,1].tick_params(axis='x', which='both', labelbottom=False)

    corrs = [corr1, corr2]
    for j in range(2):
        corr = corrs[j]
        matrix = np.zeros((len(agent_list), len(agent_list)))
        matrix[:] = corr
        img = axes2[i, j].imshow(matrix, cmap="bwr", vmin=-0.5, vmax=0.5, origin='lower')

        if corr < -0.35 or corr > 0.35:
            color = 'white'
        else:
            color = 'black'

        axes2[i,j].text(len(agent_list)/2, len(agent_list)/2, str(round(corr,2)), ha="center", va="center", fontsize=10, color=color, weight='bold')
        axes2[i,j].set_xlim([0, len(agent_list)-1])
        axes2[i,j].set_ylim([0, len(agent_list)-1])
        axes2[i,j].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    axes2[i, 0].set_ylabel(agent, rotation=0, labelpad=20)
    axes2[i, 1].tick_params(axis='y', which='both', labelleft=False)

    if i >= len(agent_list) - 1:
        axes2[i, 0].set_xlabel('ACF data')
        axes2[i, 1].set_xlabel('Nap1 data')
    else:
        axes2[i, 0].tick_params(axis='x', which='both', labelbottom=False)
        axes2[i, 1].tick_params(axis='x', which='both', labelbottom=False)

fig1.suptitle("Corrrelation betwen ACF/Nap1 data")
#plt.show()
#fig1.show(fig1)
#plt.close(fig1)

fig2.suptitle("Corrrelation betwen ACF/Nap1 data")
cbar=fig2.colorbar(img, ax=axes2, location='right', shrink=0.8)
cbar.ax.set_ylabel('Spearman correlation', rotation=-90, va="bottom")
#fig2.show()
#plt.show()
#plt.close()
#plt.close(fig2)
plt.close('all')


### Accumulation plot (AP only)
location_minfos = {}
for histone in APmutants:
    subunit, mutation = histone.split(' ')
    if '/' in mutation:
        continue
    before, pos, after = mutation[0], mutation[1:-1], mutation[-1]
    pos = int(pos)
    location = (subunit, pos)
    if location not in location_minfos:
        location_minfos[location] = []
    minfo = (before, after)
    location_minfos[location].append(minfo)

agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
location_list = sorted(location_minfos.keys())
after_color = {'K':'blue', 'A':'gray', 'D':'red', 'N':'yellow', 'Q':'yellow'}

fig, axes = plt.subplots(nrows=len(agent_list), ncols=len(location_list), sharex=True, sharey=True)

for i in range(len(agent_list)):
    agent = agent_list[i]
    for j in range(len(location_list)):
        location = location_list[j]
        subunit, pos = location
        #data = []
        #label = []
        for k in range(len(location_minfos[location])):
            minfo = sorted(location_minfos[location])[k]
            before, after = minfo
            mutation = ''.join([before, str(pos), after])
            histone = ' '.join([subunit, mutation])
            dscore = agent_histone_dscore[agent][histone]
            axes[i,j].bar(k, dscore, label=before + '->' + after, color=after_color[after])
            axes[i,j].axhline(y=0, linestyle='--', color='k', alpha=0.5)
            #data.append(dscore)
            #label.append(before + '->' + after)
        #axes[i,j].barplot(range(len(data)), data)
        #axes[i,j].legend()
        axes[i,j].tick_params(bottom='off', labelbottom='off')
        if i == len(agent_list) - 1:
            axes[i,j].set_xlabel(' '.join([subunit, str(pos)]))
        if j == 0:
            axes[i,j].set_ylabel(agent, rotation=90)

#plt.show()
plt.close()

# all single AP mutants (group by location)
agent_color = {'sp':'tab:blue', 'spd':'tab:orange', 'CoH':'tab:green', 'PEG':'tab:red', 'HP1a':'tab:purple'}
after_sign = {'K':'+', 'A':'', 'D':'-', 'N':'0', 'Q':'0'}
fig, axes = plt.subplots(figsize=(3.2, 2), nrows=1, ncols=len(location_list), sharey=True)
for i in range(len(location_list)):
    location = location_list[i]
    subunit, pos = location
    xticks = []
    for j in range(len(location_minfos[location])):
        minfo = sorted(location_minfos[location])[j]
        before, after = minfo
        xticks.append('$%s\\rightarrow{%s}^{%s}$' % (before, after, after_sign[after]))
        mutation = ''.join([before, str(pos), after])
        histone = ' '.join([subunit, mutation])
        positive, negative = 0.0, 0.0
        for k in range(len(agent_list)):
            agent = agent_list[k]
            dscore = agent_histone_dscore[agent][histone]
            if dscore > 0:
                axes[i].bar(j, dscore, bottom=positive, color=agent_color[agent], label=agent)
                positive += dscore
            elif dscore < 0:
                axes[i].bar(j, dscore, bottom=negative, color=agent_color[agent], label=agent)
                negative += dscore

    if i > 0:
        axes[i].spines['left'].set_visible(False)
        axes[i].tick_params(left=False)
    axes[i].spines['top'].set_visible(False)
    axes[i].spines['right'].set_visible(False)
    axes[i].axhline(y=0, linestyle='-', color='k', alpha=0.6)
    axes[i].set_title(' '.join([subunit, str(pos)]), fontsize=6)
    axes[i].set_xticks(range(len(xticks)))
    axes[i].set_xticklabels(xticks, fontsize=5, rotation=50, ha="right", rotation_mode="anchor")
    if i == 0:
        axes[i].set_ylabel('Accumulated Scores', fontsize=6)
        axes[i].tick_params(axis='both', which='major', labelsize=5)
        axes[i].tick_params(axis='both', which='minor', labelsize=5)

        
ax = plt.gca()
#handles, labels = ax.get_legend_handles_labels()
#ax.legend(handles[:len(agent_list)], labels[:len(agent_list)], loc='lower right', frameon=False)
#plt.legend()
#plt.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
plt.savefig("single_APmutant.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()

# all single AP mutants (group by mutation type)
after_histones = {}
for histone in APmutants:
    if '/' in histone:
        continue
    after = histone[-1]
    if after not in after_histones:
        after_histones[after] = []
    after_histones[after].append(histone)
    
after_index = {'A':0, 'K':1, 'D':2, 'Q':3, 'N':4}
fig, axes = plt.subplots(figsize=(3.2, 2), nrows=1, ncols=len(after_histones), sharey=True)
for after in after_histones:
    index = after_index[after]
    histones = sorted(after_histones[after])
    for i in range(len(histones)):
        histone = histones[i]
        positive, negative = 0.0, 0.0
        for k in range(len(agent_list)):
            agent = agent_list[k]
            dscore = agent_histone_dscore[agent][histone]
            if dscore > 0:
                axes[index].bar(i, dscore, bottom=positive, color=agent_color[agent], label=agent)
                positive += dscore
            elif dscore < 0:
                axes[index].bar(i, dscore, bottom=negative, color=agent_color[agent], label=agent)
                negative += dscore

    if index > 0:
        axes[index].spines['left'].set_visible(False)
        axes[index].tick_params(left=False)
    axes[index].spines['top'].set_visible(False)
    axes[index].spines['right'].set_visible(False)
    axes[index].axhline(y=0, linestyle='-', color='k', alpha=0.6)
    axes[index].set_title('${%s}^{%s}$ mut' % (after, after_sign[after]), fontsize=6)
    axes[index].set_xticks(range(len(histones)))
    axes[index].set_xticklabels(histones, fontsize=5, rotation=50, ha="right", rotation_mode="anchor")
    if index == 0:
        axes[index].set_ylabel('Accumulated Scores', fontsize=6)
        axes[index].tick_params(axis='both', which='major', labelsize=5)
        axes[index].tick_params(axis='both', which='minor', labelsize=5)

ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
#ax.legend(handles[:len(agent_list)], labels[:len(agent_list)], loc='lower right', frameon=False)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
#plt.tight_layout()
plt.savefig("single_APmutant2.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()

                

# check single vs mutiple AP mutants
fig = plt.figure()
histone_list = ['H2A E61A', 'H2A D90A', 'H2A E92A', 'H2A E61/D90/E92A']
for i in range(len(histone_list)):
    histone = histone_list[i]
    positive, negative = 0.0, 0.0
    for j in range(len(agent_list)):
        agent = agent_list[j]
        dscore = agent_histone_dscore[agent][histone]
        if dscore > 0:
            plt.bar(i, dscore, width=0.3, bottom=positive, color=agent_color[agent], label=agent)
            positive += dscore
        elif dscore < 0:
            plt.bar(i, dscore, width=0.3, bottom=negative, color=agent_color[agent], label=agent)
            negative += dscore
plt.axhline(y=0, linestyle='--', color='k', alpha=1.0)            
plt.ylabel('Accumulated Scores')
plt.xticks(range(len(histone_list)), histone_list, fontsize=8, rotation=40, ha="right", rotation_mode="anchor")
plt.xlim([-0.5, len(histone_list)+0.25])
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[:len(agent_list)], labels[:len(agent_list)])
plt.tight_layout()
#plt.show()
plt.close()


# Oncohistone vs PTM library AP mutant
fig, axes = plt.subplots(figsize=(3.2, 2), nrows=1, ncols=2)
histone_list = ['H2A E56A', 'H2A E61A', 'H2A D90A', 'H2A E91A', 'H2A E92A', 'H2B E113A', 'H2A E61/D90/E92A']
for i in range(len(histone_list)):
    histone = histone_list[i]
    positive, negative = 0.0, 0.0
    for j in range(len(agent_list)):
        agent = agent_list[j]
        dscore = agent_histone_dscore[agent][histone]
        if dscore > 0:
            axes[0].bar(i, dscore, width=0.4, bottom=positive, color=agent_color[agent], label=agent)
            positive += dscore
        elif dscore < 0:
            axes[0].bar(i, dscore, width=0.4, bottom=negative, color=agent_color[agent], label=agent)
            negative += dscore
axes[0].axhline(y=0, linestyle='-', color='k', alpha=0.6)            
axes[0].set_ylabel('Accumulated Scores', fontsize=6)
axes[0].set_xticks(range(len(histone_list)))
axes[0].set_xticklabels(histone_list, fontsize=5, rotation=40, ha="right", rotation_mode="anchor")
axes[0].set_title("Oncohistone library", fontsize=6)
axes[0].tick_params(axis='both', which='major', labelsize=5)
axes[0].tick_params(axis='both', which='minor', labelsize=5)


with open("PTM_dscore.pickle", "rb") as f:
    agent_ID_dscore = pickle.load(f)
    
positive, negative = 0.0, 0.0
for j in range(len(agent_list)):
    agent = agent_list[j]
    dscore = agent_ID_dscore[agent][85]
    if dscore > 0:
        axes[1].bar(0, dscore, width=0.3, bottom=positive, color=agent_color[agent], label=agent)
        positive += dscore
    elif dscore < 0:
        axes[1].bar(0, dscore, width=0.3, bottom=negative, color=agent_color[agent], label=agent)
        negative += dscore
axes[1].axhline(y=0, linestyle='-', color='k', alpha=0.6)            
axes[1].set_xticks([0])
axes[1].set_xticklabels(['AP mutant ($\\ast \\rightarrow $ A)\nH2A E56/E61/E64/D90/E91/E92A\nH2B E105/E113A'], fontsize=5, ma='left')
axes[1].set_title("PTM library", fontsize=6)
axes[1].set_xlim([-0.5, 1])
axes[1].set_ylim([-1, 8])
axes[1].tick_params(axis='both', which='major', labelsize=5)
axes[1].tick_params(axis='both', which='minor', labelsize=5)
#handles, labels = plt.gca().get_legend_handles_labels()
#plt.legend(handles[:len(agent_list)], labels[:len(agent_list)])
#plt.tight_layout()
plt.savefig("OncoVSPTMlib_APmutant.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()

# check H4 tail mutant
fig, axes = plt.subplots(figsize=(3.2, 2), nrows=1, ncols=2)
histone_list = ['H4 R17/19A']
for i in range(len(histone_list)):
    histone = histone_list[i]
    positive, negative = 0.0, 0.0
    for j in range(len(agent_list)):
        agent = agent_list[j]
        dscore = agent_histone_dscore[agent][histone]
        if dscore > 0:
            axes[0].bar(i, dscore, width=0.3, bottom=positive, color=agent_color[agent], label=agent)
            positive += dscore
        elif dscore < 0:
            axes[0].bar(i, dscore, width=0.3, bottom=negative, color=agent_color[agent], label=agent)
            negative += dscore
axes[0].axhline(y=0, linestyle='-', color='k', alpha=0.6)
axes[0].set_xlim([-0.5, 0.5])
axes[0].set_ylim([-4, 0.5])
axes[0].set_xticks([0])
axes[0].set_xticklabels(['H4 R17/19A'], fontsize=5)
axes[0].set_ylabel('Accumulated Scores', fontsize=6)
axes[0].set_title("Oncohistone library", fontsize=6)
axes[0].tick_params(axis='both', which='major', labelsize=5)
axes[0].tick_params(axis='both', which='minor', labelsize=5)


positive, negative = 0.0, 0.0
for j in range(len(agent_list)):
    agent = agent_list[j]
    dscore = agent_ID_dscore[agent][43]
    if dscore > 0:
        axes[1].bar(0, dscore, width=0.3, bottom=positive, color=agent_color[agent], label=agent)
        positive += dscore
    elif dscore < 0:
        axes[1].bar(0, dscore, width=0.3, bottom=negative, color=agent_color[agent], label=agent)
        negative += dscore
axes[1].axhline(y=0, linestyle='-', color='k', alpha=0.6)            
axes[1].set_xlim([-0.5, 0.5])
axes[1].set_ylim([-6, 0.5])
axes[1].set_xticks([0])
axes[1].set_xticklabels(['H4 R17/19A'], fontsize=5)
axes[1].set_title("PTM library", fontsize=6)
axes[1].tick_params(axis='both', which='major', labelsize=5)
axes[1].tick_params(axis='both', which='minor', labelsize=5)

ax = plt.gca()
#handles, labels = ax.get_legend_handles_labels()
#ax.legend(handles[:len(agent_list)], labels[:len(agent_list)], loc='lower right', frameon=False)
#plt.tight_layout()
plt.savefig("OncoVSPTMlib_H4tailmutant.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()
            

sys.exit(1)



### structure analysis
histone_list = cate_histones['WT+mut']
agent_outliers = WTmut_agent_outliers

subunit_chains = {'H2A':['C', 'G'], 'H2B':['D', 'H'], 'H3.1':['A', 'E'], 'H4':['B','F']}
#subunit_chains = {'H2A':['C'], 'H2B':['D'], 'H3.1':['E'], 'H4':['F']}
subunit_foldrange = {'H2A':[17, 96], 'H2B':[38, 122], 'H3.1':[45, 131], 'H4':[31, 93]}

# check tail vs fold WT mutations
tail_mutants, fold_mutants = [], []
for histone in histone_list:
    for subunit in histone_minfo[histone]:
        for pos in histone_minfo[histone][subunit]['mut']:
            if pos >= subunit_foldrange[subunit][0] and pos <= subunit_foldrange[subunit][1]:
                fold_mutants.append(histone)
            else:
                tail_mutants.append(histone)

for agent in agent_list:
    #X = [agent_histone_dscore[agent][histone] for histone in tail_mutants]
    #Y = [agent_histone_dscore[agent][histone] for histone in fold_mutants]
    X = [abs(agent_histone_dscore[agent][histone]) for histone in tail_mutants]
    Y = [abs(agent_histone_dscore[agent][histone]) for histone in fold_mutants]
    data = [X, Y]
    
    fig = plt.figure()
    plt.boxplot([X, Y], showfliers=False, notch=True, positions=[0, 0.5])
    plt.plot([random.uniform(-0.03, 0.03) for i in range(len(X))], X, '.', alpha=0.6)
    plt.plot([0.5+random.uniform(-0.03, 0.03) for i in range(len(Y))], Y, '.', alpha=0.6)
    plt.xticks([0, 0.5], ['Tail', 'Fold'])
    #plt.ylabel("Score")
    plt.ylabel("$|\Delta$ Score$|$")
    #plt.yscale('log', base=2)
    plt.title(agent_fullname[agent])
    #plt.tight_layout()
    #plt.savefig("%s_tailVSfold.png" % (agent))
    #plt.show()
    plt.close()

#sys.exit(1)

# check the WT outliers on the structure
for agent in agent_list:
    NCP = Molecules("1kx5")
    chain_resi_values_up = {}
    chain_resi_values_down = {}
    for histone in histone_list:
        if histone not in agent_outliers[agent]:
            continue
        value = agent_histone_dscore[agent][histone]
        #value = agent_histone_dChalf[agent][histone]
        #print (agent, histone, value)
        for subunit in histone_minfo[histone]:
            MW_change, pI_change, hydropathy_change = 0.0, 0.0, 0.0
            for pos in histone_minfo[histone][subunit]['mut']:
                pre_aa, post_aa = histone_minfo[histone][subunit]['mut'][pos]
                
                #if post_aa != 'A':
                #    continue
                #if pI_change >= 0:
                #    continue
                chains = subunit_chains[subunit]
                for chain in chains:
                    #if NCP.chain_seq[chain][pos-1] != pre_aa:
                    #    print (histone)

                    if value > 0:        
                        if chain not in chain_resi_values_up:
                            chain_resi_values_up[chain] = {}
                        if pos not in chain_resi_values_up[chain]:
                            chain_resi_values_up[chain][pos] = []
                        chain_resi_values_up[chain][pos].append(value)

                    else:
                        if chain not in chain_resi_values_down:
                            chain_resi_values_down[chain] = {}
                        if pos not in chain_resi_values_down[chain]:
                            chain_resi_values_down[chain][pos] = []
                        assert value < 0
                        chain_resi_values_down[chain][pos].append(value)



    #chain_resi_value = {}
    #for chain in chain_resi_values:
    #    for resi in chain_resi_values[chain]:
    #        values = chain_resi_values[chain][resi]
    #        if len(values) > 1:
    #            print (chain, resi, values)
    #        for i in range(len(values)):
    #            value = values[i]
    #            for j in range(3):
    #                new_resi = resi + i + j
    #                if new_resi not in chain_resi_values[chain]:
    #                    break
    #            if chain not in chain_resi_value:
    #                chain_resi_value[chain] = {}
    #            chain_resi_value[chain][new_resi] = value                    

    #chain_resi_mean = {}
    #for chain in chain_resi_values:
    #    for resi in chain_resi_values[chain]:
    #        if chain not in chain_resi_mean:
    #            chain_resi_mean[chain] = {}
    #        chain_resi_mean[chain][resi] = np.mean(chain_resi_values[chain][resi])

    #print (chain_resi_value)
    NCP.remove_ions()
    NCP.stylish()
    NCP.coloring(NCP.chain_resi_resn, 'white')
    NCP.make_sphere(chain_resi_values_up)
    NCP.coloring(chain_resi_values_up, 'red')
    NCP.make_sphere(chain_resi_values_down)
    NCP.coloring(chain_resi_values_down, 'blue')
    #NCP.spectrum(chain_resi_value, color_list=['white', 'yellow', 'magenta', 'red'], min=0.2, max=1)
    #NCP.spectrum(chain_resi_mean, color_list=['yellow', 'magenta', 'red'])
    #NCP.spectrum(chain_resi_mean, color_list=['blue', 'white', 'yellow', 'magenta', 'red'])
    #NCP.spectrum(chain_resi_value)
    NCP.save_session(agent)
    NCP.clear_up()
    #NCP.done()    





#sys.exit(1)
    











# Normalized count change over titration
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
agent_histone_profile = {}
agent_pick = {} # pick a titration point with max variance
for agent in agent_list:
    norm_freqs_list = []
    for num in sorted(agent_num_histone_count[agent].keys()):
        if num == 0:
            input_histone_count = agent_num_histone_count[agent][0]
            input_total = float(sum(input_histone_count.values()))
            input_freqs = [input_histone_count[histone]/input_total for histone in all_histones]
            #continue
        histone_count = agent_num_histone_count[agent][num]
        total = float(sum(histone_count.values()))
        freqs = [histone_count[histone]/total for histone in all_histones]
        norm_freqs = [float(freqs[i])/input_freqs[i] for i in range(len(freqs))]
        norm_freqs_list.append(norm_freqs)

    var_pt = sorted([[np.var(norm_freqs_list[i]), i]  for i in range(len(norm_freqs_list))], reverse=True)
    agent_pick[agent] = var_pt[0][1]
    #print (agent, agent_pick[agent])

    histone_profile = {}
    for i in range(len(all_histones)):
        histone = all_histones[i]
        profile = [ row[i] for row in norm_freqs_list ]
        histone_profile[histone] = profile

    agent_histone_profile[agent] = copy.deepcopy(histone_profile)
        
    fig = plt.figure()

    #for histone in all_histones:
    #    profile = histone_profile[histone]
    #    plt.plot(range(len(profile)), profile, '.-', alpha=0.5, label=histone)

    for cate in cate_list:
        histones = cate_histones[cate]
        color = cate_color[cate]
        for i in range(len(histones)):
            histone = histones[i]
            profile = histone_profile[histone]
            if i == 0:
                plt.plot(range(len(profile)), profile, '.-', color=color, alpha=0.5, label=cate)
            else:
                plt.plot(range(len(profile)), profile, '.-', color=color, alpha=0.5)

    plt.title(agent)
    plt.xlabel("Titration point")
    plt.ylabel("Normalized fold change")
    plt.legend()
    #plt.show()
    plt.close()


                                    
#sys.exit(1)

# PCA analysis by titration trace
agent_outliers = {}
for agent in agent_list:
    histone_profile = agent_histone_profile[agent]
    X = []
    C = []
    for i in range(len(all_histones)):
        histone = all_histones[i]
        X.append(histone_profile[histone])
        C.append(100.0 - GC_content(histone_BC[histone]))

    #scaler = preprocessing.StandardScaler().fit(X)
    #X_scaled = scaler.transform(X)

    scaler = preprocessing.StandardScaler(with_std=False).fit(X)
    X_scaled = scaler.transform(X)
    pca = PCA().fit(X_scaled)
    Xr = pca.transform(X_scaled)

    #pca = PCA().fit(X)
    #Xr = pca.transform(X)

    variance_ratio_list = 100* pca.explained_variance_ratio_

    fig = plt.figure()
    plt.plot(range(1, len(variance_ratio_list)+1), variance_ratio_list, '.-')
    plt.xlabel("PCA component")
    plt.ylabel("Variance (%)")
    plt.xticks(range(1, len(variance_ratio_list)+1))
    #plt.show()
    plt.close()

    clf = LocalOutlierFactor()
    #outcheck = clf.fit_predict(X)
    outcheck = clf.fit_predict( [row[:2] for row in Xr] )

    for i in range(len(outcheck)):
        if outcheck[i] < 0:
            outlier = all_histones[i]
            if agent not in agent_outliers:
                agent_outliers[agent] = []
            agent_outliers[agent].append(outlier)


    histone_PCA = {}
    x_list, y_list = [], []
    for i in range(len(Xr)):
        x, y = Xr[i][0], Xr[i][1]
        histone = all_histones[i]
        histone_PCA[histone] = (x, y)
        x_list.append(x)
        y_list.append(y)
        
    fig = plt.figure(figsize=(8,7))
    #plt.scatter(x_list, y_list, s=5, c=C, cmap='jet')
    for i in range(len(Xr)):
        plt.plot(Xr[i][0], Xr[i][1], 'k.')
        if outcheck[i] < 0:
            plt.plot(Xr[i][0], Xr[i][1], 'r.')
            plt.annotate(all_histones[i], (Xr[i][0], Xr[i][1]))
    plt.title("PCA plot (%s)" % (agent))
    #cbar = plt.colorbar()
    #cbar.set_label("AT content (%)")
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    #plt.tight_layout()
    #plt.show()
    plt.close()
    

    fig = plt.figure(figsize=(8,7))
    for cate in cate_list:
        histones = cate_histones[cate]
        color = cate_color[cate]
        for i in range(len(histones)):
            histone = histones[i]
            x, y = histone_PCA[histone]
            if i == 0:
                plt.plot(x, y, '.', color=color, label=cate)
            else:
                plt.plot(x, y, '.', color=color)


    plt.title("PCA plot (%s)" % (agent))
    #cbar = plt.colorbar()
    #cbar.set_label("AT content (%)")
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    #plt.tight_layout()
    #plt.show()
    plt.close()


    fig = plt.figure()
    for histone in all_histones:
        profile = histone_profile[histone]
        if histone not in agent_outliers[agent]:
            plt.plot(range(len(profile)), profile, 'k.-', alpha=0.1, label=histone)
        else:
            #plt.plot(range(len(profile)), profile, 'r.-', alpha=0.5, label=histone)
            p = plt.plot(range(len(profile)), profile, '.-', alpha=0.5, label=histone)
            plt.annotate(histone, (agent_pick[agent], profile[agent_pick[agent]]), color=p[-1].get_color())            
    
    plt.title(agent)
    plt.xlabel("Titration point")
    plt.ylabel("Normalized fold change")
    #plt.legend()
    #plt.show()
    plt.close()


# tSNE plot
for agent in agent_list:
    histone_profile = agent_histone_profile[agent]
    X = []
    for i in range(len(all_histones)):
        histone = all_histones[i]
        X.append(histone_profile[histone])

    perp=2
    tsne = sklearn.manifold.TSNE(n_components=2, perplexity=perp, init='pca', random_state=0)
    trans_data = tsne.fit_transform(X).T

    fig = plt.figure()
    for i in range(len(trans_data[0])):
        #key = num_key[i]
        #plt.plot(trans_data[0][i], trans_data[1][i], '.', color=color_list[key_cID[key]])
        plt.plot(trans_data[0][i], trans_data[1][i], '.')
    plt.title("%s tSNE plot" % (agent))
    #plt.show()
    plt.close()


# define score
agent_histone_score = {}
for agent in agent_list:
    for histone in agent_histone_profile[agent]:
        profile = agent_histone_profile[agent][histone]
        if agent == 'sp':
            profile = profile[:6] + profile[7:]
        profile = profile[1:]
        profile = np.asarray(profile)
        #profile = np.asarray(agent_histone_profile[agent][histone])
        #score = np.mean(profile)
        score = np.mean(-np.log2(profile))
        #score = profile[agent_pick[agent]]
        #score = np.log2(profile[agent_pick[agent]]+1)
        if agent not in agent_histone_score:
            agent_histone_score[agent] = {}
        agent_histone_score[agent][histone] = score

for agent in agent_list:
    origin = []
    for histone in cate_histones['WT']:
        origin.append(agent_histone_score[agent][histone])
    origin = np.mean(origin)
    for histone in all_histones:
        agent_histone_score[agent][histone] -= origin

for agent in agent_list:
    X, Y = [], []
    for histone in agent_histone_Chalf[agent]:
        Chalf = agent_histone_Chalf[agent][histone]
        score = agent_histone_score[agent][histone]
        X.append(Chalf)
        Y.append(score)

    fig = plt.figure()
    #plt.plot(X, Y, '.')
    for histone in agent_histone_Chalf[agent]:
        Chalf = agent_histone_Chalf[agent][histone]
        score = agent_histone_score[agent][histone]
        if histone not in agent_outliers[agent]:
            plt.plot(Chalf, score, 'k.')
        else:
            plt.plot(Chalf, score, 'r.')
            plt.annotate(histone, (Chalf, score))
    plt.xlabel("C half")
    plt.ylabel("Score")
    plt.title(agent)
    #plt.savefig(agent + "_scoreVSChalf.png")
    #plt.show()
    plt.close()
        
#sys.exit(1)
    
for agent in agent_list:
    score_histone = sorted([[score, histone] for histone, score in agent_histone_score[agent].items()])
    label, X, Y = [], [], []
    outlabel, outX, outY = [], [], []
    C, outC = [], []
    for i in range(len(score_histone)):
        score, histone = score_histone[i]
        if histone not in agent_outliers[agent]:
            label.append(histone)
            X.append(i)
            Y.append(score)
            C.append(100.0 - GC_content(histone_BC[histone]))
        else:
            outlabel.append(histone)
            outX.append(i)
            outY.append(score)
            outC.append(100.0 - GC_content(histone_BC[histone]))

    fig = plt.figure()
    plt.scatter(C+outC, Y+outY, s=4)
    plt.xlabel("AT content")
    plt.ylabel("Score")
    plt.title(agent)
    #plt.show()
    plt.close()

    fig = plt.figure(figsize=(18,5))
    #plt.plot(range(len(Y)), Y, 'o-')
    plt.bar(X, Y)
    plt.bar(outX, outY, color='r')
    plt.xticks(X+outX, label+outlabel, fontsize=6, rotation=80)
    for ticklabel in plt.gca().get_xticklabels():
        if ticklabel.get_text() in agent_outliers[agent]:
            ticklabel.set_color('r')
    #plt.yscale('log', base=2)
    plt.ylabel('Averaged normalized fold change')
    plt.title(agent)
    plt.tight_layout()
    #plt.show()
    plt.close()

    histone_order = {}
    for i in range(len(score_histone)):
        score, histone = score_histone[i]
        histone_order[histone] = i

    fig = plt.figure(figsize=(18,5))
    xticklocs, labels = [], []
    for cate in cate_list:
        histones = cate_histones[cate]
        color = cate_color[cate]
        X, Y = [], []
        for histone in histones:
            order = histone_order[histone]
            score = agent_histone_score[agent][histone]
            X.append(order)
            Y.append(score)
            xticklocs.append(order)
            labels.append(histone)
        plt.bar(X, Y, color=color, label=cate, align='center')
    plt.xticks(xticklocs, labels, fontsize=6, rotation=80)
    for ticklabel in plt.gca().get_xticklabels():
        cate = histone_cate[ticklabel.get_text()]
        color = cate_color[cate]
        ticklabel.set_color(color)
    plt.yscale('log', basey=2)
    plt.ylabel('Averaged normalized fold change')
    plt.title(agent)
    plt.legend()
    plt.tight_layout()
    #plt.show()
    plt.close()
    



# Clustering histones by scores
histone_state = {}
for agent in agent_list[:]:
    for histone in all_histones:
        score = agent_histone_score[agent][histone]
        if histone not in histone_state:
            histone_state[histone] = []
        histone_state[histone].append(score)


X = [histone_state[histone] for histone in all_histones]
        
#scaler = preprocessing.StandardScaler(with_std=False).fit(X)
scaler = preprocessing.StandardScaler().fit(X)
X_scaled = scaler.transform(X)
pca = PCA().fit(X_scaled)
Xr = pca.transform(X_scaled)

#pca = PCA().fit(X)
#Xr = pca.transform(X)

variance_ratio_list = 100* pca.explained_variance_ratio_

fig = plt.figure()
plt.plot(range(1, len(variance_ratio_list)+1), variance_ratio_list, '.-')
plt.xlabel("PCA component")
plt.ylabel("Variance (%)")
plt.xticks(range(1, len(variance_ratio_list)+1))
#plt.show()
plt.close()

clf = LocalOutlierFactor()
#outcheck = clf.fit_predict(X)
outcheck = clf.fit_predict( [row[:2] for row in Xr] )

histone_PCA = {}
x_list, y_list = [], []
for i in range(len(Xr)):
    histone = all_histones[i]
    x, y = Xr[i][0], Xr[i][1]
    histone_PCA[histone] = (x, y)
    x_list.append(x)
    y_list.append(y)

fig = plt.figure(figsize=(8,7))
#plt.scatter(x_list, y_list, s=5, c=C, cmap='jet')
for i in range(len(Xr)):
    plt.plot(Xr[i][0], Xr[i][1], 'k.')
    if outcheck[i] < 0:
        plt.plot(Xr[i][0], Xr[i][1], 'r.')
        plt.annotate(all_histones[i], (Xr[i][0], Xr[i][1]))
plt.title("PCA plot")
#cbar = plt.colorbar()
#cbar.set_label("AT content (%)")
plt.xlabel('PC1')
plt.ylabel('PC2')
#plt.tight_layout()
#plt.show()
plt.close()

fig = plt.figure(figsize=(8,7))
for cate in cate_list:
    histones = cate_histones[cate]
    color = cate_color[cate]
    marker = cate_marker[cate]
    if cate == 'WT+mut':
        alpha=0.5
    else:
        alpha=1.0
    for i in range(len(histones)):
        histone = histones[i]
        x, y = histone_PCA[histone]
        if i == 0:
            if cate == 'WT':
                plt.plot(x, y, '.', color=color, marker=marker, alpha=alpha, label=cate, zorder=100)
            else:
                plt.plot(x, y, '.', color=color, marker=marker, alpha=alpha, label=cate)
        else:
            if cate == 'WT':
                plt.plot(x, y, '.', color=color, marker=marker, alpha=alpha, zorder=100)
            else:
                plt.plot(x, y, '.', color=color, marker=marker, alpha=alpha)

plt.title("PCA plot")
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.legend()
#plt.show()
plt.close()




#sys.exit(1)

# Oncohistone free DNA contamination check
histone_freeDNA = read_old("Onco_freeDNA.csv")

all_histones = list(set(all_histones) - set(['H3.1 K4M']))

X = []
for histone in all_histones:
    freeDNA = histone_freeDNA[histone]
    X.append(freeDNA*100)

for agent in agent_list:
    Y = []
    for histone in all_histones:
        #score = agent_histone_score[agent][histone]
        #Y.append(score)
        Chalf = agent_histone_Chalf[agent][histone]
        Y.append(Chalf)

    corr = scipy.stats.pearsonr(X, Y)[0]
    #print ("%s VS %s: %1.2f" % ('FreeDNA', agent, corr))


    fig = plt.figure()
    #plt.plot(X, Y, '.')
    for histone in all_histones:
    #for histone in cate_histones['WT'] + cate_histones['WT+mut']:
        x = histone_freeDNA[histone]*100
        y = agent_histone_Chalf[agent][histone]
        #y = agent_histone_score[agent][histone]
        if histone not in agent_outliers[agent]:
            plt.plot(x, y, 'k.')
        else:
            plt.plot(x, y, 'r.')
            plt.annotate(histone, (x, y))
    plt.title(agent)
    plt.xlabel("Free DNA (%)")
    plt.ylabel("score")
    #plt.xlim([-10, 30])
    #plt.ylim([0, 2.5])
    #plt.yscale('log', base=2)
    #plt.xscale('log', base=2)
    #plt.show()
    plt.close()

#sys.exit(1)


#### only WT and WT+mut
# Clustering histones by titration profile
all_WT_histones = cate_histones['WT'] + cate_histones['WT+mut']
agent_outliers = {}
for agent in agent_list:
    histone_profile = agent_histone_profile[agent]
    X = []
    C = []
    for i in range(len(all_WT_histones)):
        histone = all_WT_histones[i]
        X.append(histone_profile[histone])
        C.append(100.0 - GC_content(histone_BC[histone]))

    #scaler = preprocessing.StandardScaler().fit(X)
    #X_scaled = scaler.transform(X)

    #scaler = preprocessing.StandardScaler().fit(X)
    scaler = preprocessing.StandardScaler(with_std=False).fit(X)
    X_scaled = scaler.transform(X)
    pca = PCA().fit(X_scaled)
    Xr = pca.transform(X_scaled)

    #pca = PCA().fit(X)
    #Xr = pca.transform(X)

    variance_ratio_list = 100* pca.explained_variance_ratio_

    fig = plt.figure()
    plt.plot(range(1, len(variance_ratio_list)+1), variance_ratio_list, '.-')
    plt.xlabel("PCA component")
    plt.ylabel("Variance (%)")
    plt.xticks(range(1, len(variance_ratio_list)+1))
    #plt.show()
    plt.close()

    clf = LocalOutlierFactor()
    #outcheck = clf.fit_predict(X)
    outcheck = clf.fit_predict([row[:2] for row in Xr])
    #outcheck = clf.fit_predict([row[:3] for row in Xr])

    for i in range(len(outcheck)):
        if outcheck[i] < 0:
            outlier = all_WT_histones[i]
            if agent not in agent_outliers:
                agent_outliers[agent] = []
            agent_outliers[agent].append(outlier)


    histone_PCA = {}
    x_list, y_list = [], []
    for i in range(len(Xr)):
        x, y = Xr[i][0], Xr[i][1]
        histone = all_WT_histones[i]
        histone_PCA[histone] = (x, y)
        x_list.append(x)
        y_list.append(y)
        
    fig = plt.figure(figsize=(8,7))
    #plt.scatter(x_list, y_list, s=5, c=C, cmap='jet')
    for i in range(len(Xr)):
        plt.plot(Xr[i][0], Xr[i][1], 'k.')
        if outcheck[i] < 0:
            plt.plot(Xr[i][0], Xr[i][1], 'r.')
            plt.annotate(all_WT_histones[i], (Xr[i][0], Xr[i][1]))
    plt.title("PCA plot (%s)" % (agent))
    #cbar = plt.colorbar()
    #cbar.set_label("AT content (%)")
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    #plt.tight_layout()
    #plt.show()
    plt.close()

    fig = plt.figure()
    for histone in all_WT_histones:
        profile = histone_profile[histone]
        if histone not in agent_outliers[agent]:
            plt.plot(range(len(profile)), profile, 'k.-', alpha=0.1, label=histone)
        else:
            #plt.plot(range(len(profile)), profile, 'r.-', alpha=0.5, label=histone)
            p = plt.plot(range(len(profile)), profile, '.-', alpha=0.5, label=histone)
            plt.annotate(histone, (agent_pick[agent], profile[agent_pick[agent]]), color=p[-1].get_color())            
    
    plt.title(agent)
    plt.xlabel("Titration point")
    plt.ylabel("Normalized fold change")
    #plt.legend()
    #plt.show()
    plt.close()


for agent in agent_list:
    score_histone = sorted([[agent_histone_score[agent][histone], histone] for histone in all_histones])
    label, X, Y = [], [], []
    outlabel, outX, outY = [], [], []
    C, outC = [], []
    for i in range(len(score_histone)):
        score, histone = score_histone[i]
        if histone not in agent_outliers[agent]:
            label.append(histone)
            X.append(i)
            Y.append(score)
            C.append(100.0 - GC_content(histone_BC[histone]))
        else:
            outlabel.append(histone)
            outX.append(i)
            outY.append(score)
            outC.append(100.0 - GC_content(histone_BC[histone]))

    fig = plt.figure()
    plt.scatter(C+outC, Y+outY, s=4)
    plt.xlabel("AT content")
    plt.ylabel("Score")
    plt.title(agent)
    #plt.show()
    plt.close()

    fig = plt.figure(figsize=(18,5))
    #plt.plot(range(len(Y)), Y, 'o-')
    plt.bar(X, Y)
    plt.bar(outX, outY, color='r')
    plt.xticks(X+outX, label+outlabel, fontsize=6, rotation=80)
    for ticklabel in plt.gca().get_xticklabels():
        if ticklabel.get_text() in agent_outliers[agent]:
            ticklabel.set_color('r')
    plt.yscale('log', base=2)
    plt.ylabel('Averaged normalized fold change')
    plt.title(agent)
    plt.tight_layout()
    #plt.show()
    plt.close()


# clustering by score
#agent_outliers = {}

all_histones = cate_histones['WT'] + cate_histones['WT+mut']
X = [histone_state[histone] for histone in all_histones]
        
#scaler = preprocessing.StandardScaler(with_std=False).fit(X)
scaler = preprocessing.StandardScaler().fit(X)
X_scaled = scaler.transform(X)
pca = PCA().fit(X_scaled)
Xr = pca.transform(X_scaled)

#pca = PCA().fit(X)
#Xr = pca.transform(X)

variance_ratio_list = 100* pca.explained_variance_ratio_

fig = plt.figure()
plt.plot(range(1, len(variance_ratio_list)+1), variance_ratio_list, '.-')
plt.xlabel("PCA component")
plt.ylabel("Variance (%)")
plt.xticks(range(1, len(variance_ratio_list)+1))
#plt.show()
plt.close()

clf = LocalOutlierFactor()
#outcheck = clf.fit_predict(X)
outcheck = clf.fit_predict( [row[:2] for row in Xr] )

#for i in range(len(outcheck)):
#    if outcheck[i] < 0:
#        for agent in agent_list:
#            if agent not in agent_outliers:
#                agent_outliers[agent] = []
#            agent_outliers[agent].append(all_histones[i])

histone_PCA = {}
x_list, y_list = [], []
for i in range(len(Xr)):
    histone = all_histones[i]
    x, y = Xr[i][0], Xr[i][1]
    histone_PCA[histone] = (x, y)
    x_list.append(x)
    y_list.append(y)

fig = plt.figure(figsize=(8,7))
#plt.scatter(x_list, y_list, s=5, c=C, cmap='jet')
for i in range(len(Xr)):
    plt.plot(Xr[i][0], Xr[i][1], 'k.')
    if outcheck[i] < 0:
        plt.plot(Xr[i][0], Xr[i][1], 'r.')
        plt.annotate(all_histones[i], (Xr[i][0], Xr[i][1]))
plt.title("PCA plot")
#cbar = plt.colorbar()
#cbar.set_label("AT content (%)")
plt.xlabel('PC1')
plt.ylabel('PC2')
#plt.tight_layout()
#plt.show()
plt.close()

for agent in agent_list:
    score_histone = sorted([[agent_histone_score[agent][histone], histone] for histone in all_histones])
    label, X, Y = [], [], []
    outlabel, outX, outY = [], [], []
    C, outC = [], []
    for i in range(len(score_histone)):
        score, histone = score_histone[i]
        if histone not in agent_outliers[agent]:
            label.append(histone)
            X.append(i)
            Y.append(score)
            C.append(100.0 - GC_content(histone_BC[histone]))
        else:
            outlabel.append(histone)
            outX.append(i)
            outY.append(score)
            outC.append(100.0 - GC_content(histone_BC[histone]))

    fig = plt.figure()
    plt.scatter(C+outC, Y+outY, s=4)
    plt.xlabel("AT content")
    plt.ylabel("Score")
    plt.title(agent)
    #plt.show()
    plt.close()

    fig = plt.figure(figsize=(18,5))
    #plt.plot(range(len(Y)), Y, 'o-')
    plt.bar(X, Y)
    plt.bar(outX, outY, color='r')
    plt.xticks(X+outX, label+outlabel, fontsize=6, rotation=80)
    for ticklabel in plt.gca().get_xticklabels():
        if ticklabel.get_text() in agent_outliers[agent]:
            ticklabel.set_color('r')
    plt.yscale('log', basey=2)
    plt.ylabel('Averaged normalized fold change')
    plt.title(agent)
    plt.tight_layout()
    #plt.show()
    plt.close()


# property analysis
property_list = ['MW', 'pI', 'hydropathy']

agent_histone_property_change = {}

fig1, axes1 = plt.subplots(nrows=len(agent_list), ncols=len(property_list))
fig2, axes2 = plt.subplots(nrows=len(agent_list), ncols=len(property_list))
for i in range(len(agent_list)):
    agent = agent_list[i]
    for histone in cate_histones["WT+mut"]:
        MW_change, pI_change, hydropathy_change = 0.0, 0.0, 0.0
        for subunit in histone_minfo[histone]:
            for pos in histone_minfo[histone][subunit]['mut']:
                pre_aa, post_aa = histone_minfo[histone][subunit]['mut'][pos]
                MW_change += aa_info[post_aa]['MW'] - aa_info[pre_aa]['MW']
                pI_change += aa_info[post_aa]['pI'] - aa_info[pre_aa]['pI']
                hydropathy_change += aa_info[post_aa]['KD_hydropathy'] \
                    - aa_info[pre_aa]['KD_hydropathy']

        if agent not in agent_histone_property_change:
            agent_histone_property_change[agent] = {}
            
        assert histone not in agent_histone_property_change[agent]
        agent_histone_property_change[agent][histone] = {}

        agent_histone_property_change[agent][histone]['MW'] = MW_change
        agent_histone_property_change[agent][histone]['pI'] = pI_change
        agent_histone_property_change[agent][histone]['hydropathy'] = hydropathy_change

    for j in range(len(property_list)):
        prop = property_list[j]
        X, Y = [], []
        for histone in agent_histone_property_change[agent]:
            X.append(agent_histone_property_change[agent][histone][prop])
            Y.append(agent_histone_score[agent][histone])

        axes1[i, j].plot(X, Y, '.')
        axes1[i, j].set_yscale('log', basey=2)

        if j <= 0:
            axes1[i, j].set_ylabel(agent, rotation=0, labelpad=20)
        else:
            axes1[i, j].tick_params(axis='y', which='both', labelleft=False)

        if i >= len(agent_list) - 1:
            axes1[i, j].set_xlabel('%s change' % (prop))
        else:
            axes1[i,j].tick_params(axis='x', which='both', labelbottom=False)


        #corr = scipy.stats.spearmanr(X, Y)[0]
        corr = scipy.stats.pearsonr(X, Y)[0]
        #print ("%s\t%s VS %s: %1.2f" % (agent, prop, "score", corr))

        matrix = np.zeros((len(agent_list), len(agent_list)))
        matrix[:] = corr
        img = axes2[i, j].imshow(matrix, cmap="bwr", vmin=-0.5, vmax=0.5, origin='lower')

        if corr < -0.35 or corr > 0.35:
            color = 'white'
        else:
            color = 'black'

        axes2[i,j].text(len(agent_list)/2, len(agent_list)/2, str(round(corr,2)), ha="center", va="center", fontsize=10, color=color, weight='bold')
        axes2[i,j].set_xlim([0, len(agent_list)-1])
        axes2[i,j].set_ylim([0, len(agent_list)-1])
        axes2[i,j].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labeltop=False, labelleft=False, labelright=False)

        if j <= 0:
            axes2[i, j].set_ylabel(agent, rotation=0, labelpad=20)
        else:
            axes2[i, j].tick_params(axis='y', which='both', labelleft=False)

        if i >= len(agent_list) - 1:
            axes2[i, j].set_xlabel('%s change' % (prop))
        else:
            axes2[i,j].tick_params(axis='x', which='both', labelbottom=False)



fig1.suptitle("Corrrelation betwen property change")
#plt.show()
#fig1.show(fig1)
#plt.close(fig1)

fig2.suptitle("Corrrelation betwen property change")
cbar=fig2.colorbar(img, ax=axes2, location='right', shrink=0.8)
cbar.ax.set_ylabel('Spearman correlation', rotation=-90, va="bottom")
#fig2.show()
#plt.show()
#plt.close()
#plt.close(fig2)
plt.close('all')


# categorize mutations
sub_prop_changes = {}
for histone in cate_histones["WT+mut"]:
    for subunit in histone_minfo[histone]:
        for pos in histone_minfo[histone][subunit]['mut']:
            pre_aa, post_aa = histone_minfo[histone][subunit]['mut'][pos]

            sub = '%s -> %s' % (pre_aa, post_aa)
            if sub not in sub_prop_changes:
                sub_prop_changes[sub] = {}
            else:
                continue
            
            MW_change = aa_info[post_aa]['MW'] - aa_info[pre_aa]['MW']
            pI_change = aa_info[post_aa]['pI'] - aa_info[pre_aa]['pI']
            hydropathy_change = aa_info[post_aa]['KD_hydropathy'] \
                    - aa_info[pre_aa]['KD_hydropathy']
            
            sub_prop_changes[sub]['MW'] = MW_change
            sub_prop_changes[sub]['pI'] = pI_change
            sub_prop_changes[sub]['hydropathy'] = hydropathy_change

sub_list = list(sub_prop_changes.keys())
X, Y, = [], []
C = []
state = []
for sub in sub_list:
    MW_change = sub_prop_changes[sub]['MW']
    pI_change = sub_prop_changes[sub]['pI']
    hydropathy_change = sub_prop_changes[sub]['hydropathy']
    X.append(MW_change)
    Y.append(pI_change)
    C.append(hydropathy_change)
    state.append([MW_change, pI_change, hydropathy_change])

fig = plt.figure()
plt.scatter(X, Y, s=4, c=C)
for i in range(len(sub_list)):
    sub = sub_list[i]
    plt.annotate(sub_list[i], (X[i], Y[i]))
plt.xlabel("MW change")
plt.ylabel("pI chnage")
cbar = plt.colorbar()
cbar.ax.set_ylabel('Hydropathy change', rotation=-90, va="bottom")
plt.title("Substitution mutations")
#plt.show()
plt.close()


# PCA by physical information of mutation (property change, position)
NCP = Molecules("1kx5")
subunit_chains = {'H2A':['C'], 'H2B':['D'], 'H3.1':['A'], 'H4':['B']}

histone_coord = {}
for histone in cate_histones["WT+mut"]:
    coord = [0, 0, 0]
    count = 0
    for subunit in histone_minfo[histone]:
        for pos in histone_minfo[histone][subunit]['mut']:
            chain = subunit_chains[subunit][0]
            for index in NCP.chain_resi_index_atom[chain][pos]:
                if NCP.chain_resi_index_atom[chain][pos][index]['name'] == 'CA':
                    x, y, z = NCP.chain_resi_index_atom[chain][pos][index]['coord']
                    coord[0] += x; coord[1] += y; coord[2] += z
                    count +=1
                    break
    histone_coord[histone] = [float(comp)/count for comp in coord]
    
histone_state = {}
for histone in cate_histones["WT+mut"]:
    MW_change = agent_histone_property_change['sp'][histone]['MW']
    pI_change = agent_histone_property_change['sp'][histone]['pI']
    hydropathy_change = agent_histone_property_change['sp'][histone]['hydropathy']
    coord = histone_coord[histone]
    state = [MW_change, pI_change, hydropathy_change] + coord
    histone_state[histone] = state


all_histones = cate_histones['WT+mut']
X = [histone_state[histone] for histone in all_histones]

#scaler = preprocessing.StandardScaler(with_std=False).fit(X)
scaler = preprocessing.StandardScaler().fit(X)
X_scaled = scaler.transform(X)
pca = PCA().fit(X_scaled)
Xr = pca.transform(X_scaled)

#pca = PCA().fit(X)
#Xr = pca.transform(X)

variance_ratio_list = 100* pca.explained_variance_ratio_

fig = plt.figure()
plt.plot(range(1, len(variance_ratio_list)+1), variance_ratio_list, '.-')
plt.xlabel("PCA component")
plt.ylabel("Variance (%)")
plt.xticks(range(1, len(variance_ratio_list)+1))
#plt.show()
plt.close()


for agent in agent_list[0:0]:
    histone_PCA = {}
    x_list, y_list, z_list = [], [], []
    C = []
    for i in range(len(Xr)):
        histone = all_histones[i]
        x, y, z = Xr[i][0], Xr[i][1], Xr[i][2]
        histone_PCA[histone] = (x, y)
        x_list.append(x)
        y_list.append(y)
        z_list.append(z)
        C.append(np.log2(1+agent_histone_score[agent][histone]))

    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    img = ax.scatter(x_list, y_list, z_list, c=C)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_zlabel('PC3')
    plt.title("PCA plot(%s)" % (agent))
    cbar = plt.colorbar(img)
    cbar.set_label("Score")
    plt.show()
    plt.close()

    fig = plt.figure(figsize=(8,7))
    plt.scatter(x_list, y_list, s=5, c=C, cmap='jet')
    #for i in range(len(Xr)):
    #    plt.annotate(all_histones[i], (Xr[i][0], Xr[i][1]))
    plt.title("PCA plot(%s)" % (agent))
    cbar = plt.colorbar()
    cbar.set_label("Score")
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    #plt.tight_layout()
    #plt.show()
    plt.close()


#perp=2
#tsne = sklearn.manifold.TSNE(n_components=2, perplexity=perp, init='pca', random_state=0)
tsne = sklearn.manifold.TSNE()
trans_data = tsne.fit_transform(X).T

for agent in agent_list:
    C = []
    for i in range(len(all_histones)):
        histone  = all_histones[i]
        C.append(np.log2(1+agent_histone_score[agent][histone]))

    fig = plt.figure(figsize=(8,7))
    plt.scatter(trans_data[0], trans_data[1], s=5, c=C, cmap='jet')
    plt.title("tSNE plot")
    cbar = plt.colorbar()
    cbar.set_label("%s np.log2(1+score)" % (agent))
    #plt.xlabel('PC1')
    #plt.ylabel('PC2')
    #plt.tight_layout()
    plt.show()
    plt.close()



#sys.exit(1)
