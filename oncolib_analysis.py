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
from pymol_graphics import Molecules
import random
from scipy.optimize import curve_fit
from sklearn import linear_model
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


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
    fig = plt.figure()
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

        
    plt.title(agent_fullname[agent])
    plt.xlabel("Concentration")
    plt.ylabel("Soluble fraction")
    if agent in ['HP1a']:
        plt.xscale('log', base=2)
    elif agent in ['sp', 'spd', 'CoH']:
        plt.xscale('log', base=10)
    #plt.legend()
    #plt.show()
    #plt.savefig(agent+'.png')
    plt.close()

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


        fig = plt.figure()
        for histone in histone_list:
            dChalf = agent_histone_dChalf[agent][histone]
            dscore = agent_histone_dscore[agent][histone]
            if histone not in agent_outliers[agent]:
                plt.plot(dChalf, dscore, 'k.')
            else:
                plt.plot(dChalf, dscore, 'r.')
                plt.annotate(histone, (dChalf, dscore))

        plt.axvline(x=0, linestyle='--', color='k', alpha=0.5)
        plt.axhline(y=0, linestyle='--', color='k', alpha=0.5)
        plt.title(agent_fullname[agent])
        plt.xlabel("$\Delta$ C-half")
        plt.ylabel("$\Delta$ Score")
        #plt.savefig(str(k) + '_' + agent + "_ChalfVSScore.png")
        #plt.show()
        plt.close()

    agent_outliers_list.append(copy.deepcopy(agent_outliers))
all_agent_outliers, WTmut_agent_outliers = agent_outliers_list

#sys.exit(1)


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
        plt.savefig("%s_%s_bar.png" % (agent, ylabel))
        #plt.show()
        plt.close()

sys.exit(1)

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

cate_list = ['WT', 'WT+mut', 'H2A.Z', 'H2A.Z+mut', 'H2A.V', 'H2A.V+mut', 'H3.3', 'H3.3+mut', 'tag', 'freeDNA']
histone_list = all_good_histones

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
        if histone in outliers:
            plt.annotate(histone, (x,y))
                

plt.title("PCA plot")
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.legend()
#plt.savefig("PCA_allagent.png")
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
            plt.yscale('log', base=2)
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


fig, axes = plt.subplots(nrows=len(agent_list), ncols=len(agent_list))

for i in range(len(agent_list)):
    for j in range(len(agent_list)):
        idx = len(agent_list)*i + j
        agent1, agent2 = agent_list[i], agent_list[j]
        if i > j:
            X = [ agent_histone_score[agent1][histone] for histone in histone_list ]
            Y = [ agent_histone_score[agent2][histone] for histone in histone_list ]
            axes[i,j].plot(X, Y, 'k.', markersize=1)
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
            axes[i,j].text(len(agent_list)/2, len(agent_list)/2, s, ha="center", va="center", fontsize=10, weight='bold')
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
            axes[i,j].text(len(agent_list)/2, len(agent_list)/2, str(round(value,2)), ha="center", va="center", fontsize=10, color=color, weight='bold')
            axes[i,j].set_xlim([0, len(agent_list)-1])
            axes[i,j].set_ylim([0, len(agent_list)-1])
            axes[i,j].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labeltop=False, labelleft=False, labelright=False)

plt.subplots_adjust(wspace=0.1, hspace=0.1)
cbar=fig.colorbar(img, ax=axes, location='right', shrink=0.8)
cbar.ax.set_ylabel('Spearman correlation', rotation=-90, va="bottom")
plt.suptitle("Corrrelation betwen condensing agents")
#plt.show()
plt.close()

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
    plt.savefig("%s_tailVSfold.png" % (agent))
    #plt.show()
    plt.close()

sys.exit(1)

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
        print (agent, histone, value)
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





sys.exit(1)
    











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
        
sys.exit(1)
    
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
    plt.yscale('log', base=2)
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
    print ("%s VS %s: %1.2f" % ('FreeDNA', agent, corr))


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

sys.exit(1)

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
    plt.yscale('log', base=2)
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
        axes1[i, j].set_yscale('log', base=2)

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
