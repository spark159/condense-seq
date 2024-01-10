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
import matplotlib
import pickle
from scipy.stats import norm

def GC_content(seq):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    return (num/float(len(seq)))*100

def rescale (value_list, old_st, old_ed, new_st, new_ed):
    output = []
    for value in value_list:
        assert value >= old_st and value <= old_ed
        new_value = new_st + (new_ed - new_st)*float(value-old_st)/(old_ed-old_st)
        output.append(new_value)
    return output

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


# PTM libraray NGS information
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
index_name, name_index, name_titr = read_index("PTMlib_NGS_information.csv")


# read PTM library information
# parsing the histone mutation information
def mhistone_parser (mhistone):
    hname, mutations = re.split('(H2A(?:[.]\w)?|H2B|H3(?:[.]\d)?|H4)', mhistone)[1:]
    #all_mutations = set([])
    pos_mutation = {}
    pattern = '([A-Z])(\d+(?:,\d+)*)(ac|me2[as]|me[1-3]|me|ub|ph|cr|GlcNAc|[A-Z])'  
    for find in re.findall(pattern, mutations):
        aa, pos_list, mutation = find
        assert aa in aa_info.keys()
        if mutation in aa_info.keys():
            mtype = 'mut'
        elif mutation.startswith('me'):
            mtype = 'me'
        else:
            mtype = mutation
        for pos in pos_list.split(','):
            pos = int(pos)
            assert pos not in pos_mutation
            pos_mutation[pos] = (mtype, aa, mutation)
        #all_mutations.add((hname, pos, mtype, aa, mutation))
    #return hname, all_mutations, pos_mutation
    return hname, pos_mutation

def read_table (fname):
    subunit_list = ['H2A', 'H2B', 'H3', 'H4']
    ID_BC, BC_ID = {}, {}
    ID_minfo = {}
    #ID_mutations = {}
    ID_shortname = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if First:
            First = False
            continue
        if not line:
            continue
        cols = line.split('\t')
        ID, H2A, H2B, H3, H4, DNA, BC = cols
        ID = int(ID)
        assert ID not in ID_BC
        ID_BC[ID] = BC
        assert BC not in BC_ID
        BC_ID[BC] = ID

        assert ID not in ID_minfo
        ID_minfo[ID] = {}
        mutations = []
        shortname = ""
        for subunit, mhistone in zip(subunit_list, [H2A, H2B, H3, H4]):
            mhistone = mhistone.strip()
            ID_minfo[ID][subunit] = {}
            if mhistone == 'NA':
                ID_minfo[ID][subunit]['name'] = None
                ID_minfo[ID][subunit]['mutations'] = {}
            else:
                hname, pos_mutation = mhistone_parser(mhistone)
                ID_minfo[ID][subunit]['name'] = hname
                ID_minfo[ID][subunit]['mutations'] = pos_mutation

            if 'KpolyAc' in mhistone:
                mutation = subunit + 'KpolyAc'
            elif 'Acidic Patch Mutant' in mhistone:
                mutation = subunit + ' AP mutant'
            else:
                mutation = mhistone
            if mutation != subunit:
                mutations.append(mutation)

        if set(mutations) == set(['NA']):
            shortname += 'freeDNA'
        elif len(mutations) == 0:
            shortname += 'WT'
        else:
            shortname += '/'.join(mutations)

        if DNA == 'NA':
            ID_minfo[ID]['DNA'] = None
        else:
            ID_minfo[ID]['DNA'] = DNA
            shortname += ' (' + DNA + ')'

        assert ID not in ID_shortname
        ID_shortname[ID] = shortname

    return ID_BC, BC_ID, ID_minfo, ID_shortname
ID_BC, BC_ID, ID_minfo, ID_shortname = read_table('PTMlibTable.csv')


# remove IDs already susceptible to PstI digest
# BCs with PstI site (ID 71, 78, 113) and free DNA with PstI site (ID 117)
all_IDs = sorted(ID_BC.keys())
all_good_IDs = sorted(list(set(all_IDs) - set([71, 78, 113, 117])))
agent_fullname = {'sp':'Spermine(4+)', 'spd':'Spermidine(3+)', 'CoH':'Cobalt Hexammine(3+)', 'PEG':'PEG 8000', 'HP1a':'HP1 $\\alpha$', 'sp_filter':'Spermine(4+) 100kD filter', 'spd_filter':'Spermidine(3+) 100kD filter'}


# categorize IDs
#cate_list = ['freeDNA', 'WT', 'WT+CpGme', 'WT+Nmut', 'WT+NPTM', 'Var', 'Var+Nmut']
def categorize_IDs (IDs):
    subunit_list = ['H2A', 'H2B', 'H3', 'H4']
    cate_IDs, ID_cate = {}, {}
    for ID in IDs:
        minfo = ID_minfo[ID]
        hnames, mtypes = [], []
        for subunit in subunit_list:
            hname = minfo[subunit]['name']
            hnames.append(hname)
            for pos in minfo[subunit]['mutations']:
                mtype, _, _ = minfo[subunit]['mutations'][pos]
                mtypes.append(mtype)

        cate = ""
        if set(hnames) == set([None]):
            cate += 'freeDNA'
        elif set(hnames) == set(['H2A', 'H2B', 'H3', 'H4']):
            cate += 'WT'
        else:
            #print hnames
            assert len(hnames) == 4
            varnames = sorted(list(set(hnames) - set(['H2A', 'H2B', 'H3', 'H4'])))
            cate += '/'.join(varnames)
            #cate += 'Var'

        if len(mtypes) <= 0:
            pass
        elif set(mtypes) == set(['mut']):
            cate += '+' + str(len(mtypes)) + 'mut'
        else:
            assert 'mut' not in mtypes
            cate += '+' + str(len(mtypes)) + 'PTM'

        if minfo['DNA'] == 'CpGme':
            cate += '+CpGme'

        if cate not in cate_IDs:
            cate_IDs[cate] = []
        cate_IDs[cate].append(ID)
        ID_cate[ID] = cate

    return cate_IDs, ID_cate
cate_IDs, ID_cate = categorize_IDs(all_good_IDs)


# grouping IDs
#group_list = ['H2A/2Bac','H3ac','H4ac','H3ac+H4ac','H3me','H4me','H3me+H4ac','H3ph','+ub','WT','WT+mut','Var','Var+mut','WT+CpGme','GlcNAc']
def grouping_IDs (IDs):
    subunit_list = ['H2A', 'H2B', 'H3', 'H4']
    group_IDs, ID_group = {}, {}
    for ID in IDs:
        terms = []
        minfo = ID_minfo[ID]
        hnames, mtypes = [], []
        for subunit in subunit_list:
            hname = minfo[subunit]['name']
            hnames.append(hname)
            for pos in minfo[subunit]['mutations']:
                mtype, _, _ = minfo[subunit]['mutations'][pos]
                mtypes.append(mtype)
                terms.append((subunit, mtype))

        if len(mtypes) <= 0:
            if set(hnames) == set([None]):
                group = 'freeDNA'
            elif set(hnames) == set(['H2A', 'H2B', 'H3', 'H4']):
                group = 'WT'
                if minfo['DNA'] == 'CpGme':
                    group += '+CpGme'
            else:
                group = 'Var'
            
        elif set(mtypes) == set(['mut']):
            if set(hnames) == set(['H2A', 'H2B', 'H3', 'H4']):
                if 'AP' in ID_shortname[ID]:
                    group = 'AP mutant'
                else:
                    group = 'WT+mut'
            else:
                group = 'Var+mut'
        else:
            if 'ub' in mtypes:
                group = '+ub'
            elif 'KpolyAc' in ID_shortname[ID]:
                group = 'KpolyAc'
            elif 'GlcNAc' in mtypes:
                group = 'GlcNAc'
            else:
                newterms = []
                for subunit, mtype in terms:
                    newterm = ""
                    if subunit in ['H2A', 'H2B']:
                        newterm += 'H2A/B'
                    else:
                        newterm += subunit
                    newterm += mtype
                    newterms.append(newterm)
                newterms = list(set(newterms))
                group = '+'.join(newterms)
                    
        if group not in group_IDs:
            group_IDs[group] = []
        group_IDs[group].append(ID)
        ID_group[ID] = group

    return group_IDs, ID_group
group_IDs, ID_group = grouping_IDs(all_good_IDs)

#sys.exit(1)


# read sort file
def read_sort (fname):
    validity_type_count = {}
    name_count = {}
    agent_num_ID_count = {}
    for line in open(fname):
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
                ID = BC_ID[BC]
            except:
                continue

            if name not in name_count:
                name_count[name] = 0
            name_count[name] +=1

            agent = name[:-1]
            num = int(name[-1])
            if agent not in agent_num_ID_count:
                agent_num_ID_count[agent] = {}
            if num not in agent_num_ID_count[agent]:
                agent_num_ID_count[agent][num] = {}
            if ID not in agent_num_ID_count[agent][num]:
                agent_num_ID_count[agent][num][ID] = 0
            agent_num_ID_count[agent][num][ID] +=1
    return validity_type_count, name_count, agent_num_ID_count
sort_fname1 = "Sp-Spd-CoHex-PEG-HP1a-PTMlib_S1_L001_R1_001.sort"
sort_fname2 = "Sp-Spd-PTMlib-100kdfilter_S1_L001_R1_001.sort"

validity_type_count1, name_count1, agent_num_ID_count1 = read_sort(sort_fname1)
validity_type_count2, name_count2, agent_num_ID_count2 = read_sort(sort_fname2)

agent_num_ID_count = {}
agent_num_ID_count.update(agent_num_ID_count1)
agent_num_ID_count['sp_filter'] = agent_num_ID_count2['sp']
agent_num_ID_count['spd_filter'] = agent_num_ID_count2['spd']

# sequencing QC
if False:
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
    all_IDs = sorted(ID_BC.keys())
    #agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
    agent_list = ['sp', 'spd']
    fig = plt.figure()
    for agent in agent_list:
        ID_count = agent_num_ID_count[agent][0]
        X, Y = [], []
        for ID in ID_count:
        #for ID in all_IDs:
            X.append(ID)
            Y.append(ID_count[ID])
        plt.plot(X, Y, '.-', alpha=0.8, label=agent)
    #pstI_IDs = [71, 78, 113, 117]
    #for x in pstI_IDs:
    #    plt.axvline(x=x, color='black', linestyle='--', alpha=0.8, zorder=1)
    #xtick_list = [0, 20, 40, 60, 80, 100, 120]
    #plt.xticks(xtick_list + pstI_IDs, [str(xtick) for xtick in xtick_list] + ['\n'+str(ID) for ID in pstI_IDs], rotation=10)
    plt.xlabel('PTM library IDs')
    plt.ylabel('Read counts')
    plt.title('Input read counts')
    plt.legend()
    #plt.show()
    plt.close()

#sys.exit(1)

# get fold change over titrations
#agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
#agent_list = ['sp_filter', 'spd_filter']
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a', 'sp_filter', 'spd_filter']
agent_num_ID_fraction = {}
for agent in agent_list:
    for num in agent_num_ID_count[agent]:
        total = sum(agent_num_ID_count[agent][num].values())
        ID_fraction = {}
        for ID in agent_num_ID_count[agent][num]:
            count = agent_num_ID_count[agent][num][ID]
            fraction = float(count) / total
            ID_fraction[ID] = fraction
        if agent not in agent_num_ID_fraction:
            agent_num_ID_fraction[agent] = {}
        agent_num_ID_fraction[agent][num] = copy.deepcopy(ID_fraction)

agent_ID_profile = {}
for agent in agent_num_ID_fraction:
    for num in sorted(agent_num_ID_fraction[agent]):
        for ID in agent_num_ID_fraction[agent][num]:
            try:
                control = agent_num_ID_fraction[agent][0][ID]
            except:
                control = 1.0 # temporal
            fold_change = agent_num_ID_fraction[agent][num][ID] / float(control)
            if agent not in agent_ID_profile:
                agent_ID_profile[agent] = {}
            if ID not in agent_ID_profile[agent]:
                agent_ID_profile[agent][ID] = []
            agent_ID_profile[agent][ID].append(fold_change)

# plot fold change over titrations
for agent in agent_list:
    fig = plt.figure()
    for ID in all_good_IDs:
        profile = agent_ID_profile[agent][ID]
        p = plt.plot(range(len(profile)), profile, '.-', alpha=0.5, label=ID)
        plt.annotate(str(ID), (len(profile)-1, profile[-1]), color=p[-1].get_color())
    plt.xlabel("Titration point")
    plt.ylabel("Normalized fold change")
    plt.title(agent)
    #plt.legend()
    #plt.show()
    plt.close()


# get survival probability over titrations
# read titration data
def read_titration_data (fname):
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

#agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
#agent_filename = {'sp':"PTMlib_spermine_corrected.csv", 'spd':"PTMlib_spermidine_corrected.csv", 'CoH':"PTMlib_CoHex_corrected.csv", 'PEG':"PTMlib_PEG.csv", 'HP1a':"PTMlib_HP1a.csv"}
#agent_list = ['sp_filter', 'spd_filter']
#agent_filename = {'sp_filter':"PTMlib_spermine_filter_corrected.csv", 'spd_filter':"PTMlib_spermidine_filter_corrected.csv"}
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a', 'sp_filter', 'spd_filter']
agent_filename = {'sp':"PTMlib_spermine_corrected.csv", 'spd':"PTMlib_spermidine_corrected.csv", 'CoH':"PTMlib_CoHex_corrected.csv", 'PEG':"PTMlib_PEG.csv", 'HP1a':"PTMlib_HP1a.csv", 'sp_filter':"PTMlib_spermine_filter_corrected.csv", 'spd_filter':"PTMlib_spermidine_filter_corrected.csv"}
    
agent_ID_titration = {}
agent_titration = {}
for agent in agent_list:
    full_conc_list, full_mean_list, _ = read_titration_data(agent_filename[agent])
    #print (full_conc_list)
    num_list = sorted(agent_num_ID_fraction[agent].keys())
    for num in num_list:
        if num == 0:
            conc = 0.0
            mean = 1.0
        else:
            titr = name_titr[agent.split('_')[0] + str(num)] - 2
            if agent in ['PEG', 'HP1a']:
                assert full_conc_list[0] == 0
                titr +=1
            conc = full_conc_list[titr]
            mean = full_mean_list[titr]
        if agent not in agent_titration:
            agent_titration[agent] = {}
            agent_titration[agent]['conc'] = []
            agent_titration[agent]['survival'] = []
        agent_titration[agent]['conc'].append(conc)
        agent_titration[agent]['survival'].append(mean)
        for ID, fraction in agent_num_ID_fraction[agent][num].items():
            survival = mean * fraction
            #print (survival)
            if agent not in agent_ID_titration:
                agent_ID_titration[agent] = {}
            if ID not in agent_ID_titration[agent]:
                agent_ID_titration[agent][ID] = {}
            if 'conc' not in agent_ID_titration[agent][ID]:
                agent_ID_titration[agent][ID]['conc'] = []
            if 'survival' not in agent_ID_titration[agent][ID]:
                agent_ID_titration[agent][ID]['survival'] = []
            agent_ID_titration[agent][ID]['conc'].append(conc)
            agent_ID_titration[agent][ID]['survival'].append(survival)

for agent in agent_list:
    for ID in list(agent_ID_titration[agent].keys()):
        input = float(agent_ID_titration[agent][ID]['survival'][0])
        for i in range(len(agent_ID_titration[agent][ID]['survival'])):
            agent_ID_titration[agent][ID]['survival'][i] /= input

# plot survival probabilty over titrations
for agent in agent_list:
    fig = plt.figure()
    for ID in all_good_IDs:
        X = agent_ID_titration[agent][ID]['conc']
        Y = agent_ID_titration[agent][ID]['survival']
        p = plt.plot(X[1:], Y[1:], '.-', alpha=0.5, label=ID)
        plt.annotate(str(ID), (X[1], Y[1]), color=p[-1].get_color())
    plt.title(agent)
    plt.xlabel("Concentration")
    plt.ylabel("Soluble fraction")
    if agent in ['HP1a']:
        plt.xscale('log', basex=2)
    elif agent in ['sp', 'spd', 'CoH']:
        plt.xscale('log', basex=10)
    #plt.legend()
    #plt.show()
    plt.close()

#sys.exit(1)

#save survival probabilty data
#save actual number data
if False:
    total_WT_count = 1.6*(10**12)
    total_lib_count = total_WT_count / 100.0 #1% spike-in
    agent_list = ['sp']
    ID_list = sorted(list(set(all_good_IDs) - set([116])))
    for agent in agent_list:
        #f = open("PTMlib_%s_score.txt" % (agent), 'w')
        f = open("PTMlib_%s_count.txt" % (agent), 'w')

        # write down the header
        s = ""
        s += "ID" + '\t'
        s += '\t'.join(['count (%s %f mM)' % (agent, value) for value in agent_titration[agent]['conc']])
        #s += '\t'.join(['score (%s %f mM)' % (agent, x) for x in X[1:]])
        #s += '\t' + "mean score"
        print >> f, s

        # record WT background nucleosome data
        s = ""
        s += str(0) + '\t' 
        s += '\t'.join([str(int(round(total_WT_count*value))) for value in agent_titration[agent]['survival']])
        print >> f, s

        # record library data
        for ID in ID_list:
            s = ""
            s += str(ID) + '\t'
            fraction = agent_num_ID_fraction[agent][0][ID]
            s += '\t'.join([str(int(round(total_lib_count*fraction*value))) for value in agent_ID_titration[agent][ID]['survival']])
            #s += '\t'.join([str(-np.log2(y)) for y in Y[1:]])
            #s += '\t' + str(np.mean([-np.log2(y) for y in Y[1:]]))
            print >> f, s
            
        f.close()

    #sys.exit(1)


# define the condensabiltiy metrics
# fitting with sigmoidal curve and get "C-half"
# get condensabiltiy "Score" by calculating <-log(survival)>
def sigmoid(x, L ,x0, k):
    y = L / (1 + np.exp(k*(x-x0)))
    return (y)

agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a', 'sp_filter', 'spd_filter']
agent_ID_Chalf = {}
agent_ID_fitting = {}
agent_ID_score = {}
for agent in agent_list:
    fig = plt.figure(figsize=(2, 1.4))
    #for ID in [14, 47] + [42, 111, 112, 114, 115]:
    for ID in all_good_IDs:
        #if ID not in ['WT #1', 'WT #2', 'Cuttable DNA', 'Uncuttable DNA']:
        #    continue
        #if not ID.endswith('ub'):
        #    continue

        #fig = plt.figure()
        
        X = agent_ID_titration[agent][ID]['conc']
        Y = agent_ID_titration[agent][ID]['survival']

        X, Y = X[1:], Y[1:]

        #if agent in ['sp', 'spd']:
        #    X = X[:1] + X[2:]
        #    Y = Y[:1] + Y[2:]
        #elif agent in ['CoH']:
        #    X = X[:1] + X[4:]
        #    Y = Y[:1] + Y[4:]
            
        p0 = [max(Y), np.median(X), 1]
        bounds = ([0.0, 0.0, 0.0], [max(Y)+max(Y)*0.1, np.inf, np.inf])

        popt, pcov = curve_fit(sigmoid, X, Y, p0, bounds = bounds,  method='dogbox')
        residuals = np.asarray(Y)- sigmoid(X, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((np.asarray(Y)-np.mean(Y))**2)
        r_squared = 1 - (ss_res / ss_tot)
        #pred_X = np.linspace(min(X[1:]), max(X[1:]), 1000)
        pred_X = np.linspace(min(X), max(X), 1000)
        pred_Y = sigmoid(pred_X, *popt)

        #c = plt.plot(X[1:], Y[1:], '.', alpha=0.3)
        c = plt.plot(X, Y, '.', alpha=0.3, markersize=3.5)
        plt.plot(pred_X, pred_Y, '-', color=c[0].get_color(), alpha=0.3, label=ID)
        #plt.axvline(x=popt[1], linestyle='--', color=c[0].get_color(), alpha=0.5)

        #plt.title("%s %s" % (agent, ID))
        #plt.xlabel("Concentration")
        #plt.ylabel("Soluble fraction")
        #if agent in ['HP1a']:
        #    plt.xscale('log', basex=2)
        #elif agent in ['sp', 'spd', 'CoH']:
        #    plt.xscale('log', basex=10)
        #plt.show()
        #plt.close()

        if agent not in agent_ID_Chalf:
            agent_ID_Chalf[agent] = {}
        if ID not in agent_ID_Chalf[agent]:
            agent_ID_Chalf[agent][ID] = popt[1]

        if agent not in agent_ID_fitting:
            agent_ID_fitting[agent] = {}
        if ID not in agent_ID_fitting[agent]:
            agent_ID_fitting[agent][ID] = {}
        agent_ID_fitting[agent][ID]['para'] = popt
        agent_ID_fitting[agent][ID]['r-square'] = r_squared

        #score = np.mean(-np.log2(np.asarray(Y[1:])))
        score = np.mean(-np.log2(np.asarray(Y)))
        #score = -np.log2(np.mean(np.asarray(Y)))
        #score = np.mean(np.asarray(Y))

        if agent not in agent_ID_score:
            agent_ID_score[agent] = {}
        agent_ID_score[agent][ID] = score

        #print (popt)

    plt.gca().tick_params(axis='both', which='major', labelsize=5)
    plt.gca().tick_params(axis='both', which='minor', labelsize=5)
    #plt.title(agent_fullname[agent])
    #plt.xlabel("Concentration", fontsize=8)
    plt.xlabel("Spermine concentration (mM)", fontsize=8)
    plt.title("PTM library", fontsize=8)
    plt.ylabel("Soluble fraction", fontsize=8)
    if agent in ['HP1a']:
        plt.xscale('log', basex=2)
    elif agent in ['sp', 'spd', 'CoH']:
        plt.xscale('log', basex=10)
    #plt.legend()
    #plt.show()
    #plt.savefig(agent+'.svg', format='svg', bbox_inches='tight')
    plt.close()
#sys.exit(1)

# check r-square of fitting
for agent in agent_list:
    data = []
    for ID in agent_ID_fitting[agent]:
        r_squared = agent_ID_fitting[agent][ID]['r-square']
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

    
# all good IDs except WT controls
all_good_IDs_exceptWT = list(set(all_good_IDs) - set(cate_IDs['WT']))


# get difference w.r.t. wild type
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a', 'sp_filter', 'spd_filter']
ID_list = all_good_IDs_exceptWT
agent_ID_dChalf = {}
agent_ID_dscore = {}
for agent in agent_list:
    WT_Chalf, WT_score = [], []
    for ID in cate_IDs['WT']:
        WT_Chalf.append(agent_ID_Chalf[agent][ID])
        WT_score.append(agent_ID_score[agent][ID])
    WT_Chalf = np.mean(WT_Chalf)
    WT_score = np.mean(WT_score)

    for ID in ID_list:
        dChalf = agent_ID_Chalf[agent][ID] - WT_Chalf
        dscore = agent_ID_score[agent][ID] - WT_score

        if agent not in agent_ID_dChalf:
            agent_ID_dChalf[agent] = {}
        agent_ID_dChalf[agent][ID] = dChalf

        if agent not in agent_ID_dscore:
            agent_ID_dscore[agent] = {}
        agent_ID_dscore[agent][ID] = dscore

# pickle save
f = open("PTM_dscore.pickle", "wb")
pickle.dump(agent_ID_dscore, f)
f.close()
        
# find out outliers based on two metircs
agent_outliers = {}
for agent in agent_list:
    data_list = []
    for ID in ID_list:
        dChalf = agent_ID_dChalf[agent][ID]
        dscore = agent_ID_dscore[agent][ID]
        data_list.append([dChalf, dscore])

    clf = LocalOutlierFactor()
    outcheck = clf.fit_predict(data_list)

    for i in range(len(outcheck)):
        if outcheck[i] < 0:
            outlier = ID_list[i]
            if agent not in agent_outliers:
                agent_outliers[agent] = []
            agent_outliers[agent].append(outlier)

    fig = plt.figure(figsize=(2, 1.4))
    #fig = plt.figure()
    for ID in ID_list:
        dChalf = agent_ID_dChalf[agent][ID]
        dscore = agent_ID_dscore[agent][ID]
        if ID not in agent_outliers[agent]:
            plt.plot(dChalf, dscore, 'k.', markersize=2.5)
        else:
            #pass
            #plt.plot(dChalf, dscore, 'r.')
            plt.plot(dChalf, dscore, 'k.', markersize=2.5)
            #plt.annotate(ID_shortname[ID], (dChalf, dscore), fontsize=4, annotation_clip=True)

    plt.axvline(x=0, linestyle='--', color='k', alpha=0.5)
    plt.axhline(y=0, linestyle='--', color='k', alpha=0.5)
    plt.gca().tick_params(axis='both', which='major', labelsize=5)
    plt.gca().tick_params(axis='both', which='minor', labelsize=5)
    plt.title(agent_fullname[agent], fontsize=8)
    plt.xlabel("$\Delta$ C-half", fontsize=8)
    plt.ylabel("$\Delta$ Score", fontsize=8)
    #plt.savefig(agent + "_ChalfVSScore.svg", format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()
#sys.exit(1)

#sys.exit(1)

# get pseudo p-value based on WT distributions
agent_ID_pvalue = {}
for agent in agent_list:
    WTscores = [agent_ID_score[agent][ID] for ID in cate_IDs['WT']]
    WTmean = np.mean(WTscores)
    WTstd = np.std(WTscores)
    for ID in all_good_IDs_exceptWT:
        score = agent_ID_score[agent][ID]
        re_score = float(score-WTmean)/WTstd
        if re_score >=0:
            pvalue = 1 - norm.cdf(re_score)
        else:
            pvalue = norm.cdf(re_score)
        if agent not in agent_ID_pvalue:
            agent_ID_pvalue[agent] = {}
        agent_ID_pvalue[agent][ID] = pvalue

# dscore VS p-value scatter plot
# find outliers with p-value < 0.001
cutoff = 0.01
agent_outliers = {}
for agent in agent_list:
    if agent not in agent_outliers:
        agent_outliers[agent] = []
    
    fig = plt.figure()
    for ID in cate_IDs['WT+1PTM']:        
        dscore = agent_ID_dscore[agent][ID]
        pvalue = agent_ID_pvalue[agent][ID]
        if pvalue < cutoff:
            plt.plot(dscore, -np.log10(pvalue), 'r.')
            plt.annotate(ID_shortname[ID], (dscore, -np.log10(pvalue)))
            agent_outliers[agent].append(ID)
        else:
            plt.plot(dscore, -np.log10(pvalue), 'k.')
            
    plt.axvline(x=0, linestyle='--', color='k', alpha=0.5)
    plt.axhline(y=0, linestyle='--', color='k', alpha=0.5)
    plt.xlabel("$\Delta$ Score")
    plt.ylabel("-log(p-value)")
    plt.title(agent)
    #plt.show()
    plt.close()


# condense-seq QC
if False:
    # QC1: check the correlation with freeDNA contamination
    def read_freeDNA (fname):
        ID_freeDNA = {}
        First = True
        for line in open(fname):
            cols = line.split('\t')
            if First:
                First = False
                continue
            ID, name, _, _, _, freeDNA, _ = cols
            ID = int(ID)
            freeDNA = float(freeDNA[:-1])
            ID_freeDNA[ID] = freeDNA
        return ID_freeDNA
    ID_freeDNA = read_freeDNA('PTMlib_freeDNA.csv')

    agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a', 'sp_filter', 'spd_filter']
    ID_list = all_good_IDs

    #ID_list = []
    #for ID in all_good_IDs:
    #    if ID == 116:
    #        continue

    #    pminfo = ID_pminfo[ID]
    #
    #    mutations = []
    #    for subunit in subunit_list:
    #        for mtype in pminfo[subunit]['mutations']:
    #            for pos in pminfo[subunit]['mutations'][mtype]:
    #                mutations.append((subunit, mtype, pos))

    #   if len(mutations) > 1:
    #        continue
    #    if len(mutations) <= 0:
    #        ID_list.append(ID)
    #        continue

    #    subunit, mtype, pos = mutations[0]

    #    if mtype.startswith('me'):
    #        ID_list.append(ID)

    agent_ID_value_list = [agent_ID_Chalf, agent_ID_score]
    ylabel_list = ['C-half', 'Score']
    for agent in agent_list:
        for agent_ID_value, ylabel in list(zip(agent_ID_value_list, ylabel_list)):

            X, Y = [], []
            for ID in ID_list:
                freeDNA = ID_freeDNA[ID]
                X.append(freeDNA)
                value = agent_ID_value[agent][ID]
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
            for ID in ID_list:
            #for ID in cate_IDs['WT'] + cate_IDs['WT+mut']:
                x = ID_freeDNA[ID]
                y = agent_ID_value[agent][ID]
                if ID not in agent_outliers[agent]:
                    plt.plot(x, y, 'k.')
                else:
                    plt.plot(x, y, 'r.')
                    plt.annotate(ID, (x, y))
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
    #sys.exit(1)

    # QC2: check the correlation with input count
    agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a', 'sp_filter', 'spd_filter']
    ID_list = all_good_IDs
    for agent in agent_list:
        ID_count = agent_num_ID_count[agent][0]
        X, Y = [], []
        for ID in ID_list:
            count = ID_count[ID]
            score = agent_ID_score[agent][ID]
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
        for ID in ID_list:
            x = ID_count[ID]
            y = agent_ID_score[agent][ID]
            if ID not in agent_outliers[agent]:
                plt.plot(x, y, 'k.')
            else:
                plt.plot(x, y, 'r.')
                plt.annotate(ID, (x, y))
        plt.plot([X_Ypred[0][0], X_Ypred[-1][0]], [X_Ypred[0][1], X_Ypred[-1][1]], 'b--')
        plt.title(agent_fullname[agent])
        plt.xlabel("Input count")
        plt.ylabel("Score")
        #plt.savefig("Input_%s.png" % (agent))
        #plt.show()
        plt.close()

    # QC3: check the correlation with AT content
    agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a', 'sp_filter', 'spd_filter']
    ID_list = all_good_IDs
    for agent in agent_list:
        X, Y = [], []
        for ID in ID_list:
            AT = 100.0 - GC_content(ID_BC[ID])
            score = agent_ID_score[agent][ID]
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
        for ID in ID_list:
            x = 100.0 - GC_content(ID_BC[ID])
            y = agent_ID_score[agent][ID]
            if ID not in agent_outliers[agent]:
                plt.plot(x, y, 'k.')
            else:
                plt.plot(x, y, 'r.')
                plt.annotate(ID, (x, y))
        plt.plot([X_Ypred[0][0], X_Ypred[-1][0]], [X_Ypred[0][1], X_Ypred[-1][1]], 'b--')
        plt.title(agent_fullname[agent])
        plt.xlabel("AT content(%)")
        plt.ylabel("Score")
        #plt.savefig("ATcontent_%s.png" % (agent))
        #plt.show()
        plt.close()

    # QC4: check the correlation with filtered sample
    agent_list = ['sp', 'spd']
    ID_list = all_good_IDs
    for agent in agent_list:
        X, Y = [], []
        for ID in ID_list:
            X.append(agent_ID_score[agent][ID])
            Y.append(agent_ID_score[agent + '_filter'][ID])

        corr = scipy.stats.spearmanr(X, Y)[0]
        #corr = scipy.stats.pearsonr(X, Y)[0]
        #print ("%s VS %s: %1.2f" % (agent1, agent2, corr))
            
        fig = plt.figure()
        plt.plot(X, Y, '.')
        plt.annotate("Spearman %1.2f" % (corr), xy=(0.2, 0.75), fontsize=12, xycoords='axes fraction')
        plt.title(agent_fullname[agent])
        plt.xlabel("Score (50kD filter)")
        plt.ylabel("Score (100kD filter)")
        #plt.savefig("filter_%s.png" % (agent))
        #plt.show()
        plt.close()

    # QC5: compare with Oncohistone library condense-seq data
    def read_oncolib (fname):
        agent_mhistone_score = {}
        First = True
        for line in open(fname):
            cols = line.strip().split('\t')
            if First:
                agents = cols[1:]
                First = False
                continue
            mhistone, scores = cols[0], cols[1:]
            # modify mhistone format to fit PTMlibrary
            if mhistone.startswith('H'):
                if '/' in mhistone:
                    hname, mutations = mhistone.split(' ')
                    terms = mutations.split('/')
                    for i in range(len(terms)):
                        if not terms[i][0].isalpha():
                            terms[i] = mutations[0] + terms[i]
                        if not terms[i][-1].isalpha():
                            terms[i] += mutations[-1]
                    mhistone = hname + ','.join(terms)
                mhistone = ''.join(mhistone.split(' '))
            for agent, score in zip(agents, scores):
                if agent not in agent_mhistone_score:
                    agent_mhistone_score[agent] = {}
                assert mhistone not in agent_mhistone_score[agent]
                agent_mhistone_score[agent][mhistone] = float(score)
        return agent_mhistone_score
    agent_mhistone_score = read_oncolib("Oncolib_scores.txt")
    
    for agent in agent_mhistone_score:
        WT_scores = []
        for mhistone in agent_mhistone_score[agent]:
            if mhistone.startswith('WT'):
                WT_scores.append(agent_mhistone_score[agent][mhistone])
        agent_mhistone_score[agent]['WT'] = np.mean(WT_scores)

    oncolib_shortnames = agent_mhistone_score['sp'].keys()

    shortname_ID = {}
    for ID in all_good_IDs:
        shortname = ID_shortname[ID]
        if shortname == 'WT':
            continue
        assert shortname not in shortname_ID
        shortname_ID[shortname] = ID

    agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
    shortname_list = list(set(oncolib_shortnames) & set(shortname_ID.keys())) + ['WT']
    for agent in agent_list:
        X, Y = [], []
        for shortname in shortname_list:
            X.append(agent_mhistone_score[agent][shortname])
            if shortname == 'WT':
                WTscore = np.mean([agent_ID_score[agent][ID] for ID in cate_IDs['WT']])
                Y.append(WTscore)
            else:
                Y.append(agent_ID_score[agent][shortname_ID[shortname]])

        corr = scipy.stats.spearmanr(X, Y)[0]

        feature_list = [[x] for x in X]
        test_list = [[y] for y in Y]
        reg = linear_model.Ridge(alpha=0.5)
        reg.fit (feature_list, test_list)
        Ypred = reg.predict(feature_list)
        X_Ypred = sorted([[X[i], Ypred[i][0]] for i in range(len(Ypred))])


        fig = plt.figure()
        for i in range(len(X)):
            plt.plot(X[i], Y[i], 'go', ms=3)
            plt.annotate(shortname_list[i], (X[i], Y[i]))
        plt.plot([X_Ypred[0][0], X_Ypred[-1][0]], [X_Ypred[0][1], X_Ypred[-1][1]], 'b--')
        #plt.annotate("Spearman %1.2f" % (corr), xy=(0.2, 0.75), fontsize=12, xycoords='axes fraction')
        plt.title(agent_fullname[agent])
        plt.xlabel("Score (Oncolib)")
        plt.ylabel("Score (PTMlib)")
        #plt.savefig("OncoVSPTM_%s.png" % (agent))
        #plt.show()
        plt.close()

    #sys.exit(1)


# mean pvalue for (subunit, domain (fold/tail), PTM)
ID_list = cate_IDs['WT+1PTM']
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
subunit_foldrange = {'H2A':[17, 96], 'H2B':[38, 122], 'H3.1':[45, 131], 'H4':[31, 93]}
agent_section_pvalues = {}
ptm_list = set([])
for agent in agent_list:
    for ID in ID_list:
        minfo = ID_minfo[ID]
        pos, ptm = None, None
        for subunit in minfo:
            try:
                for pos in minfo[subunit]['mutations']:
                    ptm, aa, _ = minfo[subunit]['mutations'][pos]
                    break
            except:
                pass
            
        assert pos != None
        assert ptm != None

        if pos >= subunit_foldrange[subunit][0] and pos <= subunit_foldrange[subunit][1]:
            domain = 'fold'
        else:
            domain = 'tail'

        section = (subunit, domain, ptm)
        pvalue = agent_ID_pvalue[agent][ID]
        if agent not in agent_section_pvalues:
            agent_section_pvalues[agent] = {}
        if section not in agent_section_pvalues[agent]:
            agent_section_pvalues[agent][section] = []
        agent_section_pvalues[agent][section].append(pvalue)
        ptm_list.add(ptm)

hregion_list = []
subunits = ['H2A', 'H2B', 'H3.1', 'H4']
domains = ['fold', 'tail']
for i in range(len(subunits)):
    for j in range(len(domains)):
        subunit = subunits[i]
        domain = domains[j]
        hregion_list.append((subunit, domain))

ptm_list = list(ptm_list)

for agent in agent_list:
    img = np.zeros((len(hregion_list), len(ptm_list)))
    img[:] = np.nan
    for i in range(len(ptm_list)):
        for j in range(len(hregion_list)):
            subunit, domain = hregion_list[j]
            ptm = ptm_list[i]
            section = (subunit, domain, ptm)
            try:
                mpvalue = np.mean(agent_section_pvalues[agent][section])
                img[i,j] = mpvalue
            except:
                continue

    fig = plt.figure()
    plt.imshow(img)
    #plt.show()
    plt.close()

#sys.exit(1)
            
        
# plot ranking bar with histone modification information
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a', 'sp_filter', 'spd_filter']
ID_list = list(set(all_good_IDs_exceptWT) - set([116]))
subunit_list = ['H2A', 'H2B', 'H3', 'H4']
subunit_len = {'H2A':130, 'H2B':126, 'H3':136, 'H4':103}
subunit_color = {'H2A':'tab:purple', 'H2B':'tab:olive', 'H3':'tab:green', 'H4':'tab:pink'}
#hname_color = {'H2A':'tab:purple', 'H2A.X':'purple', 'H2A.Z':'darkviolet', 'H2B':'tab:olive', 'H3':'tab:green', 'H3.3':'green', 'H4':'tab:pink'}
mtype_color = {'ac':'red', 'me':'blue', 'ub':'green', 'ph':'yellow', 'cr':'m', 'GlcNAc':'tab:brown', 'mut':'gray'}

legend_elements = [Line2D([0], [0], marker='o', color='k', label='ac', mfc='red', mew=0.5),
                   Line2D([0], [0], marker='o', color='k', label='me', mfc='blue', mew=0.5),
                   Line2D([0], [0], marker='o', color='k', label='ub', mfc='green', mew=0.5),
                   Line2D([0], [0], marker='o', color='k', label='ph', mfc='yellow', mew=0.5),
                   Line2D([0], [0], marker='o', color='k', label='cr', mfc='m', mew=0.5),
                   Line2D([0], [0], marker='o', color='k', label='GlcNAc', mfc='tab:brown', mew=0.5),
                   Line2D([0], [0], marker='o', color='k', label='mut', mfc='gray', mew=0.5),
                   #Line2D([0], [0], marker='*', color='k', label='variant') # ppt version
                   Line2D([0], [0], marker='$\\ast$', color='k', label='variant')] # paper version 
matplotlib.rcParams['legend.handlelength'] = 0
matplotlib.rcParams['legend.numpoints'] = 1


his_len = 10
his_space = 2

for agent in agent_list:

    dChalf_ID, dscore_ID = [], []
    for ID in ID_list:
        dChalf = agent_ID_dChalf[agent][ID]
        dscore = agent_ID_dscore[agent][ID]
        dChalf_ID.append([dChalf, ID])
        dscore_ID.append([dscore, ID])
    dChalf_ID = sorted(dChalf_ID, reverse=True)
    dscore_ID = sorted(dscore_ID)

    #value_ID_list = [dChalf_ID, dscore_ID]
    #ylabel_list = ["$\Delta$ C-half", "$\Delta$ Score"]

    value_ID_list = [dscore_ID]
    ylabel_list = ["$\Delta$ Score"]

    
    for value_ID, ylabel in list(zip(value_ID_list, ylabel_list)):        

        #fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(18,8), sharex=True, gridspec_kw={'height_ratios': [1, 2]}) # ppt version
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(9,7), sharex=True, gridspec_kw={'height_ratios': [0.5, 2]}) # paper version

        X, Y = [], []
        xticks = []
        X1, Y1 = [], []
        X2, Y2 = [], []
        for i in range(len(value_ID)):
            value, ID = value_ID[i]
            minfo = ID_minfo[ID]

            X.append(i)
            Y.append(value)
            xticks.append(ID_shortname[ID])

            if ylabel == "$\Delta$ C-half":
                if value < 0 :
                    X1.append(i)
                    Y1.append(value)
                else:
                    X2.append(i)
                    Y2.append(value)
            else:
                if value >=0 :
                    X1.append(i)
                    Y1.append(value)
                else:
                    X2.append(i)
                    Y2.append(value)

            xpos = i
            for j in range(len(subunit_list)):
                subunit = subunit_list[j]

                yst = (his_len + his_space)*j
                axes[1].plot([xpos, xpos], [yst, yst+his_len], color=subunit_color[subunit], lw=1, solid_capstyle='round')

                hname = minfo[subunit]['name']
                if hname != subunit:
                    #axes[1].annotate('*', (xpos, yst), color='black', ha='center', va='center') # ppt version
                    axes[1].annotate('$\\ast$', (xpos, yst), color='black', ha='center', va='center', fontsize=5) # paper version

                for pos in sorted(minfo[subunit]['mutations']):
                    mtype, aa, mutation = minfo[subunit]['mutations'][pos]
                    ypos = yst + rescale ([pos], 0, subunit_len[subunit], 0, his_len)[0]
                    #axes[1].plot([xpos], [ypos], 'o', markersize=5, mfc=mtype_color[mtype], mew=0.5, mec='k') # ppt version
                    axes[1].plot([xpos], [ypos], 'o', markersize=3, mfc=mtype_color[mtype], mew=0.5, mec='k') # paper version

        #axes[0].bar(X, Y)
        axes[0].bar(X1, Y1, color='tab:blue')
        axes[0].bar(X2, Y2, color='tab:red')
        #axes[0].set_title(agent_fullname[agent]) # ppt version only
        axes[0].set_ylabel(ylabel) # ppt version
        axes[0].set_ylabel(ylabel, fontsize=7) # paper version

        space = abs(max(Y1) - min(Y2))*0.1
        axes[0].set_ylim([min(Y2)-2*space, max(Y1)+2*space])
        axes[0].set_yticks([min(Y2), 0, max(Y1)])
        axes[0].set_yticklabels([str(round(min(Y2),1)), str(0), str(round(max(Y1),1))], ha="center", va="center", fontsize=5, rotation=90) # paper only 

        axes[1].spines['top'].set_visible(False)
        axes[1].spines['left'].set_visible(False)
        axes[1].spines['right'].set_visible(False)
        axes[1].tick_params(top='off', left='off', right='off', labelleft='off', labelbottom='on')
        axes[1].set_xticks(X)
        #axes[1].set_xticklabels(xticks, fontsize=6, rotation=80, ha="right", rotation_mode="anchor") # ppt version
        axes[1].set_xticklabels(xticks, fontsize=5, rotation=90, ha="right", va='center', rotation_mode="anchor") # paper version


        for k in range(len(subunit_list)):
            subunit = subunit_list[k]
            #axes[1].annotate(subunit, (-3, (his_len + his_space)*k+his_len/2), fontsize=14, color='black', ha='center', va='center') # ppt version
            axes[1].annotate(subunit, (-3, (his_len + his_space)*k+his_len/2), fontsize=10, color='black', ha='center', va='center', rotation=90) # paper version

        #leg = axes[1].legend(handles=legend_elements, frameon=False, bbox_to_anchor=(0.97, 0.89), loc='upper left') # ppt version only

        plt.subplots_adjust(left=0.04, bottom=0.3, right=0.96, top=0.95, wspace=None, hspace=0.02)
        #plt.savefig("%s_%s_ladder_bar.png" % (agent, ylabel), dpi=300)
        plt.savefig("%s_%s_ladder_bar.svg" % (agent, ylabel), bbox_inches='tight')
        #plt.show()
        plt.close()
sys.exit(1)


# correlation with physical properties of modifications (TO DO)

# PCA analysis of histones by combining all condensing agents scores
# nucleosome with single PTMs
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
#ID_list = all_good_IDs
ID_list = cate_IDs['WT'] + cate_IDs['WT+1PTM']
#cate_IDs['WT+1mut'] + cate_IDs['WT+CpGme'] + cate_IDs['H2A.Z'] + cate_IDs['H2A.X'] + cate_IDs['H3.3'] 
#ID_list = cate_IDs['WT+1PTM']

mtype_IDs = {}
for ID in ID_list:
    shortname = ID_shortname[ID]
    if shortname != 'WT':
        hname, pos_mutation = mhistone_parser(shortname)
        mtype, _, _ = pos_mutation.values()[0]
    else:
        mtype = 'WT'
    if mtype not in mtype_IDs:
        mtype_IDs[mtype] = []
    mtype_IDs[mtype].append(ID)

mtype_color = {'WT':'black',
               'ac':'red',
               'me':'blue',
               'ub':'green',
               'ph':'yellow',
               'cr':'m',
               'GlcNAc':'tab:brown',
               'mut':'gray'}

ID_state = {}
for agent in agent_list:
    for ID in ID_list:
        score = agent_ID_score[agent][ID]
        if ID not in ID_state:
            ID_state[ID] = []
        ID_state[ID].append(score)


X = [ID_state[ID] for ID in ID_list]
        
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
ID_PCA = {}
x_list, y_list = [], []
for i in range(len(Xr)):
    ID = ID_list[i]
    x, y = Xr[i][0], Xr[i][1]
    ID_PCA[ID] = (x, y)
    x_list.append(x)
    y_list.append(y)
    if outcheck[i] < 0:
        outliers.append(ID)

#fig = plt.figure(figsize=(8,7))
fig = plt.figure(figsize=(6,5))
mtype_list = ['WT', 'ac','me', 'ub', 'ph', 'cr', 'GlcNAc']
for mtype in mtype_list:
    IDs = mtype_IDs[mtype]
    for i in range(len(IDs)):
        ID = IDs[i]
        x, y = ID_PCA[ID]
        if i == 0:
            label = mtype
        else:
            label = None
        color = mtype_color[mtype]
        plt.plot(x, y, 'o', color=color, mew=0.5, mec='k', label=label)
        if ID in outliers:
            plt.annotate(ID_shortname[ID], (x, y), fontsize=6)

plt.title("PCA plot", fontsize=10)
plt.xlabel('PC1', fontsize=8)
plt.ylabel('PC2', fontsize=8)
plt.gca().tick_params(axis='both', which='major', labelsize=5)
plt.gca().tick_params(axis='both', which='minor', labelsize=5)
plt.legend(fontsize=6)
#plt.savefig("PCA_allagent.svg", format='svg', bbox_inches='tight')
#plt.tight_layout()
#plt.show()
plt.close()

#sys.exit(1)

# all PTM library
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
ID_list = list(set(all_good_IDs) - set([116]))

ID_state = {}
for agent in agent_list:
    for ID in ID_list:
        score = agent_ID_score[agent][ID]
        if ID not in ID_state:
            ID_state[ID] = []
        ID_state[ID].append(score)


X = [ID_state[ID] for ID in ID_list]
        
#scaler = preprocessing.StandardScaler(with_std=False).fit(X)
scaler = preprocessing.StandardScaler().fit(X)
X_scaled = scaler.transform(X)
pca = PCA().fit(X_scaled)
Xr = pca.transform(X_scaled)

#pca = PCA().fit(X)
#Xr = pca.transform(X)

clf = LocalOutlierFactor()
#outcheck = clf.fit_predict(X)
outcheck = clf.fit_predict( [row[:2] for row in Xr] )

outliers = []
ID_PCA = {}
x_list, y_list = [], []
for i in range(len(Xr)):
    ID = ID_list[i]
    x, y = Xr[i][0], Xr[i][1]
    ID_PCA[ID] = (x, y)
    x_list.append(x)
    y_list.append(y)
    if outcheck[i] < 0:
        outliers.append(ID)

group_color = {"H2A/Bac":"magenta",
               "H3ac":"tab:pink",
               "H4ac":"tab:red",
               "KpolyAc":"red",
               "H4ac+H3me":"saddlebrown",
               "H3me":"blue",
               "H4me":"tab:blue",
               "H3cr":"blueviolet",
               "GlcNAc":"purple",
               'WT':"black",
               "WT+CpGme":"dimgray",
               'WT+mut':"tab:gray",
               "+ub":"green",
               'Var':"darkslategray",
               'Var+mut':"teal",
               "H3ph":"tab:olive",
               "AP mutant":"cyan"}

group_marker = {"H2A/Bac":"o",
                "H3ac":"o",
                "H4ac":"h",
                "KpolyAc":"X",
                "H4ac+H3me":"<",
                "H3me":"s",
                "H4me":"s",
                "H3cr":"^",
                "GlcNAc":"h",
                'WT':"8",
                "WT+CpGme":"X",
                'WT+mut':"8",
                "+ub":"D",
                'Var':"p",
                'Var+mut':"p",
                "H3ph":"P",
                "AP mutant":"X"}

group_list = ['WT', "WT+CpGme", 'WT+mut', 'Var', 'Var+mut', "AP mutant", "H2A/Bac", "H3ac", "H4ac", "KpolyAc", "H3me", "H4me", "H4ac+H3me", "H3cr","GlcNAc","H3ph", "+ub"]

#fig = plt.figure(figsize=(12,12))
fig = plt.figure(figsize=(6,5))
for group in group_list:
    IDs = group_IDs[group]
    color = group_color[group]
    marker = group_marker[group]
    alpha = 1.0
    mew = 0.5 
    for i in range(len(IDs)):
        ID = IDs[i]
        x, y = ID_PCA[ID]
        if i == 0:
            plt.plot(x, y, '.', color=color, marker=marker, mew=mew, mec='k', alpha=alpha, label=group)
        plt.plot(x, y, '.', color=color, marker=marker, mew=mew, mec='k', alpha=alpha)
        #if ID in outliers:
        #    plt.annotate(ID_shortname[ID], (x,y), fontsize=6)

#plt.xlim([-6, 15])
#plt.ylim([-2, 4])
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

# scatter plot of scores for each group
group_list = ['WT', "WT+CpGme", 'WT+mut', 'Var', 'Var+mut', "AP mutant", "H2A/Bac", "H3ac", "H4ac", "KpolyAc", "H3me", "H4me", "H4ac+H3me", "H3cr","GlcNAc","H3ph", "+ub"]
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
agent_color = {'sp':'tab:red',
               'spd':'tab:green',
               'CoH':'tab:blue',
               'PEG':'tab:olive',
               'HP1a':'tab:purple'}

fig, axes = plt.subplots(nrows=len(agent_list), ncols=1, figsize=(15,10), sharex=True)
for agent, ax in zip(agent_list, axes):
    WTscore = np.mean([agent_ID_score[agent][ID] for ID in group_IDs['WT']])
    X, Y = [], []
    Ymean = []
    for i in range(len(group_list)):
        group = group_list[i]
        IDs = group_IDs[group]
        values = []
        for ID in IDs:
            X.append(i)
            Y.append(agent_ID_score[agent][ID]-WTscore)
            values.append(agent_ID_score[agent][ID]-WTscore)
        Ymean.append(np.mean(values))
    ax.plot(X, Y, 'o', color=agent_color[agent], markersize=3, alpha=0.5, label=agent_fullname[agent])
    ax.plot(range(len(Ymean)), Ymean, 'x-', color=agent_color[agent], markersize=5, mec='k', mew=3)
    #ax.set_title(agent_fullname[agent])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    leg = ax.legend(prop={'size': 12}, loc='upper right', frameon=False, handlelength=0, handletextpad=0)
    for item in leg.legendHandles:
        item.set_visible(False)
ax.set_xticks(range(len(group_list)))
ax.set_xticklabels(group_list, fontsize=8, rotation=45, ha="right", rotation_mode="anchor")
fig.text(0.08, 0.5, 'Score (A.U.)', va='center', rotation='vertical')
#plt.savefig("score_group.png", bbox_inches='tight')
#plt.show()
plt.close()

#sys.exit(1)
    

# compare between condensing agents
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
ID_list = all_good_IDs
pair_corr = {}
corr_matrix = [ [np.nan for i in range(len(agent_list))] for i in range(len(agent_list)) ]
for i in range(len(agent_list)-1):
    for j in range(i+1, len(agent_list)):
        agent1 = agent_list[i]
        agent2 = agent_list[j]

        X, Y = [], []
        for ID in ID_list:
            X.append(agent_ID_score[agent1][ID])
            Y.append(agent_ID_score[agent2][ID])

        corr = scipy.stats.spearmanr(X, Y)[0]
        #corr = scipy.stats.pearsonr(X, Y)[0]
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
            X = [ agent_ID_score[agent1][ID] for ID in ID_list ]
            Y = [ agent_ID_score[agent2][ID] for ID in ID_list ]
            axes[i,j].plot(X, Y, 'k.', markersize=0.5)
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
            axes[i,j].text(len(agent_list)/2, len(agent_list)/2, s, ha="center", va="center", fontsize=6, weight='bold')
            axes[i,j].set_xlim([0, len(agent_list)-1])
            axes[i,j].set_ylim([0, len(agent_list)-1])
            axes[i,j].set_axis_off()
        else:
            assert i < j
            value = pair_corr[(agent1, agent2)]
            matrix = np.zeros((len(agent_list), len(agent_list)))
            matrix[:] = value
            img = axes[i,j].imshow(matrix, cmap="jet", vmin=0.0, vmax=1.0, origin='lower')
            if value < 0.3 or value > 0.7:
                color = "white"
            else:
                color = "black"
            axes[i,j].text(len(agent_list)/2, len(agent_list)/2, str(round(value,2)), ha="center", va="center", fontsize=6, color=color, weight='bold')
            axes[i,j].set_xlim([0, len(agent_list)-1])
            axes[i,j].set_ylim([0, len(agent_list)-1])
            axes[i,j].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labeltop=False, labelleft=False, labelright=False)

plt.subplots_adjust(wspace=0.2, hspace=0.2)
cbar=fig.colorbar(img, ax=axes, location='right', shrink=0.6, aspect=30, ticks=[0, 1.0])
cbar.ax.set_yticklabels([str(0), str(1)], fontsize=5)
cbar.ax.set_ylabel('Spearman correlation', rotation=-90, va="bottom", fontsize=5, labelpad=-5)
#plt.suptitle("Corrrelation betwen condensing agents")
#plt.show()
#plt.savefig("PTM_corr_btw_agent.svg", format='svg', bbox_inches='tight')
plt.close()

#sys.exit(1)

# compare with chromatin remodeling data
def read_remodellers (fname):
    enzyme_ID_rate = {}
    #enzyme_ID_foldchange = {}
    line_count = 0
    for line in open(fname):
        line = line.strip()
        if not line:
            break
        cols = line.split('\t')
        if line_count == 0:
            enzyme_list = []
            for col in cols:
                if col.strip():
                    enzyme_list.append(col)
                    enzyme_ID_rate[col] = {}
                    #enzyme_ID_foldchange[col] = {}
            line_count +=1
            #print enzyme_list
            continue
        if line_count < 2:
            line_count +=1
            continue
        for i in range(len(enzyme_list)):
            enzyme = enzyme_list[i]
            ID = cols[11*i]
            try:
                ID = int(ID)
            except:
                if ID == 'DNA standard 2':
                    ID = 117
            rate = float(cols[11*i + 3])
            foldchange = float(cols[11*i + 9])

            assert ID not in enzyme_ID_rate[enzyme]
            #enzyme_ID_rate[enzyme][ID] = rate
            enzyme_ID_rate[enzyme][ID] = foldchange
            #assert ID not in enzyme_ID_foldchange[enzyme]
            #enzyme_ID_foldchange[enzyme][ID] = foldchange
        line_count +=1
    return enzyme_list, enzyme_ID_rate
    #return enzyme_list, enzyme_ID_rate, enzyme_ID_foldchange
enzyme_list, enzyme_ID_rate  = read_remodellers("PTMlib_remodellers.csv")

agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
ID_list = list(set(all_good_IDs) - set([116]))

fig, axes = plt.subplots(nrows=len(agent_list), ncols=len(enzyme_list))
for i in range(len(agent_list)):
    for j in range(len(enzyme_list)):
        agent = agent_list[i]
        enzyme = enzyme_list[j]
                
        X = [enzyme_ID_rate[enzyme][ID] for ID in ID_list]
        Y = [agent_ID_score[agent][ID] for ID in ID_list]

        corr = scipy.stats.spearmanr(X, Y)[0]
        #print ("%s VS %s: %1.2f" % (enzyme, agent, corr))

        axes[i, j].plot(X, Y, 'k.', markersize=1.0)
        #axes[i, j].set_yscale('log', basey=2)

        if j <= 0:
            axes[i, j].set_ylabel(agent, rotation=0, labelpad=20)
        else:
            axes[i, j].tick_params(axis='y', which='both', labelleft=False)

        if i >= len(agent_list) - 1:
            axes[i, j].set_xlabel('%s' % (enzyme))
        else:
            axes[i, j].tick_params(axis='x', which='both', labelbottom=False)

#plt.tight_layout()
fig.suptitle("Corrrelation with remodellers rate")
#plt.show()
plt.close()

fig, axes = plt.subplots(nrows=len(agent_list), ncols=len(enzyme_list))
for i in range(len(agent_list)):
    for j in range(len(enzyme_list)):

        agent = agent_list[i]
        enzyme = enzyme_list[j]
                
        X = [enzyme_ID_rate[enzyme][ID] for ID in ID_list]
        Y = [agent_ID_score[agent][ID] for ID in ID_list]

        corr = scipy.stats.spearmanr(X, Y)[0]

        matrix = np.zeros((len(agent_list), len(agent_list)))
        matrix[:] = corr
        img = axes[i, j].imshow(matrix, cmap="jet", vmin=0.0, vmax=1.0, origin='lower')

        if corr < 0.3 or corr > 0.7:
            color = 'white'
        else:
            color = 'black'

        axes[i,j].text(len(agent_list)/2, len(agent_list)/2, str(round(corr,2)), ha="center", va="center", fontsize=10, color=color, weight='bold')
        axes[i,j].set_xlim([0, len(agent_list)-1])
        axes[i,j].set_ylim([0, len(agent_list)-1])
        axes[i,j].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labeltop=False, labelleft=False, labelright=False)

        if j <= 0:
            axes[i, j].set_ylabel(agent, rotation=0, labelpad=20)
        else:
            axes[i, j].tick_params(axis='y', which='both', labelleft=False)

        if i >= len(agent_list) - 1:
            axes[i, j].set_xlabel('%s' % (enzyme))
        else:
            axes[i, j].tick_params(axis='x', which='both', labelbottom=False)

#plt.tight_layout()
fig.suptitle("Corrrelation with remodellers rate")
cbar=fig.colorbar(img, ax=axes, location='right', shrink=0.8)
cbar.ax.set_ylabel('Spearman correlation', rotation=-90, va="bottom")
#plt.show()
plt.close()

#sys.exit(1)


# synergy analysis
#ID_list = list(set(all_good_IDs)&set(cate_IDs['WT+PTM']))
ID_list = list(set(all_good_IDs)&set(cate_IDs['WT+1PTM']))

ID_mutations = {}
mutations_ID = {}
mutnum_IDs = {}
for ID in ID_list:
    pminfo = ID_pminfo[ID]
    mutations = set([])
    for subunit in subunit_list:
        hname = pminfo[subunit]['name']
        for mtype in pminfo[subunit]['mutations']:
            for pos in pminfo[subunit]['mutations'][mtype]:
                aa, mutation = pminfo[subunit]['mutations'][mtype][pos]
                mutations.add((hname, pos, mtype, aa, mutation))
    mutations = frozenset(mutations)
    
    ID_mutations[ID] = mutations
    mutations_ID[mutations] = ID

    mutnum = len(mutations)
    if mutnum not in mutnum_IDs:
        mutnum_IDs[mutnum] = []
    mutnum_IDs[mutnum].append(ID)


IDs = set([])
for i in range(len(mutnum_IDs[1])-1):
    for j in range(i+1, len(mutnum_IDs[1])):
        ID1, ID2 = mutnum_IDs[1][i], mutnum_IDs[1][j]
        tmutations = frozenset(set(ID_mutations[ID1]) | set(ID_mutations[ID2]))
        try:
            mutations_ID[tmutations]
            IDs.add(ID1)
            IDs.add(ID2)
        except:
            continue
IDs = list(IDs)

for agent in agent_list:
    ID_score = agent_ID_score[agent]
    matrix = np.zeros((len(IDs), len(IDs)))
    for i in range(len(IDs)):
        for j in range(i, len(IDs)):
            if i == j:
                synergy = np.nan
            else:
                ID1, ID2 = IDs[i], IDs[j]
                tmutations = frozenset(set(ID_mutations[ID1]) | set(ID_mutations[ID2]))
                try:
                    tID = mutations_ID[tmutations]
                    synergy = ID_score[tID] - 0.5*(ID_score[ID1] + ID_score[ID2])
                except:
                    synergy = np.nan
            matrix[i][j] = synergy
            matrix[j][i] = synergy

    fig = plt.figure()
    cmap = cm.bwr_r
    cmap.set_bad(color='black')
    plt.imshow(matrix, cmap=cmap, vmin=-2.0, vmax=2.0)
    plt.xticks(range(len(IDs)), [str(get_mname(ID_minfo[ID])) for ID in IDs], rotation=45)
    plt.yticks(range(len(IDs)), [str(get_mname(ID_minfo[ID])) for ID in IDs])
    plt.title(agent_fullname[agent])
    plt.colorbar()
    #plt.show()
    plt.close()

#sys.exit(1)
    

### structure analysis
# property analsyis
ID_list = list(set(all_good_IDs) - set([116]))
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
subunit_list = ['H2A', 'H2B', 'H3', 'H4']
subunit_chains = {'H2A':['C'], 'H2B':['D'], 'H3.1':['A'], 'H4':['B']}

for ID in ID_list:
    NCP = Molecules("1kx5")
    mutation_chain_resi = {}
    minfo = ID_minfo[ID]
    for subunit in subunit_list:
        hname = ID_minfo[ID][subunit]['name']
        pos_mutation = ID_minfo[ID][subunit]['mutations']
        for pos in pos_mutation:
            mtype, aa, mutation = pos_mutation[pos]
            
            if mutation not in mutation_chain_resi:
                mutation_chain_resi[mutation] = {}
            for chain in subunit_chains[subunit]:
                if chain not in mutation_chain_resi[mutation]:
                    mutation_chain_resi[mutation][chain] = []
                resi = pos_resi[pos]
                mutation_chain_resi[mutation][chain].append(resi)

    for mutation in mutation_chain_resi:
        chain_resi = mutation_chain_resi[mutation]
        mut, num = re.findall('([A-Z]+)(\d+(?:,\d+)*)', mutations)
                
# plot histone cartoon
#ID_list = list(set(all_good_IDs) & set(cate_IDs['WT+PTM']))
ID_list = list(set(all_good_IDs) & set(cate_IDs['WT+1PTM']))
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
subunit_list = ['H2A', 'H2B', 'H3', 'H4']
location_list = ['fold', 'tail']

#subunit_chains = {'H2A':['C', 'G'], 'H2B':['D', 'H'], 'H3.1':['A', 'E'], 'H4':['B','F']}
#subunit_chains = {'H2A':['C'], 'H2B':['D'], 'H3.1':['E'], 'H4':['F']}
subunit_foldrange = {'H2A':[17, 96], 'H2B':[38, 122], 'H3':[45, 131], 'H4':[31, 93]}

# rescale the data in old range (old_st, old_ed) into new range (new_st, new_ed)
def rescale (value_list, old_st, old_ed, new_st, new_ed):
    output = []
    for value in value_list:
        assert value >= old_st and value <= old_ed
        new_value = new_st + (new_ed - new_st)*float(value-old_st)/(old_ed-old_st)
        output.append(new_value)
    return outpu

# find IDs with single PTM
subunit_infos = {}
for ID in ID_list:
    pminfo = ID_pminfo[ID]
    mutations = []
    for subunit in subunit_list:
        for mtype in pminfo[subunit]['mutations']:
            for pos in pminfo[subunit]['mutations'][mtype]:
                mutations.append((subunit, mtype, pos))

    if len(mutations) > 1:
        continue

    subunit, mtype, pos = mutations[0]

    if subunit not in subunit_infos:
        subunit_infos[subunit] = []
    subunit_infos[subunit].append((ID, pos, mtype))
    

agent_pdscores = {}
for agent in agent_list:
    for subunit in subunit_infos:
        for ID, pos, mtype in subunit_infos[subunit]:
            dscore = agent_ID_dscore[agent][ID]
            if agent not in agent_pdscores:
                agent_pdscores[agent] = []
            agent_pdscores[agent].append(abs(dscore))


mtype_marker = {'ac':'o', 'me':'s', 'ub':'D', 'ph':'P', 'cr':'p', 'GlcNAc':'*'}

legend_elements = [Line2D([0], [0], marker='o', color='k', label='acetylation', ms=10, markerfacecolor='w'),
                   Line2D([0], [0], marker='s', color='k', label='methylation', ms=10, markerfacecolor='w'),
                   Line2D([0], [0], marker='D', color='k', label='ubiquitylation', ms=10, markerfacecolor='w'),
                   Line2D([0], [0], marker='P', color='k', label='phosphorylation', ms=10, markerfacecolor='w'),
                   Line2D([0], [0], marker='p', color='k', label='crotonylation', ms=10, markerfacecolor='w'),
                   Line2D([0], [0], marker='*', color='k', label='N-acetylglucosamine', ms=10, markerfacecolor='w')]
matplotlib.rcParams['legend.handlelength'] = 0
matplotlib.rcParams['legend.numpoints'] = 1
    
fig = plt.figure()
ax = plt.gca()
leg = ax.legend(handles=legend_elements, frameon=False, fontsize=20)
#for lh in leg.legendHandles: 
#    lh.set_markersize(1000.0)
#plt.savefig('legend.svg', format='svg', bbox_inches='tight')
plt.show()
plt.close()

sys.exit(1)

for agent in []:
    for subunit in subunit_infos:
        for ID, pos, mtype in subunit_infos[subunit]:
            aa, mutation = ID_pminfo[ID][subunit]['mutations'][mtype][pos]
            dscore = agent_ID_dscore[agent][ID]
            min_absdscore, max_absdscore = min(agent_pdscores[agent]), max(agent_pdscores[agent])
            s = rescale([abs(dscore)], min_absdscore, max_absdscore, 5, 100)[0]

            if dscore > 0:
                color = 'blue'
                fontcolor = 'white'
            else:
                color = 'red'
                fontcolor = 'yellow'
            
            fig = plt.figure()
            plt.plot([0,0],[0,0], marker=mtype_marker[mtype], markersize=s, mfc=color, markeredgewidth=1.5, mec='k')
            if len(mutation) > 3:
                fontsize = int(0.3*s)
            else:
                fontsize = int(0.4*s)
            if fontsize >= 9:
                plt.annotate(mutation.title(), (0, 0), fontsize=fontsize, color=fontcolor, weight='bold', ha='center', va='center')
            #plt.annotate(aa+str(pos)+mutation, (0, 0), fontsize=int(0.3*s), ha='center', va='center')
            #plt.annotate(aa+str(pos)+mutation, (0, 0), ha='center', va='center')
            #plt.xlim([-0.2, 0.2])
            #plt.ylim([-0.2, 0.2])
            plt.gca().axis('off')
            #plt.savefig(agent + '_' + get_mname(ID_minfo[ID])+".png", bbox_inches='tight', pad_inches=0, transparent=True)
            plt.close()
            
sys.exit(1)
    

for agent in agent_list:
    fig = plt.figure()
    for i in range(len(subunit_list)):
        subunit = subunit_list[i]
        mtype_IDs = subunit_location_mtype_IDs[subunit]['tail']
        for mtype in mtype_IDs:
            if mtype == 'ac':
                color = 'red'
            elif mtype == 'me':
                color = 'blue'
            else:
                color = 'black'
            for ID in mtype_IDs[mtype]:
                pos = ID_pminfo[ID][subunit]['mutations'][mtype].keys()[0]
                if pos < subunit_foldrange[subunit][0]:
                    terminal = 'N'
                    tail_len = subunit_foldrange[subunit][0]
                else:
                    assert pos > subunit_foldrange[subunit][1]
                    terminal = 'C'
                    tail_len = None
                r = get_r (pos, tail_len, 0, terminal)
                s = agent_ID_score[agent][ID]
                plt.polar(np.pi/4.0 + i*np.pi/2.0, r, 'o', color=color, markersize=3, alpha=0.5)

    plt.show()
    plt.close()
                


sys.exit(1)
    
for agent in agent_list:
    fig = plt.figure()
    for i in range(len(subunit_list)):
        for j in range(len(location_list)):
            subunit = subunit_list[i]
            location = location_list[j]
            mtype_IDs = subunit_location_mtype_IDs[subunit][location]
            for mtype, IDs in mtype_IDs.items():
                if mtype == 'me':
                    color = 'blue'
                elif mtype == 'ac':
                    color = 'red'
                else:
                    color = 'black'
                for ID in IDs:
                    score = agent_ID_score[agent][ID]
                    x = i + (j-0.5)
                    y = score
                    plt.plot([x], [y], '.', makercolor=color)
    plt.show()
    plt.close()
    
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
