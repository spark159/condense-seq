import math
import re
import numpy as np
import os, sys, subprocess, re
import Interval_dict
import load_file
import matplotlib.pyplot as plt

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
all_good_IDs = sorted(list(set(ID_BC.keys()) - set([71, 78, 113, 116, 117])))

# categorize IDs
#['WT+2PTM', 'WT+10PTM', 'WT+4PTM', 'H2A.X', 'WT+1PTM', 'WT+1mut', 'WT+2mut', 'WT+5PTM', 'H2A.Z', 'WT+3mut', 'H3.3+1mut', 'H2A.Z/H3.3', 'WT', 'WT+6PTM', 'WT+8mut', 'WT+19PTM', 'H3.3', 'WT+14PTM', 'WT+CpGme']
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

# read condensation score file
def read_PTM_score(fname):
    ID_score = {}
    First = True
    for line in open(fname):
        cols = line.strip().split('\t')
        if First:
            assert cols[-1] == 'mean score'
            First = False
            continue
        ID = int(cols[0])
        score = float(cols[-1])
        assert ID not in ID_score
        ID_score[ID] = score
    return ID_score
ID_score = read_PTM_score("PTMlib_sp_score.txt")

# get dscores
WT_score = np.mean([ID_score[ID] for ID in cate_IDs['WT']])
all_good_IDs_exceptWT = list(set(all_good_IDs) - set(cate_IDs['WT']))
ID_dscore = {}
for ID in all_good_IDs_exceptWT:
    dscore = ID_score[ID] - WT_score
    assert ID not in ID_dscore
    ID_dscore[ID] = dscore
    
# make new IDs
target_cates = ['WT+1PTM', 'WT+2PTM', 'WT+4PTM', 'WT+5PTM', 'WT+6PTM', 'WT+10PTM', 'WT+14PTM', 'WT+19PTM', 'H2A.X', 'H2A.Z', 'H3.3', 'H2A.Z/H3.3']
subunit_list = ['H2A', 'H2B', 'H3', 'H4']

newID_dscore = {}
for cate in target_cates:
    for ID in cate_IDs[cate]:
        mutations = []
        minfo = ID_minfo[ID]
        for subunit in subunit_list:
            hname = minfo[subunit]['name']
            mutation = minfo[subunit]['mutations']
            if hname != subunit:
                mutations.append(hname)
            for pos in mutation:
                mtype, aa, ptm = mutation[pos]
                mname = ''.join([subunit, aa.upper(), str(pos), ptm])
                mutations.append(mname)
        assert len(mutations) > 0
        if len(mutations) < 2:
            newID = mutations[0]
        else:
            newID = tuple(sorted(mutations))
        assert newID not in newID_dscore
        newID_dscore[newID] = ID_dscore[ID]
newID_dscore['WT'] = 0.0

# read binary PTM file and write predicted score
def combine_scores (scores, metric):
    if metric == 'average':
        return np.mean(scores)
    elif metric == 'override':
        return sorted([(abs(score), score) for score in scores])[-1][1]
    elif metric == 'sum':
        return sum(scores)
        
chr = 'chr1'
binsize = 200
metric = 'average'

infname = 'H1_%s_binary.txt' % (chr)
outfname = 'H1_sp_%s_PredScore_%s' % (chr, metric)

f = open(outfname + '_anot.cn', 'w')
print >> f, 'SNP\tChromosome\tPhysicalPosition\tPredScore'

linecount = 0
k = 0
for line in open(infname):
    cols = line.strip().split('\t')
    linecount +=1
    if linecount < 2:
        continue
    if linecount == 2:
        marklist = []
        skiplist = []
        for mark in cols:
            if mark == 'H2AFZ':
                mark = 'H2A.Z'
            try:
                newID_dscore[mark]
            except:
                print 'no data for %s' % (mark)
                skiplist.append(mark)
            marklist.append(mark)
        continue
    
    mutations = []
    for i in range(len(cols)):
        value = int(cols[i])
        if value == 1:
            mark = marklist[i]
            if mark not in skiplist:
                mutations.append(mark)
            
    if len(mutations) < 1:
        PredScore = 0.0
    elif len(mutations) == 1:
        newID = mutations[0]
        PredScore = newID_dscore[newID]
    else:
        try:
            newID = tuple(sorted(mutations))
            PredScore = newID_dscore[newID]
        except:
            scores = []
            for mark in mutations:
                    scores.append(newID_dscore[mark])
            PredSCore = combine_scores(scores, metric=metric)

    pos = k*binsize + binsize/2
    print >> f, '%d\t%s\t%d\t%.5f' % (k, chr, pos, PredScore)
    k +=1

f.close()                


# check the predicted scores for each chroHMM states
chr = 'chr1'
metric = 'average'
fname = 'H1_sp_%s_PredScore_%s_anot.cn' % (chr, metric)
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(fname)

def read_chromHMM(fname, chr_target, change=False):
    state_intervals = {}
    for line in open(fname):
        cols = line.strip().split()
        chr, st, ed, state = cols[:4]
        if chr != chr_target:
            continue
        st, ed = int(st), int(ed)
        if change:
            state = change[state]
        if state not in state_intervals:
            state_intervals[state] = []
        state_intervals[state].append((st,ed))
    return state_intervals

name_dict = {"E1":"Polycomb repressed",
             "E2":"Poised promoter",
             "E3":"Weak promoter",
             "E4":"Strong enhancer",
             "E5":"Active promoter",
             "E6":"Weak enhancer",
             "E7":"Quiescence1",
             "E8":"Quiescence2",
             "E9":"Heterochromatin",
             "E10":"Tx elongation",
             "E11":"Weak Tx",
             "E12":"Insulator"}

state_intervals = read_chromHMM("H1_12_segments.bed", chr_target='chr1', change=name_dict)

dID_interval = {}
for state in state_intervals:
    intervals = state_intervals[state]
    for i in range(len(intervals)):
        dID = state + ':' + str(i)
        assert dID not in dID_interval
        dID_interval[dID] = intervals[i]

dinterval_dict = Interval_dict.double_hash(dID_interval, 10000, 250000000)

state_name_values = {}
name_state_values = {}
for ID in ID_pos:
    pos = ID_pos[ID]
    dIDs = dinterval_dict.find(pos)
    if not dIDs:
        continue
    for dID in dIDs:
        state, _ = dID.split(':')
        if state not in state_name_values:
            state_name_values[state] = {}
        for name in name_ID_value:
            if name == 'Sequence':
                continue
            if name not in state_name_values:
                state_name_values[state][name] = []
            if name not in name_state_values:
                name_state_values[name] = {}
            if state not in name_state_values[name]:
                name_state_values[name][state] = []
            value = name_ID_value[name][ID]
            if np.isnan(value):
                continue
            state_name_values[state][name].append(value)
            name_state_values[name][state].append(value)

states = ["Active promoter", "Weak promoter", "Poised promoter", "Strong enhancer", "Weak enhancer", "Tx elongation", "Weak Tx", "Insulator", "Polycomb repressed", "Heterochromatin", "Quiescence1", "Quiescence2"]

states_label = [state.split('_')[-1] for state in states]
for name in name_state_values:
    print name
    state_values = name_state_values[name]
    boxdata = [state_values[state] for state in states]
    fig = plt.figure(figsize=(8,6))
    plt.boxplot(boxdata, 0, '')
    #plt.ylabel(name.split('/')[-1], fontsize=15)
    plt.ylabel('Condensability (A.U.)', fontsize=15)
    plt.xticks(range(1, len(boxdata)+1), states_label, rotation=75, fontsize=15)
    #if name.startswith('work'):
    #    plt.ylim([-3.4, 2.8])
    #plt.savefig("box_" + 'HMM_'+name.split('/')[-1] + ".png",bbox_inches='tight')
    plt.savefig("box_" + 'HMM_'+name.split('/')[-1] + ".svg", format='svg', bbox_inches='tight')
    #plt.show()
    plt.close('all')

for name in name_state_values:
    print name
    state_values = name_state_values[name]
    data = [np.median(state_values[state]) for state in states]
    fig = plt.figure()
    plt.plot(range(1, len(data)+1), data, 'o-')
    plt.ylabel('Condensability (A.U.)', fontsize=15)
    plt.xticks(range(1, len(boxdata)+1), states_label, rotation=75, fontsize=15)
    plt.savefig("Mean_" + 'HMM_'+name.split('/')[-1] + ".svg", format='svg', bbox_inches='tight')
    plt.close('all')

