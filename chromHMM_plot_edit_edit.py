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

## for H1
H1_name_dict = {"E1":"Polycomb repressed",
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

# state for H1
H1_states = ["Active promoter", "Weak promoter", "Poised promoter", "Strong enhancer", "Weak enhancer", "Tx elongation", "Weak Tx", "Insulator", "Polycomb repressed", "Heterochromatin", "Quiescence1", "Quiescence2"]

H1_HMM_fname = "H1_12_segments.bed"

# for GM12878
GM_name_dict = {"E1":"Polycomb repressed",
             "E2":"Quiescence",
             "E3":"Heterochromatin",
             "E4":"Weak Tx",
             "E5":"Tx elongation",
             "E6":"Weak enhancer",
             "E7":"Active enhancer",
             "E8":"Strong enhancer",
             "E9":"Active promoter",
             "E10":"Weak promoter",
             "E11":"Poised promoter",
             "E12":"Insulator"}

# state for GM
GM_states = ["Active promoter", "Weak promoter", "Poised promoter", "Strong enhancer", "Active enhancer", "Weak enhancer", "Tx elongation", "Weak Tx", "Insulator", "Polycomb repressed", "Heterochromatin", "Quiescence"]

GM_HMM_fname = "GM12878_12_segments.bed"

# for mouse CD8 T cell
mCD8T_name_dict = {"E1":"Weak Tx",
             "E2":"Tx elongation",
             "E3":"Weak enhancer2",
             "E4":"Strong enhancer2",
             "E5":"Strong enhancer1",
             "E6":"Weak enhancer1",
             "E7":"Active promoter",
             "E8":"Poised promoter",
             "E9":"Polycomb repressed1",
             "E10":"Polycomb repressed2",
             "E11":"Quiescence",
             "E12":"Heterochromatin"}

# state for mouse CD8 T cell
mCD8T_states = ["Active promoter", "Poised promoter", "Strong enhancer1", "Strong enhancer2", "Weak enhancer1", "Weak enhancer2", "Tx elongation", "Weak Tx", "Polycomb repressed1", "Polycomb repressed2", "Heterochromatin", "Quiescence"]

mCD8T_HMM_fname = "Mouse CD8 T cell (invitro activated)_12_segments.bed"

# read num file
def read_num_file (fname, chr_choice=None):
    ID_pos = {}
    ID_score1, ID_score2 = {}, {}
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        if cols[-1] == '*':
            #print "here"
            continue
        chr = cols[1]
        if chr_choice and cols[1] != chr_choice:
            continue
        #ID = chr + ':' + cols[0]
        ID = cols[0]
        pos = int(cols[2])
        sup1 = float(cols[3])
        sup2 = float(cols[4])
        control = float(cols[5])
        if sup1 * sup2 * control <= 0:
            continue
        score1 = -np.log(sup1/control)
        score2 = -np.log(sup2/control)
        ID_pos[ID] = pos
        ID_score1[ID] = score1
        ID_score2[ID] = score2
    return ID_pos, ID_score1, ID_score2


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

def single_boxplot (key_values, ylabel='', title=None, keys=None,
                    rotation=None, ylim=[None, None], note=""):
    if keys:
        assert set(keys) <= set(key_values.keys())
    else:
        keys = key_values.keys()
    fig, ax1 = plt.subplots()
    pos_list = [i for i in range(len(keys))]
    bp1 = ax1.boxplot([key_values[key] for key in keys],
                      0, "",
                      positions=pos_list,
                      widths=0.5, patch_artist=True, boxprops=dict(facecolor="pink"))
    ax1.set_ylabel(ylabel, color='r')
    ax1.tick_params('y', colors='r')
    if title:
        plt.title(title)
    ax1.set_xticks(range(len(keys)))
    ax1.set_xticklabels(keys, rotation=rotation)
    plt.xlim([-0.5, len(keys)-0.5])
    plt.ylim(ylim)
    fig.tight_layout()
    plt.savefig("singlebox_" + note + ".png",bbox_inches='tight')
    #plt.show()
    plt.close()


def pair_boxplot (key_values1, key_values2, ylabel1='', ylabel2='Condensability (A.U.)', title=None, keys=None, rotation=None, sharey = False, note=""):
    if keys:
        assert set(keys) <= set(key_values1.keys())
        assert set(keys) <= set(key_values2.keys())
    else:
        assert set(key_values1,keys()) == set(key_values2.keys())
        keys = key_values1.keys()
    fig, ax1 = plt.subplots()
    pos_list1 = [i - 0.2 for i in range(len(keys))]
    bp1 = ax1.boxplot([key_values1[key] for key in keys], 0, "", positions=pos_list1, widths=0.3, patch_artist=True, boxprops=dict(facecolor="pink"))
    ax1.set_ylabel(ylabel1, color='r')
    ax1.tick_params('y', colors='r')
    if sharey:
        ax2 = ax1
    else:
        ax2 = ax1.twinx()
    pos_list2 = [i + 0.2 for i in range(len(keys))]
    bp2 = ax2.boxplot([key_values2[key] for key in keys], 0, "", positions=pos_list2, widths=0.3, patch_artist=True, boxprops=dict(facecolor="lightblue"))
    ax2.set_ylabel(ylabel2, color='b')
    ax2.tick_params('y', colors='b')
    plt.legend([bp1["boxes"][0], bp2["boxes"][0]], [ylabel1, ylabel2], loc='upper right')
    if title:
        plt.title(title)
    ax1.set_xticks(range(len(keys)))
    ax1.set_xticklabels(keys, rotation=rotation)
    ax2.set_xticks(range(len(keys)))
    ax2.set_xticklabels(keys, rotation=rotation)
    plt.xlim([-0.5, len(keys)-0.5])
    fig.tight_layout()
    plt.savefig("pairbox_" + note + ".png",bbox_inches='tight')
    #plt.show()
    plt.close('all')

def standardization (ID_value):
    ID_newvalue = {}
    mean = np.mean(ID_value.values())
    std = np.std(ID_value.values())

    for ID, value in ID_value.items():
        newvalue = float(value - mean) / std
        ID_newvalue[ID] = newvalue
    return ID_newvalue

def categorize (state_intervals, ID_pos, ID_value):
    
    dID_interval = {}
    for state in state_intervals:
        intervals = state_intervals[state]
        for i in range(len(intervals)):
            dID = state + ':' + str(i)
            assert dID not in dID_interval
            dID_interval[dID] = intervals[i]

    dinterval_dict = Interval_dict.double_hash(dID_interval, 10000, 250000000)

    state_values = {}
    for ID in ID_pos:
        pos = ID_pos[ID]
        dIDs = dinterval_dict.find(pos)
        if not dIDs:
            continue
        for dID in dIDs:
            state, _ = dID.split(':')
            if state not in state_values:
                state_values[state] = []
            value = ID_value[ID]
            if np.isnan(value):
                continue
            state_values[state].append(value)

    return state_values

path = "/home/spark159/../../media/spark159/sw/"

#cell1, cell2 = "H1", "GM"
cell1, cell2 = "mCD8T", "mCD8T"
#sample1, sample2 = "NCP", "NCP"
sample1, sample2 = "WT-NCP", "inht-NCP"
agent = "sp"
chr_choice = 'chr1'

num_fname1 = "_".join([cell1, sample1, agent, chr_choice]) + "_num.cn"
num_fname2 = "_".join([cell2, sample2, agent, chr_choice]) + "_num.cn"

ID_pos1, _ , ID_score1 = read_num_file(path+num_fname1)
ID_pos2, _ , ID_score2 = read_num_file(path+num_fname2)

ID_score1 = standardization(ID_score1)
ID_score2 = standardization(ID_score2)

#state_intervals1 = read_chromHMM(H1_HMM_fname, chr_target=chr_choice, change=H1_name_dict)
#state_intervals2 = read_chromHMM(GM_HMM_fname, chr_target=chr_choice, change=GM_name_dict)
state_intervals = read_chromHMM(mCD8T_HMM_fname, chr_target=chr_choice, change=mCD8T_name_dict)

state_scores1 = categorize (state_intervals, ID_pos1, ID_score1)
state_scores2 = categorize (state_intervals, ID_pos2, ID_score2)

pair_boxplot (state_scores1, state_scores2, ylabel1="WT", ylabel2="+ODCiht", title=None, keys=mCD8T_states, note='HMM_diff_', rotation=75, sharey=True)
#single_boxplot (state_scores1, keys=H1_states, rotation=75, ylim=[0, 8], note='H1')
#single_boxplot (state_scores2, keys=GM_states, rotation=75, ylim=[0, 8], note='GM')




"""
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

#states = ['TssA', 'EnhA', 'Tx', 'TxWk', 'TssBiv', 'me2Het', 'me3Het', 'PcWk', 'Pc', 'Quies']
#states = state_intervals.keys()

# state for H1
#states = ["Active promoter", "Weak promoter", "Poised promoter", "Strong enhancer", "Weak enhancer", "Tx elongation", "Weak Tx", "Insulator", "Polycomb repressed", "Heterochromatin", "Quiescence1", "Quiescence2"]

# state for GM
#states = ["Active promoter", "Weak promoter", "Poised promoter", "Strong enhancer", "Active enhancer", "Weak enhancer", "Tx elongation", "Weak Tx", "Insulator", "Polycomb repressed", "Heterochromatin", "Quiescence"]

# state for mouse CD8 T cell
states = ["Active promoter", "Poised promoter", "Strong enhancer1", "Strong enhancer2", "Weak enhancer1", "Weak enhancer2", "Tx elongation", "Weak Tx", "Polycomb repressed1", "Polycomb repressed2", "Heterochromatin", "Quiescence"]



states_label = [state.split('_')[-1] for state in states]
#state_scores1 = name_state_values['/home/spark159/scratch4-tha4/sangwoo/2022_09_08_GM_sp_H1_HP1a_deep/GM-NCP-sp-4']
#state_scores2 = name_state_values['/home/spark159/scratch4-tha4/sangwoo/2022_09_08_GM_sp_H1_HP1a_deep/GM-NCP-sp-8']
#for name in ["/home/spark159/scratch4-tha4/sangwoo/MouseCD8Tcell_detail/mCD8T-inht-NCP-sp-8"]:
for name in name_state_values:
    print name
    state_values = name_state_values[name]
    boxdata = [state_values[state] for state in states]
    fig = plt.figure(figsize=(8,6))
    plt.boxplot(boxdata, 0, '')
    #plt.ylabel(name.split('/')[-1], fontsize=15)
    plt.ylabel('Condensability (A.U.)', fontsize=15)
    plt.xticks(range(1, len(boxdata)+1), states_label, rotation=75, fontsize=15)
    if name.startswith('work'):
        plt.ylim([-3.4, 2.8])
    #plt.savefig("box_" + 'HMM_'+name.split('/')[-1] + ".png",bbox_inches='tight')
    plt.savefig("box_" + 'HMM_'+name.split('/')[-1] + ".svg", format='svg', bbox_inches='tight')
    #plt.show()
    plt.close('all')
    if name in ['work/2021_06_07_H1_sp_detail/H1-NCP-sp-4', 'work/2021_06_07_H1_sp_detail/H1-NCP-sp-8']:
        continue
    #pair_boxplot (state_values, state_scores1, ylabel1=name.split('/')[-1], ylabel2='Condensability (A.U.)', title=None, keys=states, note='HMM_'+name.split('/')[-1], rotation=75)

fig = plt.figure(figsize=(8,6))
plt.boxplot(boxdata, 0, '')
#plt.ylabel(name.split('/')[-1], fontsize=15)
plt.ylabel('Condensability (A.U.)', fontsize=15)
plt.gca().tick_params(axis='both', which='major', labelsize=15)
plt.gca().tick_params(axis='both', which='minor', labelsize=15)
plt.xticks(range(1, len(boxdata)+1), states_label, rotation=75, fontsize=15)
if name.startswith('work'):
    plt.ylim([-3.4, 2.8])
#plt.savefig("box_" + 'HMM_'+name.split('/')[-1] + ".png",bbox_inches='tight')
plt.savefig("box_" + 'HMM_'+name.split('/')[-1] + ".svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close('all')
"""
