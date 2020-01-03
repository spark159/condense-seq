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

def pair_boxplot (key_values1, key_values2, ylabel1='', ylabel2='Condensability (A.U.)', title=None, keys=None, rotation=None, note=""):
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

path = "/home/spark159/../../media/spark159/sw/sp_spd_tests_detail/"

ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path + "hg19_chr1_167win25step_anot.cn", jump=10)
#ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path + "hg19_chr1_NCP_anot.cn")

ID_meGC = name_ID_value['meGCNumber']
ID_CpG = name_ID_value['CpGNumber']
ID_meCpG = {}
for ID in ID_meGC:
    meGC = ID_meGC[ID]
    CpG = ID_CpG[ID]
    if CpG <= 0:
        meCpG = 0.0
    else:
        meCpG = float(meGC) / CpG
    ID_meCpG[ID] = meCpG
name_ID_value['meCpGfrac'] = ID_meCpG

#name_dict = {'E1':'TssBiv', 'E2':'TssA', 'E3':'EnhA', 'E4':'TxWk', 'E5':'Tx', 'E6':'me3Het', 'E7':'Quies', 'E8':'me2Het', 'E9':'PcWk', 'E10':'Pc'}
#state_intervals = read_chromHMM("/home/spark159/../../media/spark159/sw/dataforcondense/38-Per_10_segments.bed", chr_target='chr1', change=name_dict)

state_intervals = read_chromHMM("/home/spark159/../../media/spark159/sw/dataforcondense/wgEncodeAwgSegmentationCombinedGm12878.bed", chr_target='chr1', change=False)

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

states = ['TssA', 'EnhA', 'Tx', 'TxWk', 'TssBiv', 'me2Het', 'me3Het', 'PcWk', 'Pc', 'Quies']
states = state_intervals.keys()

states_label = [state.split('_')[-1] for state in states]
state_scores1 = name_state_values['data/sp_spd_tests_detail/sp7']
state_scores2 = name_state_values['data/sp_spd_tests_detail/sp8']
for name in name_state_values:
    print name
    state_values = name_state_values[name]
    boxdata = [state_values[state] for state in states]
    fig = plt.figure()
    plt.boxplot(boxdata, 0, '')
    plt.ylabel(name.split('/')[-1])
    plt.xticks(range(1, len(boxdata)+1), states_label, rotation=75)
    plt.savefig("box_" + 'HMM_'+name.split('/')[-1] + ".png",bbox_inches='tight')
    plt.show()
    plt.close('all')
    if name in ['data/sp_spd_tests_detail/sp7', 'data/sp_spd_tests_detail/sp8']:
        continue
    pair_boxplot (state_values, state_scores1, ylabel1=name.split('/')[-1], ylabel2='Condensability (A.U.)', title=None, keys=states, note='HMM_'+name.split('/')[-1], rotation=75)
