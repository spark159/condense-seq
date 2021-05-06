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

def state_cmp (a, b):
    num1 = int(a.split('_')[0])
    num2 = int(b.split('_')[0])
    if num1 < num2:
        return -1
    return 1

def tuple_cmp (a, b):
    if a[0] <= b[0]:
        return -1
    return 1
    
def read_bincountfile (fname, chr_list=None):
    First = True
    for line in open(fname):
        if First:
            cols = line.strip().split()
            names = [name.rsplit('.')[-2] for name in cols[4:-1]]
            chr_binID_counts = [{} for i in range(len(names))]
            chr_binID_range = {}
            chr_binID_GC = {}
            First = False
            continue
        line = line.strip()
        if not line:
            continue
        cols = line.strip().split()
        ID, chr, st, ed = cols[:4]
        if chr_list != None and chr not in chr_list:
            continue
        ID = int(ID)
        st, end = int(st), int(ed)
        GC = float(cols[-1])
        if chr not in chr_binID_range:
            chr_binID_range[chr] = []
        chr_binID_range[chr].append((st, ed))
        if chr not in chr_binID_GC:
            chr_binID_GC[chr] = []
        chr_binID_GC[chr].append(GC)
        datas = cols[4:-1]
        for i in range(len(datas)):
            data = float(datas[i])
            chr_binID_count = chr_binID_counts[i]
            if chr not in chr_binID_count:
                chr_binID_count[chr] = []
            chr_binID_count[chr].append(data)
    return names, chr_binID_counts, chr_binID_range, chr_binID_GC

def read_chromHMM(fname, chr_list=None, change=False):
    chr_state_intervals = {}
    for line in open(fname):
        cols = line.strip().split()
        chr, st, ed, state = cols[:4]
        if chr_list!=None and chr not in chr_list:
            continue
        if chr not in chr_state_intervals:
            chr_state_intervals[chr] = {}
        st, ed = int(st), int(ed)
        if change:
            state = change[state]
        if state not in chr_state_intervals[chr]:
            chr_state_intervals[chr][state] = []
        chr_state_intervals[chr][state].append((st,ed))
    return chr_state_intervals

def boxplot (key_values1, ylabel1='', title=None, keys=None, rotation=None, note=""):
    if keys:
        assert set(keys) <= set(key_values1.keys())
    else:
        keys = key_values1.keys()
    fig, ax1 = plt.subplots()
    pos_list1 = [i for i in range(len(keys))]
    bp1 = ax1.boxplot([key_values1[key] for key in keys], 0, "", positions=pos_list1, widths=0.3, patch_artist=True, boxprops=dict(facecolor="pink"))
    ax1.set_ylabel(ylabel1, color='r')
    ax1.tick_params('y', colors='r')
    #plt.legend([bp1["boxes"][0], [ylabel1], loc='upper right')
    if title:
        plt.title(title)
    ax1.set_xticks(range(len(keys)))
    ax1.set_xticklabels(keys, rotation=rotation)
    plt.xlim([-0.5, len(keys)-0.5])
    fig.tight_layout()
    #plt.savefig("box_" + note + ".png",bbox_inches='tight')
    plt.show()
    plt.close('all')


def pair_boxplot (key_values1, key_values2, ylabel1='', ylabel2='Condensability (A.U.)', title=None, keys=None, rotation=None, note=""):
    if keys:
        assert set(keys) <= set(key_values1.keys())
        assert set(keys) <= set(key_values2.keys())
    else:
        assert set(key_values1.keys()) == set(key_values2.keys())
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
    #plt.savefig("pairbox_" + note + ".png",bbox_inches='tight')
    plt.show()
    plt.close('all')

#fname = "NCP_Spermine(4+)_1kb"
#bin_size = 1000
#names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile("/home/spark159/Downloads/" + fname + "_bin.cn")
#chr_binID_control = chr_binID_counts[-1]

fname = "H1_NCP-new_spd_10kb_bin.cn"
bin_size = 10000
names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)
chr_binID_control = chr_binID_counts[-1]

#name_dict = {'E1':'TssBiv', 'E2':'TssA', 'E3':'EnhA', 'E4':'TxWk', 'E5':'Tx', 'E6':'me3Het', 'E7':'Quies', 'E8':'me2Het', 'E9':'PcWk', 'E10':'Pc'}
#chr_state_intervals = read_chromHMM("/home/spark159/../../media/spark159/sw/dataforcondense/38-Per_10_segments.bed", change=name_dict)

#chr_state_intervals = read_chromHMM("wgEncodeBroadHmmH1hescHMM.bed", change=False)
#chr_state_intervals = read_chromHMM("wgEncodeAwgSegmentationChromhmmH1hesc.bed", change=False)
chr_state_intervals = read_chromHMM("wgEncodeAwgSegmentationCombinedH1hesc.bed", change=False)

 
chr_binID_states = {}
for chr in chr_state_intervals:
    state_intervals = chr_state_intervals[chr]
    if chr not in chr_binID_states:
        chr_binID_states[chr] = {}
    for state in state_intervals:
        for interval in state_intervals[state]:
            st, ed = interval
            st_binID = int(st) / bin_size
            ed_binID = int(ed) / bin_size
            if st_binID == ed_binID:
                value = ed - st
                if st_binID not in chr_binID_states[chr]:
                    chr_binID_states[chr][st_binID] = []
                chr_binID_states[chr][st_binID].append((value, state))
            else:
                for k in range(st_binID, ed_binID+1):
                    if k == st_binID:
                        value = ((st_binID + 1) * bin_size) - st
                    elif k == ed_binID:
                        value = ed - ed_binID*bin_size
                    else:
                        value = bin_size
                    if k not in chr_binID_states[chr]:
                        chr_binID_states[chr][k] = []
                    chr_binID_states[chr][k].append((value, state))

state_GCs = {}
for chr in chr_binID_states:
    for binID in chr_binID_states[chr]:
        state = sorted(chr_binID_states[chr][binID], cmp=tuple_cmp, reverse=True)[0][1]
        if state not in state_GCs:
            state_GCs[state] = []
        GC = chr_binID_GC[chr][binID]
        state_GCs[state].append(GC*100)
                    
state_rcounts_list = []
for i in range(len(names)-1):
    chr_binID_count = chr_binID_counts[i]
    state_rcounts = {}
    for chr in chr_binID_states:
        for binID in chr_binID_states[chr]:
            state = sorted(chr_binID_states[chr][binID], cmp=tuple_cmp, reverse=True)[0][1]
            if state not in state_rcounts:
                state_rcounts[state] = []
            control = chr_binID_control[chr][binID]
            if control <= 0:
                continue
            test = chr_binID_count[chr][binID]
            rcount = float(test) / control
            #rcount = float(test)
            state_rcounts[state].append(rcount)
    state_rcounts_list.append(state_rcounts)

#states = ['TssA', 'EnhA', 'Tx', 'TxWk', 'TssBiv', 'me2Het', 'me3Het', 'PcWk', 'Pc', 'Quies']
#states = sorted(state_rcounts_list[0].keys(), cmp=state_cmp)
states = sorted(state_rcounts_list[0].keys())
for i in range(len(state_rcounts_list)):
    name = names[i]
    state_rcounts = state_rcounts_list[i]
    pair_boxplot (state_GCs, state_rcounts, ylabel1='GC content(%)', ylabel2='Normalized Counts', title='Titration ' +str(i+1), keys=states, note='HMM_bin_' + fname + '_' + str(i+1), rotation=75)
    #pair_boxplot (state_GCs, state_rcounts, ylabel1='GC content(%)', ylabel2='Normalized Counts', title='Titration ' +str(i+1), keys=None, note='HMM_bin_' + fname + '_' + str(i+1), rotation=75)

"""
dID_interval = {}
for state in state_intervals:
    intervals = state_intervals[state]
    for i in range(len(intervals)):
        dID = state + ':' + str(i)
        assert dID not in dID_interval
        dID_interval[dID] = intervals[i]

dinterval_dict = Interval_dict.double_hash(dID_interval, 10000, 250000000)

ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("data/hg19_chr1_171_everything_anot.cn")

name_ID_value['sp10-ATcontent'] = statis.neutralize_score_by_target(name_ID_value['work/condense_seq/sp10_hg19_chr1'], name_ID_value['ATcontent'])

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

name_ID_value['sp10-ATcontent'] = statis.neutralize_score_by_target(name_ID_value['work/condense_seq/sp10_hg19_chr1'], name_ID_value['ATcontent'])

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
            state_name_values[state][name].append(value)
            name_state_values[name][state].append(value)

#states = ['1_TssA', '2_TssAFlnk', '3_TxFlnk', '4_Tx', '5_TxWk', '6_EnhG', '7_Enh', '8_ZNF/Rpts', '9_Het', '10_TssBiv', '11_BivFlnk', '12_EnhBiv', '13_ReprPC', '14_ReprPCWk', '15_Quies']
#states = ['TssA', 'Enh', 'Tx', 'me2Het', 'me3Het', 'Polycomb', 'Quies']
states = ['TssA', 'EnhA', 'Tx', 'TxWk', 'TssBiv', 'me2Het', 'me3Het', 'PcWk', 'Pc', 'Quies']


states_label = [state.split('_')[-1] for state in states]
state_scores2 = name_state_values['work/condense_seq/sp10_hg19_chr1']
for name in name_state_values:
    state_values = name_state_values[name]
    boxdata = [state_values[state] for state in states]
    fig = plt.figure()
    plt.boxplot(boxdata, 0, '')
    plt.ylabel(name.split('/')[-1])
    plt.xticks(range(1, len(boxdata)+1), states_label, rotation=75)
    plt.savefig("box_" + 'HMM_'+name.split('/')[-1] + ".png",bbox_inches='tight')
    #plt.show()
    plt.close('all')
    if name == 'work/condense_seq/sp10_hg19_chr1':
        continue
    pair_boxplot (state_values, state_scores2, ylabel1=name.split('/')[-1], ylabel2='Condensability (A.U.)', title=None, keys=states, note='HMM_'+name.split('/')[-1], rotation=75)
"""
