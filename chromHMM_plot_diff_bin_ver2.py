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

# chromHMM data set
cell_chromHMM = {'H1':[H1_HMM_fname, H1_name_dict, H1_states],
                 'GM':[GM_HMM_fname, GM_name_dict, GM_states],
                 'mCD8T':[mCD8T_HMM_fname, mCD8T_name_dict, mCD8T_states]}

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

# read bin score file
def read_bin_score (fname):
    name_chr_binID_score = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        if First:
            names = cols[4:]
            First = False
            continue
        _, chr, st, ed = cols[:4]
        st, ed = int(st), int(ed)
        scores = [float(score) for score in cols[4:]]
        binID = (st, ed)
        for name, score in zip(names, scores):
            if name not in name_chr_binID_score:
                name_chr_binID_score[name] = {}
            if chr not in name_chr_binID_score[name]:
                name_chr_binID_score[name][chr] = {}
            name_chr_binID_score[name][chr][binID] = score
    return name_chr_binID_score

# read bin score file
def read_bin_score_new (fname, bin_size):
    name_chr_binID_score = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        if First:
            names = cols[4:]
            First = False
            continue
        _, chr, st, ed = cols[:4]
        st, ed = int(st), int(ed)
        scores = [float(score) for score in cols[4:]]
        binID = st / int(bin_size)
        for name, score in zip(names, scores):
            if name not in name_chr_binID_score:
                name_chr_binID_score[name] = {}
            if chr not in name_chr_binID_score[name]:
                name_chr_binID_score[name][chr] = {}
            name_chr_binID_score[name][chr][binID] = score
    return name_chr_binID_score


# read bin Chalf
def read_bin_Chalf (fname, chr_choices=None):
    ID_pos = {}
    ID_value = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        if First:
            First = False
            continue
        binID, chr, st, ed, Chalf = cols
        if chr_choices != None and chr not in chr_choices:
            continue
        st, ed = int(st), int(ed)
        ID = (st, ed)
        st, ed = int(st), int(ed)
        Chalf = float(Chalf)
        pos = int(float(st + ed)/2)
        ID_pos[ID] = pos
        ID_value[ID] = Chalf
    return ID_pos, ID_value


def read_chromHMM(fname, chr_choice, change=False):
    chr_state_intervals = {}
    for line in open(fname):
        cols = line.strip().split()
        chr, st, ed, state = cols[:4]
        if chr not in chr_choice:
            continue
        st, ed = int(st), int(ed)
        if change:
            state = change[state]
        if chr not in chr_state_intervals:
            chr_state_intervals[chr] = {}
        if state not in chr_state_intervals[chr]:
            chr_state_intervals[chr][state] = []
        chr_state_intervals[chr][state].append((st,ed))
    return chr_state_intervals

def single_boxplot (key_values, ylabel='', title=None, keys=None,
                    rotation=None, ylim=[None, None], note=""):

    if keys:
        new_keys = []
        for key in keys:
            if key in key_values:
                new_keys.append(key)
        keys = new_keys

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

def multi_boxplot (key_values_list,
                   axhline=None,
                   ylabel='Condensability (A.U.)',
                   title=None,
                   keys=None,
                   labels=[],
                   colors=[],
                   rotation=None,
                   legend_loc='best',
                   note=""):

    common_keys = set([])
    for i in range(len(key_values_list)):
        key_values = key_values_list[i]
        if i == 0:
            common_keys |= set(key_values.keys())
            continue
        common_keys &= set(key_values.keys())
    common_keys = list(common_keys)

    if keys:
        new_keys = []
        for key in keys:
            if key in common_keys:
                new_keys.append(key)
        keys = new_keys
    else:
        keys = common_keys

    if not labels:
        labels = [None]*len(key_values_list)

    if not colors:
        colors = ['white']*len(key_values_list)

    offset = 0.6/len(key_values_list)
    
    bp_list = []
    fig = plt.figure()
    if axhline != None:
        plt.axhline(y=axhline, linestyle='--', color='k', alpha=0.5)
    for i in range(len(key_values_list)):
        key_values = key_values_list[i]
        pos_list = [k + offset*i for k in range(len(keys))]
        bp = plt.boxplot([key_values[key] for key in keys], positions=pos_list, showfliers=False,
                         notch=True, widths=offset*0.7, patch_artist=True,
                         boxprops=dict(facecolor=colors[i]))
        bp_list.append(bp)

    for bp in bp_list:
        for median in bp['medians']:
            median.set_color('red')
    xtick_locs =[k + 0.5*offset*(len(key_values_list)-1) for k in range(len(keys))]
    plt.xticks(xtick_locs, keys, rotation=rotation)
    plt.xlim([xtick_locs[0]-offset*3, xtick_locs[-1]+offset*3])
    plt.legend([bp["boxes"][0] for bp in bp_list], labels, loc=legend_loc)
    plt.ylabel(ylabel)
    if title:
        plt.title(title)
    plt.tight_layout()
    #plt.savefig("multibox_" + note + ".png",bbox_inches='tight')
    plt.savefig("multibox_" + note + ".svg", format='svg', bbox_inches='tight')
    #plt.show()
    plt.close('all')

def multi_boxplot2 (key_values_list,
                    figsize=None,
                    axhline=None,
                    ylabel='Condensability (A.U.)',
                    ylims=[None, None],
                    title=None,
                    keys=None,
                    labels=[],
                    colors=[],
                    rotation=None,
                    legend_loc='best',
                    note=""):

    common_keys = set([])
    for i in range(len(key_values_list)):
        key_values = key_values_list[i]
        if i == 0:
            common_keys |= set(key_values.keys())
            continue
        common_keys &= set(key_values.keys())
    common_keys = list(common_keys)

    if keys:
        new_keys = []
        for key in keys:
            if key in common_keys:
                new_keys.append(key)
        keys = new_keys
    else:
        keys = common_keys

    if not labels:
        labels = [None]*len(key_values_list)

    if not colors:
        colors = ['white']*len(key_values_list)

    bp_list = []
    #width = 0.2*len(keys)
    #height = 0.2*len(key_values_list)

    if figsize !=None:
        figsize = tuple(figsize)
    
    fig, axes = plt.subplots(figsize=figsize, nrows=len(key_values_list), ncols=1)
    for i in range(len(key_values_list)):
        if axhline != None:
            axes[i].axhline(y=axhline, linestyle='--', color='k', alpha=0.5)

        key_values = key_values_list[i]
        bp = axes[i].boxplot([key_values[key] for key in keys],
                             positions=range(len(keys)),
                             showfliers=False,
                             notch=True,
                             widths=0.3,
                             patch_artist=True,
                             boxprops=dict(facecolor=colors[i]))
        bp_list.append(bp)

        axes[i].set_xlim([-0.5, len(keys)-0.5])
        axes[i].set_ylim(ylims)
        axes[i].set_ylabel(ylabel)

        axes[i].spines['top'].set_visible(False)
        axes[i].spines['bottom'].set_visible(False)
        axes[i].spines['left'].set_visible(True)
        axes[i].spines['right'].set_visible(False)

        if i < len(key_values_list) - 1:        
            axes[i].set_xticks([])
            axes[i].set_xticklabels([])
                            
        else:
            assert i == len(key_values_list) - 1
            axes[i].spines['bottom'].set_visible(True)
            axes[i].set_xticks(range(len(keys)))
            axes[i].set_xticklabels(keys, rotation=rotation, ha="right", va='center', rotation_mode='anchor')
    
    for bp in bp_list:
        for median in bp['medians']:
            median.set_color('red')
            
    plt.legend([bp["boxes"][0] for bp in bp_list], labels, loc=legend_loc)
    if title:
        plt.title(title)
    plt.tight_layout()
    plt.savefig("multibox2_" + note + ".svg", format='svg', bbox_inches='tight')
    plt.close('all')

def standardization (ID_value):
    ID_newvalue = {}
    mean = np.mean(ID_value.values())
    std = np.std(ID_value.values())

    for ID, value in ID_value.items():
        newvalue = float(value - mean) / std
        ID_newvalue[ID] = newvalue
    return ID_newvalue

def categorize (chr_state_intervals, chr_binID_value, chr_choice=None):

    if chr_choice == None:
        chr_choice = list(set(chr_state_intervals.keys()) & set(chr_binID_value.keys()))
    
    state_values = {}
    for chr in chr_choice:
        try:
            state_intervals = chr_state_intervals[chr]
            binID_value = chr_binID_value[chr]
        except:
            continue

        print "processing %s" % (chr)

        dID_interval = {}
        for state in state_intervals:
            intervals = state_intervals[state]
            for i in range(len(intervals)):
                dID = state + ':' + str(i)
                assert dID not in dID_interval
                dID_interval[dID] = intervals[i]

        dinterval_dict = Interval_dict.double_hash(dID_interval, 10000, 250000000)

        for binID in binID_value:
            bst, bed = binID
            dIDs = dinterval_dict.find_range(bst, bed)
            if not dIDs:
                continue
            overlap_dID = []
            for dID in dIDs:
                dst, ded = dinterval_dict.ID_interval[dID]
                a, b = max(dst, bst), min(ded, bed)
                length = b - a
                overlap_dID.append((length, dID))
            overlap_dID = sorted(overlap_dID, reverse=True)
            best_dID = overlap_dID[0][1]
            state, _ = best_dID.split(':')
            if state not in state_values:
                state_values[state] = []
            value = binID_value[binID]
            if np.isnan(value):
                continue
            state_values[state].append(value)

    return state_values


def categorize_bin (chr_state_intervals, chr_binID_value, bin_size, chr_choice=None):

    if chr_choice == None:
        chr_choice = list(set(chr_state_intervals.keys()) & set(chr_binID_value.keys()))
    
    chr_binID_states = {}
    for chr in chr_choice:
        try:
            state_intervals = chr_state_intervals[chr]
            binID_value = chr_binID_value[chr]
        except:
            continue

        for state in state_intervals:
            for interval in state_intervals[state]:
                st, ed = interval
                st_binID = int(st) / bin_size
                ed_binID = int(ed) / bin_size
                if st_binID == ed_binID:
                    value = ed - st
                    if chr not in chr_binID_states:
                        chr_binID_states[chr] = {}
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
                        if chr not in chr_binID_states:
                            chr_binID_states[chr] = {}
                        if k not in chr_binID_states[chr]:
                            chr_binID_states[chr][k] = []
                        chr_binID_states[chr][k].append((value, state))

    state_values = {}
    for chr in chr_binID_states:
        for binID in chr_binID_states[chr]:
            state = sorted(chr_binID_states[chr][binID], reverse=True)[0][1]
            if state not in state_values:
                state_values[state] = []
            try:
                value = chr_binID_value[chr][binID]
                state_values[state].append(value)
            except:
                pass

    return state_values


#path = "/home/spark159/../../media/spark159/sw/"
path = "/home/spark159/../../storage/"

agent_fullname = {'sp':'Spermine(4+)',
                  'spd':'Spermidine(3+)',
                  'CoH':'Cobalt Hexammine(3+)',
                  'PEG':'PEG 8000',
                  'HP1a':'HP1$\\alpha$',
                  'HP1bSUV':'HP1$\\beta$+tSUV',
                  'LKH':'Linker histone1',
                  'Ki67':'Ki67',
                  'FUS':'FUS',
                  'Mg':'Magnesium',
                  'Ca':'Calcium'}

# experiment list (cell, sample, agent, tnum)
# it should be same cell line
cell = 'H1'
#exp_list = [(cell, 'NCP', 'sp'),
#            (cell, 'NCP', 'spd'),
#            (cell, 'NCP', 'CoH'),
#            (cell, 'NCP', 'PEG'),
#            (cell, 'NCP', 'Ca'),
#            (cell, 'NCP', 'Mg'),
#            (cell, 'NCP', 'HP1a'),
#            (cell, 'NCP', 'HP1bSUV'),
#            (cell, 'NCP', 'LKH'),
#            (cell, 'NCP', 'Ki67'),
#            (cell, 'NCP', 'FUS')]

#exp_list = [(cell, 'NCP', 'sp'),
#            (cell, 'NCP', 'spd')]

#exp_list = [(cell, 'NCP', 'sp', 8),
#            (cell, 'NCP', 'HP1a', 3)]

#exp_list = [(cell, 'NCP', 'sp', 8)]
#exp_list = [(cell, 'NCP', 'HP1a', 3)]

#exp_list = [(cell, 'NCP', 'sp', 8),
#            (cell, 'NCP', 'spd', 6),
#            (cell, 'NCP', 'CoH', 5),
#            (cell, 'NCP', 'PEG', 6),
#            (cell, 'NCP', 'Ca', 5),
#            (cell, 'NCP', 'Mg', 5),
#            (cell, 'NCP', 'HP1a', 3),
#            (cell, 'NCP', 'HP1bSUV', 3),
#            (cell, 'NCP', 'LKH', 1),
#            (cell, 'NCP', 'Ki67', 1),
#            (cell, 'NCP', 'FUS', 1)]

#exp_list = [(cell, 'NCP', 'sp', 8),
#            (cell, 'NCP', 'spd', 6),
#            (cell, 'NCP', 'CoH', 5),
#            (cell, 'NCP', 'PEG', 6),
#            (cell, 'NCP', 'Ca', 5),
#            (cell, 'NCP', 'Mg', 5),
#            (cell, 'NCP', 'HP1a', 3),
#            (cell, 'NCP', 'HP1bSUV', 3)]

#exp_list = [(cell, 'NCP', 'sp', 8),
#            (cell, 'NCP', 'HP1a', 3)]

#exp_list = [(cell, 'NCP', 'sp', 8),
#            (cell, 'NCP', 'spd', 6),
#            (cell, 'NCP', 'CoH', 5),
#            (cell, 'NCP', 'PEG', 6),
#            (cell, 'NCP', 'Ca', 5),
#            (cell, 'NCP', 'Mg', 5),
#            (cell, 'NCP', 'HP1a', 3),
#            (cell, 'NCP', 'HP1bSUV', 4),
#            (cell, 'NCP', 'LKH', 3),
#            (cell, 'NCP', 'Ki67', 4),
#            (cell, 'NCP', 'FUS', 5)]

#exp_list = [(cell, 'NCP', 'sp', 8),
#            (cell, 'NCP', 'HP1a', 3),
#            (cell, 'NCP', 'LKH', 3),
#            (cell, 'NCP', 'Ki67', 4)]

exp_list = [(cell, 'NCP', 'sp', 8),
            (cell, 'NCP', 'spd', 6),
            (cell, 'NCP', 'CoH', 5),
            (cell, 'NCP', 'PEG', 6),
            (cell, 'NCP', 'Ca', 5),
            (cell, 'NCP', 'HP1a', 3),
            (cell, 'NCP', 'HP1bSUV', 4)]


# binsize of input data
#bin_size = 10000
bin_size = 5000
#bin_size = 1000

# set species and gender
if cell in ['H1', 'GM']:
    species = 'human'
elif cell in ['mCD8T']:
    species = 'mouse'

if cell in ['H1']:
    gender = 'male'
elif cell in ['GM', 'mCD8T']:
    gender = 'female'

# set chromosome list
if species == 'human':
    chr_list = ['chr' + str(i) for i in range(1, 23)]
elif species == 'mouse':
    chr_list = ['chr' + str(i) for i in range(1, 20)]
chr_list += ['chrX']

if gender == 'male':
    chr_list += ['chrY']

#chr_list = ['chr1']

# other parameters
dtype = 'zscore'
note = 'HMM_diff_%skb' % (int(bin_size/1000.0))
ylabel = 'zscore'
labels = [agent_fullname[agent] for cell, sample, agent, tnum in exp_list]
#colors = [None for i in range(len(exp_list))]
colors = ['tab:blue', 'tab:red']
colors = ['tab:blue', 'tab:red', 'tab:green', 'tab:purple'] * 10


# load chromHMM genome segmentation file
HMM_fname, name_dict, states = cell_chromHMM[cell]
chr_state_intervals = read_chromHMM(HMM_fname, chr_choice=chr_list, change=name_dict)


# load data and sort them according to chromatin-state
state_scores_list = []
for cell, sample, agent, tnum in exp_list:
    fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', dtype]) + '.cn'

    if dtype in ['Chalf', 'zChalf']:
        field_name = 'Chalf'
    else:
        field_name = "%s-%s-%s-%d.bam" % (cell, sample, agent, tnum)
    

    chr_binID_score = read_bin_score_new(fname, bin_size)[field_name]
    state_scores = categorize_bin (chr_state_intervals, chr_binID_score, bin_size)

    single_boxplot (state_scores,
                    ylabel=ylabel,
                    title=None,
                    keys=states,
                    rotation=75,
                    ylim=[None, None],
                    note=note + "%s_%s_%s" % (cell, sample, agent))

    state_scores_list.append(state_scores)

    
# plot box plot
multi_boxplot2 (state_scores_list,
                figsize=(3.5, 6),
                axhline=0,
                ylabel=ylabel,
                ylims=[-4, 4],
                title=None,
                keys=states,
                labels = None, 
                colors = colors, 
                rotation=75,
                note=note + '_' + ylabel)

#multi_boxplot (state_scores_list,
#               ylabel=ylabel,
#               title=None,
#               keys=states,
#               labels = labels, 
#               colors = colors, 
#               rotation=75,
#               note=note + '_' + ylabel)


# plot mean score heatmap
states.remove('Insulator')
img = []
for state_scores in state_scores_list:
    row = []
    for state in states:
        scores = state_scores[state]
        mean = np.mean(scores)
        row.append(mean)
    img.append(row)


vmin = -2.5
vmax = 2.5
cmap = 'jet_r'
#cmap = 'bwr_r'
#cmap = 'Spectral'
#cmap = 'seismic_r'
width = 0.4*len(states)
height = 0.4*len(exp_list)
fig = plt.figure(figsize=(width, height))
im = plt.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)
plt.xticks(range(len(states)), states, ha='right', va='center', rotation_mode='anchor', rotation=45)
plt.yticks(range(len(exp_list)), [exp[2]  for exp in exp_list])
plt.savefig("chrom_meanscore.svg", format='svg', bbox_inches='tight')
plt.close()

# plot colorbar only
fig = plt.figure(figsize=(1.2,1))
plt.subplot(1,2,1)
cbar = plt.colorbar(im, cax=plt.gca(), ticks=[vmin, vmax])
cbar.ax.set_yticklabels([vmin, vmax], fontsize=8)
cbar.ax.set_ylabel('Mean zscore', rotation=-90, va="bottom", fontsize=8)
plt.tight_layout()
plt.savefig('chrom_mean_cbar.svg', format='svg', bbox_inches='tight')
plt.close()

