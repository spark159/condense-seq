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


def read_bin_file (fname, chr_choices=None):
    First = True
    ID_pos = {}
    name_ID_value = {}
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        if First:
            names = []
            for name in cols[4:]:
                names.append(name.rsplit('.', 1)[0])
            First = False
            continue
        ID, chr, start, end = cols[:4]
        start, end = int(start), int(end)
        if chr_choices!=None and chr not in chr_choices:
            continue
        pos = int(float(start + end)/2)
        ID = (chr, start, end)
        assert ID not in ID_pos
        ID_pos[ID] = pos
        values = [float(value) for value in cols[4:]]
        for name, value in zip(names, values):
            if name not in name_ID_value:
                name_ID_value[name] = {}
            assert ID not in name_ID_value[name]
            name_ID_value[name][ID] = value
    return ID_pos, name_ID_value


path = "/home/spark159/../../storage/"

# experiment list (cell, sample, agent, tnum)
exp_list = [('H1', 'NCP', 'sp', 8),
            ('H1', 'NCP', 'spd', 6),
            ('H1', 'NCP', 'CoH', 5),
            ('H1', 'NCP', 'PEG', 6),
            ('H1', 'NCP', 'Ca', 5),
            ('H1', 'NCP', 'Mg', 5),
            ('H1', 'NCP', 'HP1a', 3),
            ('H1', 'NCP', 'HP1bSUV', 4),
            ('H1', 'NCP', 'LKH', 3),
            ('H1', 'NCP', 'Ki67', 4),
            ('H1', 'NCP', 'FUS', 5)]

#exp_list = [('H1', 'NCP', 'sp', 8),
#            ('H1', 'NCP', 'HP1a', 3),
#            ('H1', 'NCP', 'LKH', 3),
#            ('H1', 'NCP', 'Ki67', 4)]

exp_list = [('H1', 'NCP', 'sp', 8),
            ('H1', 'NCP', 'Ki67', 4)]


cell = 'H1'


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



# bin size
#bin_size = 10000
bin_size = 1000

# chromosome choice
chr_choices = ['chr1']

# data type
dtype = 'zscore'
#dtype = 'zChalf'

# read data
exp_ID_pos = {}
exp_ID_score = {}
for exp in exp_list:
    cell, sample, agent, tnum = exp

    print "reading %s %s %s" % (cell, sample, agent)

    fname = '_'.join([cell, sample, agent,
                      str(int(bin_size/1000.0)) + 'kb',
                      dtype]) + '.cn'

    if dtype in ['score', 'zscore']:
        field_name = '-'.join([cell, sample, agent, str(tnum)])
        
    elif dtype in ['Chalf', 'zChalf']:
        field_name = 'Chalf'
        
    ID_pos, field_ID_value = read_bin_file(path + fname) 
    ID_score = field_ID_value[field_name]
        
    exp_ID_pos[exp] = ID_pos
    exp_ID_score[exp] = ID_score


# sort by chromosome
exp_chr_scores = {}
for exp in exp_ID_score:
    for ID, score in exp_ID_score[exp].items():
        chr = ID[0]
        if exp not in exp_chr_scores:
            exp_chr_scores[exp] = {}
        if chr not in exp_chr_scores[exp]:
            exp_chr_scores[exp][chr] = []
        exp_chr_scores[exp][chr].append(score)


# plot box plot
chr_scores_list = []
labels = []
for exp in exp_list:
    cell, sample, agent, tnum = exp
    chr_scores_list.append(exp_chr_scores[exp])
    labels.append(agent)

multi_boxplot (chr_scores_list,
               axhline=None,
               ylabel='Condensability (A.U.)',
               title=None,
               keys=chr_list,
               labels=labels,
               colors=[],
               rotation=75,
               legend_loc='best',
               note="")
