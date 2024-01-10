import load_file
import graphics
import statis
import sys
import copy
import pickle
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
import Interval_dict
import seaborn as sns

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

exp_list = [(cell, 'NCP', 'sp', 8)]



# binsize of input data
#bin_size = 10000
#bin_size = 5000
bin_size = 1000

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

chr_list = ['chr1']

# set data names to check
target_names = []
for cell, sample, agent, tnum in exp_list:
    target_names.append('-'.join([cell, sample, agent, str(tnum)]))

# set domain names
feature_names = {'genetic':['ATcontent'],
                 'methylation':['meCpG density', 'meCHG density', 'meCHH density'],
                 'PTMs':['H2AFZ',
                         'H2AK5ac',
                         'H2BK120ac',
                         'H2BK12ac',
                         'H2BK15ac',
                         'H2BK20ac',
                         'H2BK5ac',
                         'H3K14ac',
                         'H3K18ac',
                         'H3K23ac',
                         'H3K23me2',
                         'H3K27ac',
                         'H3K27me3',
                         'H3K36me3',
                         'H3K4ac',
                         'H3K4me1',
                         'H3K4me2',
                         'H3K4me3',
                         'H3K56ac',
                         'H3K79me1',
                         'H3K79me2',
                         'H3K9ac',
                         'H3K9me3',
                         'H4K20me1',
                         'H4K5ac',
                         'H4K8ac',
                         'H4K91ac'],
                 'subPTMs':['H2AFZ',
                            'H3K14ac',
                            'H3K18ac',
                            'H3K23ac',
                            'H3K23me2',
                            'H3K36me3',
                            'H3K4ac',
                            'H3K4me1',
                            'H3K4me2',
                            'H3K56ac',
                            'H3K79me1',
                            'H3K79me2',
                            'H4K20me1'],
                 'NSpeckle':['SON'],
                 'Trnx':['POLR2A', 'POLR2AphosphoS5', 'H3K9ac', 'H3K4me3', 'H3K27ac'],
                 'Polycomb':['CBX8', 'EZH2', 'RNF2', 'SUZ12'],
                 'Hetero':['H3K9me3', 'CBX5'],
                 'Nucleolus':['Nucleolar'],
                 'Lamin':['LaminB1'],
                 'CTCF-Cohesin': ['CTCF', 'Rad21'],
                 'Compartment':['eigen']}

names = []
for feature in feature_names:
    names += feature_names[feature]


# other parameters
dtype = 'zscore'
note = '%skb' % (int(bin_size/1000.0))


# HMM parameters
state_num = 10
HMM_type = 'Gaussian'


# load HMM data
HMM_fname = '_'.join([str(state_num), 'ConHMMstate', HMM_type, str(int(bin_size/1000.0))+'kb'])
chr_state_intervals = read_chromHMM(HMM_fname + '.bed', chr_choice=chr_list)
HMM_transm = pickle.load(open(HMM_fname + "_transm.pickle", "rb"))
HMM_means = pickle.load(open(HMM_fname + "_means.pickle", "rb"))
HMM_covars = pickle.load(open(HMM_fname + "_covars.pickle", "rb"))


# set HMM new name
if state_num ==10:
    HMM_name_dict = {"E4":"GC-rich intermediate",
                     "E1":"Strong condensate",
                     "E6":"HP1$\\alpha$ condensate",
                     "E9":"AT-rich soluble",
                     "E8":"GC-rich soluble",
                     "E5":"Polyamine soluble",
                     "E7":"Strong soluble",
                     "E2":"HP1$\\alpha$ soluble",
                     "E3":"AT-rich intermediate",
                     "E10":"Ki67 condensate"}

    #HMM_states = ["E4", "E1", "E6", "E9", "E8", "E5", "E7", "E2", "E3", "E10"]
    HMM_states = ["E1", "E6", "E10", "E3", "E4", "E2", "E9", "E8", "E5", "E7"]


# plot HMM emission matrix
img = []
for state in HMM_states:
    idx = int(state[1:]) - 1
    img.append(HMM_means[idx])

xticklabels = ['sp', 'HP1$\\alpha$', 'LKH', 'Ki67']
yticklabels = HMM_states
cmap = 'bwr_r'
vmin, vmax = -2, 2

row_num, col_num = HMM_means.shape
height = 0.4*row_num
width = 0.4*col_num

fig = plt.figure(figsize=(width, height))
plt.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
plt.xticks(range(len(xticklabels)), xticklabels,
           ha='right', va='center', rotation_mode='anchor', rotation=45)
plt.yticks(range(len(yticklabels)), [])
plt.savefig(HMM_fname + "_emission.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()


# load annotation data
name_binID_value = {}
binID_name_value = {}
binID_start, binID_end = {}, {}
for cell, sample, agent, tnum in exp_list:
    print "reading %s %s %s" % (cell, sample, agent)

    # check methylation density need to be computed
    search_names = []
    density_names = []
    for name in names:
        if name.endswith('density'):
            motif, _ = name.split()
            motif = motif[2:]
            search_names.append('CNumber' + '(' + motif + ')')
            search_names.append('meCNumber' + '(' + motif + ')')
            density_names.append(name)

        else:
            search_names.append(name)

    # load file
    fname = '_'.join([cell, sample, agent,
                      str(int(bin_size/1000.0)) + 'kb',
                      dtype, 'anot']) + '.txt'

    ID_chr, ID_start, ID_end, name_ID_value = load_file.read_anot_file (path + fname,
                                                                        chr_choice=chr_list,
                                                                        target_names=search_names)

    # compute density data
    for name in density_names:
        motif, _ = name.split()
        motif = motif[2:]
        ID_control = name_ID_value['CNumber' + '(' + motif + ')']
        ID_test = name_ID_value['meCNumber' + '(' + motif + ')']

        ID_density = {}
        for ID in ID_control:
            control = ID_control[ID]
            test = ID_test[ID]
            if control <=0:
                density = np.nan
            else:
                density = float(test)/control
            ID_density[ID] = density
                
        del name_ID_value['CNumber' + '(' + motif + ')']
        del name_ID_value['meCNumber' + '(' + motif + ')']
        name_ID_value[name] = ID_density


    # change ID to binID
    for ID in ID_chr:
        chr = ID_chr[ID]
        start = ID_start[ID]
        end = ID_end[ID]
        binID = (chr, start, end)
        
        for name in names:            
            if name not in name_binID_value:
                name_binID_value[name] = {}
            name_binID_value[name][binID] = name_ID_value[name][ID]

        if binID not in binID_name_value:
            binID_name_value[binID] = {}
        for name in names:
            binID_name_value[binID][name] = name_ID_value[name][ID]

        binID_start[binID] = int(start)
        binID_end[binID] = int(end)

    del name_ID_value


# find common IDs for names
common_binIDs = set([])
for name in name_binID_value:
    binIDs = set(name_binID_value[name].keys())
    if len(common_binIDs) == 0:
        common_binIDs |= binIDs
        continue
    common_binIDs &= binIDs
common_binIDs = list(common_binIDs)

chr_binIDs = {}
for binID in common_binIDs:
    chr, start, end = binID
    if chr not in chr_binIDs:
        chr_binIDs[chr] = []
    chr_binIDs[chr].append(binID)

    
# sort each bins by HMM state
chr_state_binIDs = {}
for chr_name in chr_list:

    state_intervals = chr_state_intervals[chr_name]

    sID_interval = {}
    for state in state_intervals:
        intervals = state_intervals[state]
        for k in range(len(intervals)):
            sID = (state, k)
            sID_interval[sID] = intervals[k]

    state_dict = Interval_dict.double_hash(sID_interval,
                                           domain_size=10000,
                                           max_pos=5*10**8)

    for binID in chr_binIDs[chr_name]:
        bst = binID_start[binID]
        bed = binID_end[binID]
        find_sIDs = state_dict.insert_range(bst, bed, 1)
        value_sID = []
        for sID in find_sIDs:
            value = state_dict.ID_value[sID]
            value_sID.append((value, sID))
            del state_dict.ID_value[sID]

        sID = sorted(value_sID, reverse=True)[0][1]
        state = sID[0]
        
        if chr not in chr_state_binIDs:
            chr_state_binIDs[chr] = {}
        if state not in chr_state_binIDs[chr]:
            chr_state_binIDs[chr][state] = []
        chr_state_binIDs[chr][state].append(binID)
    

# get values for corresponding HMM state
name_state_values = {}
state_name_values = {}
for chr in chr_state_binIDs:
    for state in chr_state_binIDs[chr]:
        binIDs = chr_state_binIDs[chr][state]
        for name in names:
            values = [binID_name_value[binID][name] for binID in binIDs]
            if name not in name_state_values:
                name_state_values[name] = {}
            name_state_values[name][state] =  values
            if state not in state_name_values:
                state_name_values[state] = {}
            state_name_values[state][name] = values

            
# get mean value for corresponding HMM state
name_state_mvalue = {}
state_name_mvalue = {}
for name in names:
    for state in name_state_values[name]:
        mvalue = np.nanmean(name_state_values[name][state])
        if name not in name_state_mvalue:
            name_state_mvalue[name] = {}
        name_state_mvalue[name][state] = mvalue
        if state not in state_name_mvalue:
            state_name_mvalue[state] = {}
        state_name_mvalue[state][name] = mvalue
        

# plot boxplot
for name in names:
    state_values = name_state_values[name]
    #single_boxplot (state_values,
    #                ylabel='A.U.',
    #                title=name,
    #                keys=HMM_states,
    #                rotation=75,
    #                ylim=[None, None],
    #                note=note + name)


# plot heatmap
features = ['genetic',
            'methylation',
            'PTMs',
            'CTCF-Cohesin',
            'NSpeckle',
            'Trnx',
            'Polycomb',
            'Hetero',
            'Nucleolus',
            'Lamin',
            'Compartment']

name_cmap = {#'ATcontent':'Spectral_r',
             'SON':'Spectral_r', 
             'Nucleolar':'Spectral_r',
             'LaminB1':'Spectral_r',
             'eigen':'Spectral_r'}


#cmap_list = ['Reds', 'Blues', 'Greens', 'Purples', 'Oranges',  'Greys']*10
cmap_list = ['Greys']*50

for feature, cmap in zip(features, cmap_list):
    subnames = feature_names[feature]
    img_list = []
    for j in range(len(subnames)):
        img = np.zeros((len(HMM_states), len(subnames)))
        img[:] = np.nan
        for i in range(len(HMM_states)):
            state = HMM_states[i]
            name = subnames[j]
            mvalue = name_state_mvalue[name][state]
            img[i][j] = mvalue
        img_list.append(img)
    width = 0.4*len(subnames)
    height = 0.4*state_num
    fig = plt.figure(figsize=(width, height))
    for img, name in zip(img_list, subnames):
        try:
            cmap = name_cmap[name]
        except:
            cmap = 'Greys'
        if name == 'ATcontent':
            center = 0.5
            vmin = center-max([abs(np.nanmin(img)-center), abs(np.nanmax(img)-center)])
            vmax = center+max([abs(np.nanmin(img)-center), abs(np.nanmax(img)-center)])
        elif name in ['eigen', 'Nucleolar', 'LaminB1']:
            center = 0
            vmin = -max([abs(np.nanmin(img)-center), abs(np.nanmax(img)-center)])
            vmax = max([abs(np.nanmin(img)-center), abs(np.nanmax(img)-center)])
        else:
            vmin, vmax = None, None
        plt.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
    plt.xticks(range(len(subnames)), subnames, ha='right', va='center', rotation_mode='anchor', rotation=45)
    #plt.yticks(range(len(HMM_states)), HMM_states)
    plt.yticks(range(len(HMM_states)), [])
    #plt.ylabel("State #")
    #plt.gca().tick_params(left='off')
    plt.savefig(HMM_fname + "_%s_heat.svg" % (feature), format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

# check each state size
state_size = {}
for chr_name in chr_list:
    state_intervals = chr_state_intervals[chr_name]
    for state in state_intervals:
        for interval in state_intervals[state]:
            st, ed = interval
            if state not in state_size:
                state_size[state] = 0
            state_size[state] += ed - st

for state, size in state_size.items():
    size = float(size) / 10000000
    radius = np.sqrt(float(size)/(4*np.pi))
    fig = plt.figure(figsize=(4,4))
    circle = plt.Circle((0.5, 0.5), radius=radius, Fill=False, clip_on=False)
    plt.gca().add_patch(circle)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().tick_params(top='off',bottom='off',left='off',right='off', labelleft='off',labelbottom='off')
    plt.xlim([-0.5, 1.5])
    plt.ylim([-0.5, 1.5])
    #plt.savefig(HMM_fname + "_%s_circle.svg" % (state), format='svg', bbox_inches='tight')
    plt.close()
    


sys.exit(1)
# plot HMM mean matrix
xticklabels = ['sp', 'HP1$\\alpha$', 'LKH', 'Ki67']
yticklabels = HMM_states
cmap = 'bwr_r'
vmin, vmax = -2, 2
note = ""

figsize = (4, 10)
hmap = sns.clustermap(HMM_means,
                      #method='ward',
                      #metric='euclidean',
                      metric='cosine',
                      figsize=figsize,
                      cbar_kws=None,
                      row_cluster=True,
                      col_cluster=False,
                      dendrogram_ratio=0.3,
                      #colors_ratio=0.03,
                      tree_kws={'linewidths':1.8},
                      cmap=cmap,
                      center=0,
                      vmin=vmin,
                      vmax=vmax,
                      xticklabels=xticklabels,
                      yticklabels=yticklabels,
                      linewidth=0,
                      rasterized=True,
                      edgecolor=None,
                      cbar_pos=None)

plt.gca().tick_params(left='on', right='off', bottom='on')
plt.gca().xaxis.tick_top()
for _, spine in plt.gca().spines.items():
    spine.set(visible=True, lw=1.5, edgecolor="black")
plt.gca().xaxis.tick_bottom() # x axis on top
plt.gca().xaxis.set_label_position('bottom')
plt.xticks(rotation=45, fontsize=18, ha='right', va='center', rotation_mode='anchor')
plt.yticks(rotation=0, ha='left', va='center', rotation_mode='anchor')
plt.savefig(HMM_fname + "_means.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()

# got new order of HMM states
HMM_states = []
for tick_label in hmap.ax_heatmap.axes.get_yticklabels():
    HMM_states.append(tick_label.get_text())    
