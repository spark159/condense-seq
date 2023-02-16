import load_file_python3
import Interval_dict_python3
import sys
import copy
import pickle
import matplotlib.pyplot as plt
import numpy as np
import math
import networkx as nx
import pandas as pd
import upsetplot
from upsetplot import UpSet

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

#exp_list = [(cell, 'NCP', 'sp', 8)]

exp_list = [(cell, 'NCP', 'sp', 8),
            (cell, 'NCP', 'HP1a', 3)]

#exp_list = [(cell, 'NCP', 'sp', 8),
#            (cell, 'NCP', 'spd', 6),
#            (cell, 'NCP', 'CoH', 5),
#            (cell, 'NCP', 'PEG', 6),
#            (cell, 'NCP', 'Ca', 5),
#            (cell, 'NCP', 'HP1a', 3),
#            (cell, 'NCP', 'HP1bSUV', 4)]






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
features = []
for feature in features:
    names += feature_names[feature]

for cell, sample, agent, tnum in exp_list:
    names.append('-'.join([cell, sample, agent, str(tnum)]))


# other parameters
dtype = 'zscore'
note = '%skb' % (int(bin_size/1000.0))


# nuclear domain files
dname_fname = {"SPAD (NSpeckle)":"SPAD_selected.bed",
               "NAD (Nucleolus)":"NAD_selected.bed",
               "LAD (Lamin)":"LAD_selected.bed",
               "B-compartment":"eigen_H1_100kb_state.bed",
               "chrom":"H1_12_segments.bed"}

# change HMM state name
dname_change = {"SPAD (NSpeckle)":{'E2':1},
                "NAD (Nucleolus)":{'E1':1},
                "LAD (Lamin)":{'E2':1},
                "B-compartment":{'E2':1}}

#chrom_change = {"H3K9me3": {"E9":1},
#                "H3K27me3":{"E1":1},
#                "Promoter":{"E3":1,"E5":1},
#                "Enhancer":{"E4":1, "E6":1},
#                "Gene body":{"E10":1, "E11":1},
#                "Poised promoter":{"E2":1}}
#chroms = ["Promoter", "Enhancer", "Gene body", "Poised promoter", "H3K27me3", "H3K9me3"]

chrom_change = {"H3K9me3 enriched": {"E9":1},
                "Promoter/Enhancer":{"E2":1, "E3":1,"E4":1, "E5":1, "E6":1},
                "Tx elongation":{"E10":1, "E11":1}}
chroms = ["Promoter/Enhancer", "H3K9me3 enriched", "Tx elongation"]





# select domain names to be analyzed
#domains = ['chrom', 'SPAD', 'NAD', 'LAD', 'eigen']
#domains = ['eigen', 'SPAD', 'NAD', 'LAD', "chrom"]
#domains = ['SPAD', 'LAD', 'NAD', 'B-compt']
#domains = ['chrom', 'SPAD (NSpeckle)', 'LAD (Lamin)', 'NAD (Nucleolus)', 'B-compartment']
domains = ['chrom', 'SPAD (NSpeckle)', 'LAD (Lamin)', 'NAD (Nucleolus)']
 


# load domain files
dname_chr_state_intervals = {}
for dname in domains:
    chr_state_intervals = read_chromHMM(dname_fname[dname], chr_choice=chr_list)
    dname_chr_state_intervals[dname] = chr_state_intervals


# load annotation data
name_binID_value = {}
binID_name_value = {}
for cell, sample, agent, tnum in exp_list:
    print ("reading %s %s %s" % (cell, sample, agent))

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

    ID_chr, ID_start, ID_end, name_ID_value = load_file_python3.read_anot_file (path + fname,
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


    # change ID to binID (chr, st, ed)
    for ID in ID_chr:
        chr = ID_chr[ID]
        start = ID_start[ID]
        end = ID_end[ID]
        binID = (chr, start, end)
        
        for name in name_ID_value:            
            if name not in name_binID_value:
                name_binID_value[name] = {}
            name_binID_value[name][binID] = name_ID_value[name][ID]

        if binID not in binID_name_value:
            binID_name_value[binID] = {}
        for name in name_ID_value:
            binID_name_value[binID][name] = name_ID_value[name][ID]

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
states_binIDs = {}
binID_states = {}
for chr_name in chr_list:
    dname_dict = {}
    for dname in domains:        
        chr_state_intervals = dname_chr_state_intervals[dname]
        state_intervals = chr_state_intervals[chr_name]

        sID_interval = {}
        for state in state_intervals:
            intervals = state_intervals[state]
            for k in range(len(intervals)):
                sID = (state, k)
                sID_interval[sID] = intervals[k]

        state_dict = Interval_dict_python3.double_hash(sID_interval,
                                                       domain_size=10000,
                                                       max_pos=5*10**8)
        dname_dict[dname] = state_dict

    for binID in chr_binIDs[chr_name]:
        _, bst, bed = binID
        
        states = []
        for dname in domains:
            state_dict = dname_dict[dname]
            find_sIDs = state_dict.insert_range(bst, bed, 1)

            if len(find_sIDs) > 0:
                value_sID = []
                for sID in find_sIDs:
                    value = state_dict.ID_value[sID]
                    value_sID.append((value, sID))
                    del state_dict.ID_value[sID]

                sID = sorted(value_sID, reverse=True)[0][1]
                state = sID[0]
            else:
                state = 'E0'

            if dname == 'chrom':
                for chrom in chroms:
                    try:
                        new_state = chrom_change[chrom][state]
                    except:
                        new_state = 0
                    states.append(new_state)
            else:
                try:
                    new_state = dname_change[dname][state]
                except:
                    new_state = 0
                states.append(new_state)

        states = tuple(states)

        if states not in states_binIDs:
            states_binIDs[states] = []
        states_binIDs[states].append(binID)

        assert binID not in binID_states
        binID_states[binID] = states


# get values for corresponding HMM state
name_states_values = {}
for states in states_binIDs:
    binIDs = states_binIDs[states]
    for name in names:
        values = [name_binID_value[name][binID] for binID in binIDs]
        if name not in name_states_values:
            name_states_values[name] = {}
        name_states_values[name][states] =  values

            
# get mean value for corresponding HMM state
name_states_mvalue = {}
for name in names:
    for states in name_states_values[name]:
        mvalue = np.nanmean(name_states_values[name][states])
        if name not in name_states_mvalue:
            name_states_mvalue[name] = {}
        name_states_mvalue[name][states] = mvalue

### UpSet plot
size_states = [(len(binIDs), states) for states, binIDs in states_binIDs.items()]
size_states = sorted(size_states, reverse=True)
cutoff = int(size_states[0][0]*0.01)
#cutoff = 0
states_list = []
for size, states in size_states:
    if size > cutoff:
        states_list.append(states)

sizes = []
index = []
for states in states_list:
    row = []
    for state in states:
        if state == 0:
            row.append(False)
        elif state != 0:
            row.append(True)
    index.append(tuple(row))
    sizes.append(len(states_binIDs[states]))

labels = []
for dname in domains:
    if dname == 'chrom':
        for chrom in chroms:
            labels.append(chrom)
    else:
        labels.append(dname)
            
index = pd.MultiIndex.from_tuples(index, names=labels)
sizes = pd.Series(sizes, index=index)
upset = UpSet(sizes, sort_by='input', sort_categories_by='-input')

# size (upset plot)
fig = plt.figure()
upset.plot()
plt.savefig("UpSet_size.svg", format="svg", bbox_inches='tight')
#plt.show()
plt.close()

# scores (box plot)
nrows = len(names)+2
ncols = 2
colors = ['tab:blue','tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan'] * 10
fig, axes = plt.subplots(figsize=(8,10), nrows=nrows, ncols=ncols,
                         gridspec_kw={'height_ratios':[2]*len(names)+[2]+[3],
                                      'width_ratios':[0.5, 5]})
fig.subplots_adjust(wspace=0.65)
for i in range(nrows):
    for j in range(ncols):
        if i < len(names):
            if j == 0:
                axes[i][j].spines['top'].set_visible(False)
                axes[i][j].spines['bottom'].set_visible(False)
                axes[i][j].spines['left'].set_visible(False)
                axes[i][j].spines['right'].set_visible(False)
                axes[i][j].set_xticks([])
                axes[i][j].set_xticklabels([])
                axes[i][j].set_yticks([])
                axes[i][j].set_yticklabels([])
            else:
                assert j == 1
                name = names[i]
                states_values = name_states_values[name]
                box_data = [states_values[states] for states in states_list]
                bplot = axes[i][j].boxplot(box_data, 0, '', patch_artist=True)
                for patch, color in zip(bplot['boxes'], colors):
                    patch.set_facecolor(color)
                axes[i][j].set_ylim([-5, 5])
                axes[i][j].spines['top'].set_visible(False)
                axes[i][j].spines['bottom'].set_visible(False)
                axes[i][j].spines['left'].set_visible(True)
                axes[i][j].spines['right'].set_visible(False)
                axes[i][j].axhline(y=0, color='k', linestyle='--', alpha=0.3, zorder=1)
                axes[i][j].set_xticks([])
                axes[i][j].set_xticklabels([])
                cell, sample, agent, tnum = name.split('-')
                axes[i][j].set_ylabel(agent_fullname[agent], va='center', ha='right', rotation=0)
        else:
            if i == len(names):
                if j == 0:
                    axes[i][j].spines['top'].set_visible(False)
                    axes[i][j].spines['bottom'].set_visible(False)
                    axes[i][j].spines['left'].set_visible(False)
                    axes[i][j].spines['right'].set_visible(False)
                    axes[i][j].set_xticks([])
                    axes[i][j].set_xticklabels([])
                    axes[i][j].set_yticks([])
                    axes[i][j].set_yticklabels([])
                else:
                    assert j ==1
                    upset.plot_intersections(axes[i][j])
                    axes[i][j].set_xlim([-0.5, len(states_list)-0.5])
                    #axes[i][j].set_ylabel('intersection size', va='center', ha='right', rotation=0)
            else:
                assert i == len(names) +1
                if j == 0:
                    pass
                    upset.plot_totals(axes[i][j])
                    axes[i][j].set_ylim([-0.5, len(labels)-0.5])
                else:
                    assert j ==1
                    upset.plot_matrix(axes[i][j])
                    axes[i][j].set_xlim([-0.5, len(states_list)-0.5])
                    axes[i][j].set_ylim([-0.5, len(labels)-0.5])
plt.savefig("box_upset.svg", format="svg", bbox_inches='tight')
#plt.show()
plt.close()

# heat map for coloring each intersection mvalue
img = []
for states in states_list:
    row = []
    for name in names:
        mvalue = name_states_mvalue[name][states]
        row.append(mvalue)
    img.append(row)
    
fig = plt.figure()
plt.imshow(img, cmap='bwr_r', vmin=-2.5, vmax=2.5)
plt.xlabel([name.split('-')[2] for name in names])
plt.ylabel(states)
plt.colorbar(shrink=0.25)
plt.savefig("venn_color.svg", format="svg", bbox_inches='tight')
#plt.show()
plt.close()






sys.exit(1)



















# Network plot
#def check_edge (states1, states2):
#    check = False
#    for i in range(len(states1)):
#        if states1[i] == states2[i] and states1[i] != 0:
#            check = True
#            break
#    return check

    
#nodes = sorted(states_binIDs.keys())
#node_sizes = [len(states_binIDs[node]) for node in nodes]
#edges = []
#for i in range(len(nodes)-1):
#    for j in range(i+1, len(nodes)):
#        node1 = nodes[i]
#        node2 = nodes[j]
#        if check_edge(node1, node2):
#            edges.append((node1, node2))
        

#G = nx.Graph()
#G.add_nodes_from(nodes)
#G.add_edges_from(edges)
#nx.draw_networkx(G, pos=nx.spring_layout(G), node_size = node_sizes)
#plt.axis('off')
#plt.show()
#plt.close()



# PCA plot
#data = []
#for binID in sorted(binID_states):
#    data.append(list(binID_states[binID]))

#from sklearn.manifold import Isomap
#embedding = Isomap(n_components=2)
#trans_data = embedding.fit_transform(data).T
#fig = plt.figure()
#plt.plot(trans_data[0], trans_data[1], '.')
#plt.title("Isomap")
#plt.savefig("Isomap.svg", format='svg', bbox_inches='tight')
#plt.show()
#plt.close()

#pca = PCA().fit(data)
#trans_data = pca.transform(data).T

#fig = plt.figure()
#plt.plot(trans_data[0], trans_data[1], '.')
#plt.title("PCA plot")
#plt.savefig("PCA.svg", format='svg', bbox_inches='tight')
#plt.show()
#plt.close()

# tSNE plot
#data = []
#for binID in sorted(binID_states):
#    data.append(list(binID_states[binID]))

#perp=10
#tsne = sklearn.manifold.TSNE(n_components=2, perplexity=perp, init='pca', random_state=0)
#trans_data = tsne.fit_transform(data).T
#fig = plt.figure()
#plt.plot(trans_data[0], trans_data[1], '.')
#plt.title("tSNE plot")
#plt.savefig("tSNE.svg", format='svg', bbox_inches='tight')
#plt.show()
#plt.close()

#size_states = [(len(binIDs), states) for states, binIDs in states_binIDs.items()]
#size_states = sorted(size_states, reverse=True)
#states_list = [states for size, states in size_states]

#sizes = []
#index = []
#for states in states_list:
#    row = []
#    for state in states:
#        if state == 0:
#            row.append(False)
#        elif state != 0:
#            row.append(True)
#    index.append(tuple(row))
#    sizes.append(len(states_binIDs[states]))

#labels = []
#for dname in domains:
#    if dname == 'chrom':
#        for chrom in chroms:
#            labels.append(chrom)
#    else:
#        labels.append(dname)
            
#index = pd.MultiIndex.from_tuples(index, names=labels)
#cutoff = int(max(sizes)*0.1)
#sizes = pd.Series(sizes, index=index)
#upset = UpSet(sizes, sort_by='input', sort_categories_by='-input', min_subset_size=cutoff)

#fig = plt.figure()
#upset.plot()
#plt.savefig("UpSet_size.svg", format="svg", bbox_inches='tight')
#plt.show()
#plt.close()

#fig, axes = plt.subplots(nrows=1+len(names), ncols=1, sharex=True)
#for i in range(len(names)):
#    states_values = name_states_values[name]
#    box_data = [states_values[states] for states in states_list]
#    axes[i].boxplot(box_data, 0, '')
#upset.plot()
#upset.plot()
#plt.savefig("UpSet_all.svg", format="svg", bbox_inches='tight')
#plt.show()
#plt.close()


#index = []
#data = []
#for binID in sorted(binID_states):
#    states = binID_states[binID]
#    row = []
#    for state in states:
#        if state == 0:
#            row.append(False)
#        elif state != 0:
#            row.append(True)
#    index.append(tuple(row))
#    row = []
#    for name in names:
#        value = name_binID_value[name][binID]
#        row.append(value)
#    data.append(row)

#labels = []
#for dname in domains:
#    if dname == 'chrom':
#        for chrom in chroms:
#            labels.append(chrom)
#    else:
#        labels.append(dname)

#index = pd.MultiIndex.from_tuples(index, names=labels)
#data = pd.DataFrame(data, index=index, columns=names)

#upset = UpSet(data, subset_size='count', sort_by='cardinality', sort_categories_by='-input')
#for name in names[::-1]:
#    upset.add_catplot(value=name, kind='box', showfliers=False)

#fig = plt.figure(figsize=(10,5))
#upset.plot()
#plt.savefig("UpSet.svg", format="svg", bbox_inches='tight')
##plt.show()
#plt.close()

#for name in names:
#    states_mvalue = name_states_mvalue[name]
#    mvalues = []
#    for states in sorted(states_mvalue):
#        mvalues.append(states_mvalue[states])
#    fig = plt.figure()
#    data = pd.Series(mvalues, index=index)
#    upsetplot.plot(data)
#    plt.savefig("UpSet_%s.svg" % (name), format="svg", bbox_inches='tight')
#    #plt.show()
#    plt.close()
