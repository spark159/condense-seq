import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
from sklearn.decomposition import NMF
import scipy.sparse as sparse
import pickle
import random
import sklearn.cluster
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import seaborn as sns


def all_path (N, states='ATCG'):
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output

def get_density (ID_CNum, ID_meNum):
    ID_mefrac = {}
    for ID in ID_CNum:
        CNum = ID_CNum[ID]
        if CNum <= 0:
            ID_mefrac[ID] = 0.0
            continue
        meNum = ID_meNum[ID]
        mefrac = float(meNum) / (CNum)
        ID_mefrac[ID] = mefrac
    return ID_mefrac

def Euc_dist (vect1, vect2):
    return np.sqrt(sum([(vect2[k]-vect1[k])**2 for k in range(len(vect1))]))

# Hierarchial clustering
def Hierarchial_clustering (dist_matrix):
    y = squareform(dist_matrix)
    Z = linkage(y, 'ward', optimal_ordering=True)
    #Z = linkage(y, 'average', optimal_ordering=False)

    #idx_cID = [cID-1 for cID in fcluster(Z, t=0.01, criterion='distance')]
    idx_cID = [cID-1 for cID in fcluster(Z, t=3, criterion='maxclust')]
    cID_idxs = {}
    for i in range(len(idx_cID)):
        cID = idx_cID[i]
        if cID not in cID_idxs:
            cID_idxs[cID] = []
        cID_idxs[cID].append(i)
    return Z, idx_cID, cID_idxs


def plot_dendrogram(Z, idx_name=None, node_color=None, name_color=None):
    fig = plt.figure(figsize=(8, 8))
    if node_color != None:
        dendrogram(Z, link_color_func=lambda k: node_color[k], orientation='right')
    else:
        dendrogram(Z, orientation='right')

    if idx_name:
        ax = plt.gca()
        new_labels = []
        old_labels = ax.get_yticklabels()
        for label in old_labels:
            idx = int(label.get_text())
            name = idx_name[idx]
            label.set_text(name)
            if name_color:
                label.set_color(name_color[name])
            new_labels.append(label)
        ax.set_yticklabels(new_labels, weight='bold')

    plt.savefig("dendrogram.svg", format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()


# load annotation data
print "Data reading start"

path = "/home/spark159/../../storage/"

cell = 'H1'

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
#            (cell, 'NCP', 'spd', 6),
#            (cell, 'NCP', 'CoH', 5),
#            (cell, 'NCP', 'PEG', 6),
#            (cell, 'NCP', 'Ca', 5),
#            (cell, 'NCP', 'HP1a', 3),
#            (cell, 'NCP', 'LKH', 3),
#            (cell, 'NCP', 'Ki67', 4)]


exp_list = [(cell, 'NCP', 'sp', 8),
            (cell, 'NCP', 'HP1a', 3),
            (cell, 'NCP', 'LKH', 3),
            (cell, 'NCP', 'Ki67', 4)]


# set chromosome
chr_list = ['chr1']

# set bin size
bin_size= 10000

# set domain names
domain_names = {'NSpeckle':['SON'],
                 'Trnx':['POLR2A', 'POLR2AphosphoS5', 'H3K9ac', 'H3K4me3', 'H3K27ac'],
                 'Polycomb':['CBX8', 'EZH2', 'RNF2', 'SUZ12'],
                 'Hetero':['H3K9me3', 'CBX5'],
                 'Nucleolus':['Nucleolar'],
                 'Lamin':['LaminB1'],
                 'Compartment':['eigen'],
                 'other':['ATcontent']}

domains = ['NSpeckle', 'Trnx', 'Polycomb', 'Hetero', 'Nucleolus', 'Lamin']

names = []
for domain in domains:
    names += domain_names[domain]

# load data
agent_list = []
agent_ID_score = {}
for cell, sample, agent, tnum in exp_list:
    fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', 'anot']) + '.txt'
    score_field = "%s-%s-%s-%d" % (cell, sample, agent, tnum)
    
    field_ID_value = load_file.read_tabular_file(fname, mode='col')
    ID_chr = field_ID_value['Chromosome']

    ID_score = {}
    for ID, chr in ID_chr.items():
        if chr not in chr_list:
            continue
        ID_score[ID] = field_ID_value[score_field][ID]

    agent_list.append(agent)
    agent_ID_score[agent] = ID_score


# select the common IDs
IDs = set([])
for cell, sample, agent, tnum in exp_list:
    if len(IDs) <= 0:
        IDs = set(agent_ID_score[agent].keys())
        continue
    IDs &= set(agent_ID_score[agent].keys())


# get target data
ID_name_value = {}
for ID in IDs:
    if ID not in ID_name_value:
        ID_name_value[ID] = {}
    for agent in agent_ID_score:
        ID_name_value[ID][agent] = agent_ID_score[agent][ID]
    for name in names:
        ID_name_value[ID][name] = field_ID_value[name][ID]

del field_ID_value
del agent_ID_score
del ID_chr
del ID_score


# make data matrix
data = [[] for i in range(len(IDs))]
#for name in names + agent_list:
for name in names:
    values = [ID_name_value[ID][name] for ID in IDs]
    mean = np.mean(values)
    std = np.std(values)
    for i in range(len(IDs)):
        value = values[i]
        re_value = float(value-mean)/std
        #if name == 'eigen':
        #    re_value = - re_value
        data[i].append(re_value)
    del values
    

#X = sparse.csr_matrix(X)

print "Data reading is done"

figsize = (10, 20)
#xticklabels = [agent+'$^{%s}$' % (agent_charge[agent]) for cell, sample, agent, tnum in exp_list]
#xticklabels = names + agent_list
xticklabels = names
cmap = 'bwr'
#cmap = 'Spectral'
#cmap = 'binary'
vmin, vmax = -1.8, 1.8
note = ""

hmap = sns.clustermap(data,
                      #method='ward',
                      #metric='euclidean',
                      metric='correlation',
                      figsize=figsize,
                      cbar_kws=None,
                      row_cluster=True,
                      col_cluster=True,
                      dendrogram_ratio=0.2,
                      colors_ratio=0.03,
                      tree_kws=None,
                      cmap=cmap,
                      center=0,
                      vmin=vmin,
                      vmax=vmax,
                      xticklabels=xticklabels,
                      cbar_pos=None)

plt.gca().set_yticklabels([])
plt.gca().tick_params(right='off')
plt.gca().xaxis.tick_top()
plt.xticks(rotation=90, ha='left', va='center', rotation_mode='anchor', weight='bold')

plt.savefig("hmap_" + note + ".png", dpi=500, bbox_inches='tight')
#plt.show()
plt.close()

sys.exit(1)

# hierarchial clustering
dist_matrix = np.zeros((len(IDs), len(IDs)))

for i in range(len(IDs)-1):
    for j in range(i+1, len(IDs)):
        dist = Euc_dist(X[i], X[j])
        dist_matrix[i][j] = dist
        dist_matrix[j][i] = dist

Z, idx_cID, cID_idxs = Hierarchial_clustering (dist_matrix)
plot_dendrogram(Z)
