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
import random
from scipy.optimize import curve_fit
from sklearn import linear_model
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib
import pickle
from scipy.stats import norm
import glob
import statis
import sklearn.cluster
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform


# rescale the data in old range (old_st, old_ed) into new range (new_st, new_ed)
def rescale (value, old_st, old_ed, new_st, new_ed):
    if value < old_st:
        return new_st
    elif value > old_ed:
        return new_ed
    else:
        return new_st + (new_ed - new_st)*float(value-old_st)/(old_ed-old_st)
    

def read_rank (fname):
    gene_value = {}
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        gene, value = cols
        value = float(value)
        gene_value[gene] = value
    return gene_value
gene_value = read_rank('KO-WT.rnk')

def read_GSEA (path):
    def read_report (fname, cutoff=100):
        gs_list = []
        First = True
        for line in open(fname):
            if First:
                First = False
                continue
            cols = line.strip().split('\t')
            name, size, nes = cols[0], cols[3], cols[5]
            size = int(size)
            nes = float(nes)
            gs_list.append({'name':name, 'size':size, 'nes':nes})
            if len(gs_list) >= cutoff:
                break
        return gs_list

    def read_gs (fname):
        gene_info = {}
        First = True
        for line in open(fname):
            if First:
                First = False
                continue
            cols = line.strip().split('\t')
            gene, rank, core, es = cols[1], cols[2], cols[-1], cols[-2]
            rank = int(rank)
            es = float(es)
            if core == 'Yes':
                core = True
            else:
                core = False
            if gene not in gene_info:
                gene_info[gene] = {}
            gene_info[gene]['rank'] = rank
            gene_info[gene]['core'] = core
            gene_info[gene]['es'] = es
        return gene_info

    pos_fname = glob.glob(path + "/gsea_report_for_na_pos_*.tsv")[0]
    pos_gs_list = []
    for gs in read_report(pos_fname):
        try:
            gs_name = gs['name']
            gene_info = read_gs(path + '/' + gs_name + '.tsv')
            gs['genes'] = gene_info
            pos_gs_list.append(gs)
        except:
            break

    neg_fname = glob.glob(path + "/gsea_report_for_na_neg_*.tsv")[0]
    neg_gs_list = []
    for gs in read_report(neg_fname):
        try:
            gs_name = gs['name']
            gene_info = read_gs(path + '/' + gs_name + '.tsv')
            gs['genes'] = gene_info
            neg_gs_list.append(gs)
        except:
            break

    return pos_gs_list, neg_gs_list


def select_best (gs_list, num):
    nes_gs = []
    for gs in gs_list:
        nes = gs['nes']
        nes_gs.append((abs(nes), gs))
    nes_gs = sorted(nes_gs, reverse=True)
    return [gs for _, gs in nes_gs[:num]]

def get_kappa (bvec1, bvec2):
    assert len(bvec1) == len(bvec2)
    bpair_count = {(0,0):0, (1,0):0, (0,1):0, (1,1):0 }
    for i in range(len(bvec1)):
        bpair = (bvec1[i], bvec2[i])
        bpair_count[bpair] +=1

    sum0_ = bpair_count[(0,0)] + bpair_count[(0,1)]
    sum1_ = bpair_count[(1,0)] + bpair_count[(1,1)]
    sum_0 = bpair_count[(0,0)] + bpair_count[(1,0)]
    sum_1 = bpair_count[(0,1)] + bpair_count[(1,1)]
    total = sum(bpair_count.values())
    obs = float(bpair_count[(1,1)]+bpair_count[(0,0)])/total
    exp = float(sum_1*sum1_+sum_0*sum0_)/(total*total)
    kappa = float(obs-exp)/(1-exp)
    assert kappa <= 1
    return kappa


def get_kappa_matrix (gs_list, gene_list=None):
    # make union gene_list to be considered
    if gene_list == None:
        gene_list = set([])
        for gs in gs_list:
            gene_list |= set(gs['genes'].keys())
        gene_list = list(gene_list)

    # make binary vector for each gene-set
    idx_bvec = []
    for i in range(len(gs_list)):
        gs = gs_list[i]
        bvec = []
        for gene in gene_list:
            try:
                gs['genes'][gene]
                bvec.append(1)
            except:
                bvec.append(0)
        idx_bvec.append(bvec)

    # make kappa matrxi
    kappa_matrix = np.zeros((len(gs_list), len(gs_list)))
    kappa_matrix[:] = np.nan
    for i in range(len(gs_list)):
        for j in range(i, len(gs_list)):
            bvec1, bvec2 = idx_bvec[i], idx_bvec[j]
            kappa = get_kappa(bvec1, bvec2)
            kappa_matrix[i][j] = kappa
            kappa_matrix[j][i] = kappa

    return gene_list, kappa_matrix

# convert kappa matrix to score matrix
def kappa_to_score (kappa_matrix, scale=1.0):
    nrow, ncol = kappa_matrix.shape
    score_matrix = np.zeros((nrow, ncol))
    for i in range(nrow):
        for j in range(ncol):
            kappa = kappa_matrix[i][j]
            score_matrix[i][j] = np.exp(-scale*(1.0-kappa))
    return score_matrix

# convert kappa matrix to distance matrix
def kappa_to_dist (kappa_matrix, scale=1.0):
    nrow, ncol = kappa_matrix.shape
    dist_matrix = np.zeros((nrow, ncol))
    for i in range(nrow):
        for j in range(ncol):
            if i == j:
                dist_matrix[i][j] = 0
                continue
            kappa = kappa_matrix[i][j]
            dist = scale*(1.0 - kappa)
            #dist_matrix[i][j] = np.exp(dist)
            dist_matrix[i][j] = dist
    return dist_matrix

# Spectral clustering
def Spectral_clustering (score_matrix, cluster_num):
    idx_cID = sklearn.cluster.spectral_clustering(affinity=score_matrix,
                                                  n_clusters=cluster_num,
                                                  random_state=0)
    cID_idxs = {}
    for i in range(len(idx_cID)):
        cID = idx_cID[i]
        if cID not in cID_idxs:
            cID_idxs[cID] = []
        cID_idxs[cID].append(i)
    return idx_cID, cID_idxs

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

# OPTICS clustering
def OPTICS_clustering (dist_matrix):
    idx_cID = sklearn.cluster.OPTICS(min_samples=2, metric='precomputed').fit_predict(dist_matrix)

    cID_idxs = {}
    for i in range(len(idx_cID)):
        cID = idx_cID[i]
        if cID not in cID_idxs:
            cID_idxs[cID] = []
        cID_idxs[cID].append(i)
    return idx_cID, cID_idxs

# decode Z (linkage) information
def decode_Z (Z, idx_key):
    node_children = {i:{} for i in range(len(idx_key))}
    node_dist = {i:None for i in range(len(idx_key))}
    node_keys = {i:{idx_key[i]} for i in range(len(idx_key))}
    for i in range(len(Z)):
        node1, node2 = int(Z[i][0]), int(Z[i][1])
        new_node = max(node_keys.keys()) + 1
        node_children[new_node] = set([node1, node2])
        node_dist[new_node] = float(Z[i][2])
        node_keys[new_node] = node_keys[node1] | node_keys[node2]
    return node_children, node_dist, node_keys

# encode Z (linkage)
def encode_Z (node_children, node_dist, node_keys):
    idx_node = []
    for node in node_children:
        idx_node.append(node)
    idx_node = sorted(idx_node)

    node_idx = {}
    for i in range(len(idx_node)):
        node = idx_node[i]
        node_idx[node] = i

    dist_row = []
    for node in node_children:
        if len(node_children[node]) <= 0:
            continue
        cnode1, cnode2 = sorted(list(node_children[node]))
        dist, size = node_dist[node], len(node_keys[node])
        row = [node_idx[cnode1], node_idx[cnode2], dist, size]
        dist_row.append((dist, row))
    dist_row = sorted(dist_row)
    Z = np.asarray([row for dist, row in dist_row])
    return idx_node, Z

# plot dendrogram
def plot_dendrogram(Z, idx_name, node_color=None, name_color=None):
    fig = plt.figure(figsize=(8, 8))
    if node_color != None:
        dendrogram(Z, link_color_func=lambda k: node_color[k], orientation='right')
    else:
        dendrogram(Z, orientation='right')
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
    plt.savefig("dendrogram.png", dpi=500, bbox_inches='tight')
    #plt.show()
    plt.close()


# get geneset name
def gs_name(gs):
    words = gs['name'].split('_')[1:]
    for k in range(len(words)):
        if words[k] in ['B', 'STAT']:
            continue
        words[k] = words[k].lower()
    words[0] = words[0][0].upper() + words[0][1:]
    name = ' '.join(words)
    return name



# load data
path1 = "NEW_GSEA_inht_GOBP"
path2 = "NEW_GSEA_KO_GOBP"

pos_gs_list1, neg_gs_list1 = read_GSEA(path1)
pos_gs_list2, neg_gs_list2 = read_GSEA(path2)


# select gene-sets with the biggest |nes|
num = 20
gs_list1 = select_best(pos_gs_list1 + neg_gs_list1, num)
gs_list2 = select_best(pos_gs_list2 + neg_gs_list2, num)


# clustering gene-sets based on the degree of sharing genes using Kappa-metric
#algorithm = 'Spectral'
#cluster_num = 3
#algorithm = 'OPTICS'

algorithm = 'Hierarchial'
cID_color = ['tab:orange', 'tab:green', 'tab:red']

gs_list =  gs_list1 + gs_list2

#gs_list = list(set(gs_list1) | set(gs_list2))
gene_list, kappa_matrix = get_kappa_matrix (gs_list)
#gene_list, kappa_matrix = get_kappa_matrix (gs_list, gene_list=gene_value.keys())

if algorithm == 'Spectral':
    score_matrix = kappa_to_score (kappa_matrix)
    idx_cID, cID_idxs = Spectral_clustering (score_matrix, cluster_num=cluster_num)

elif algorithm == 'Hierarchial':
    dist_matrix = kappa_to_dist (kappa_matrix)
    Z, idx_cID, cID_idxs = Hierarchial_clustering (dist_matrix)

elif algorithm == 'OPTICS':
    dist_matrix = kappa_to_dist (kappa_matrix)
    idx_cID, cID_idxs = OPTICS_clustering (dist_matrix)



name_list = []
name_color = {}
for gs, cID in zip(gs_list, idx_cID):
    name = gs_name(gs)
    name_list.append(name)
    name_color[name] = cID_color[cID]

# plot dendrogram after remvoing redundant nodes
if algorithm == 'Hierarchial':
    
    node_children, node_dist, node_names = decode_Z(Z, name_list)

    # prune redundant leaves of binary tree
    node_newnode = {}
    for node in sorted(node_children):
        # skip leaf nodes
        if len(node_children[node]) <= 0:
            continue
        cnode1, cnode2 = sorted(list(node_children[node]))
        # level-1 nodes
        try:
            name1, name2 = name_list[cnode1], name_list[cnode2]
            # the case of having redundant leaves: remove the root and one of leaves
            if name1 == name2:
                assert node_dist[node] == 0
                assert len(node_names[node]) == 1
                del node_children[node]
                del node_children[cnode2]
                del node_dist[node]
                del node_dist[cnode2]
                del node_names[node]
                del node_names[cnode2]
                node_newnode[node] = cnode1
        # level > 1 nodes
        except:
            # reconnect upper node to survived leaves
            try:
                newcnode1 = node_newnode[cnode1]
            except:
                newcnode1 = cnode1
            try:
                newcnode2 = node_newnode[cnode2]
            except:
                newcnode2 = cnode2
            del node_children[node]
            node_children[node] = set([newcnode1, newcnode2])
    
    # make new Z and plot dendrogram
    idx_node, newZ = encode_Z (node_children, node_dist, node_names)

    # set leaf name and node color
    new_name_list = []
    node_color = {}
    for k in range(len(idx_node)):
        node = idx_node[k]
        names = list(node_names[node])
        if len(names) == 1:
            name = names[0]
            new_name_list.append(name)

        colors = set([])
        for name in names:
            color = name_color[name]
            colors.add(color)
        colors = list(colors)

        if len(colors) == 1:
            node_color[k] = colors[0]
        else:
            node_color[k] = 'blue'
    
    plot_dendrogram(newZ, new_name_list, node_color, name_color)
    

    
#print ("%s clustering" % (algorithm))
#for cID, idxs in sorted(cID_idxs.items()):
#    print ([gs_list[i]['name'] for i in idxs])
#    print ()
#print ()
#print ()

sys.exit(1)

for gs_list, note in zip([gs_list1, gs_list2], [path1, path2]):

    min_value = min(gene_value.values())
    max_value = max(gene_value.values())
    lim_value = min(abs(min_value), abs(max_value))
    nes_list = [gs['nes'] for gs in gs_list]
    min_nes = min(nes_list)
    max_nes = max(nes_list)

    cmap = mpl.cm.get_cmap("binary")
    bg_cmap = mpl.cm.get_cmap("Spectral_r")
    #cmap = mpl.cm.get_cmap("bwr")
    #cmap = mpl.cm.get_cmap("seismic")
    #cmap = mpl.cm.get_cmap("coolwarm")
    #cmap = mpl.cm.get_cmap("Spectral_r")

    nrows = len(gs_list)+1
    ncols=2
    fig, axes = plt.subplots(nrows=nrows,
                             ncols=ncols,
                             #ncols=3,
                             #figsize=(5,5),
                             #figsize=(5,10),
                             #figsize=(3,10),
                             figsize=(3,5),
                             sharey=True,
                             #gridspec_kw={'width_ratios': [0.5, 1]})
                             gridspec_kw={'width_ratios': [0.6, 0.4]})
                             #gridspec_kw={'width_ratios': [0.5, 1, 0.5]})


    for i in range(nrows):
        for j in range(ncols):

            #if i == 0 and j == 1:
            if i == 0 and j == 0:
                # plot ranked values
                axes[i][j].bar(range(len(gene_value)), sorted(gene_value.values(), reverse=True), color='k')
                axes[i][j].set_xlim([0, len(gene_value)-1])

            if i > 0:
                gs = gs_list[i-1]
                nes = gs['nes']
                name = gs_name(gs)
                color = name_color[name]
                #cID = idx_cID[i-1]
                #color = cID_color[cID]

                #if j == 0:
                if j == 1:
                    # put annotation
                    #axes[i][j].annotate(name, (0,0), ha='right', va='center', weight='bold', fontsize=10)
                    #axes[i][j].set_xlim([-max_nes,0])
                    axes[i][j].annotate(name, (0,0), ha='left', va='center', weight='bold', color=color, fontsize=10)
                    axes[i][j].set_xlim([0, max_nes])

                    axes[i][j].set_ylim([-1,1])

                #elif j == 1:
                elif j == 0:

                    # plot background color gradient
                    binnum = 10
                    binsize = len(gene_value)/10
                    for k in range(binnum):
                        color = float(k+0.5)/binnum
                        axes[i][j].axvspan(k*binsize, (k+1)*binsize-1, color=bg_cmap(color), alpha=0.3, lw=0, zorder=1)
                    axes[i][j].axhspan(0, 1, color='white', alpha=1, lw=0, zorder=1)

                    # plot gene set ranks
                    gene_info = gs['genes']
                    pos_list = []
                    color_list = []
                    alpha_list = []
                    for gene in gene_info:
                        rank = gene_info[gene]['rank']
                        value = gene_value[gene]
                        pos_list.append(rank-1)
                        #color = rescale(value, old_st=-lim_value, old_ed=lim_value, new_st=0, new_ed=1) #color by value
                        color = rescale(rank, old_st=1, old_ed=len(gene_value), new_st=0, new_ed=1) #color by rank
                        color_list.append(color)
                        es = gene_info[gene]['es']
                        if es < 0:
                            alpha_list.append(-es)
                        else:
                            alpha_list.append(es)

                    alpha_list = statis.rescale(alpha_list, old_st=0, old_ed=max(alpha_list), new_st=0.05, new_ed=1)

                    for pos, color, alpha in zip(pos_list, color_list, alpha_list):
                        axes[i][j].axvline(x=pos, color=cmap(alpha), linewidth=0.25, alpha=1)
                        axes[i][j].set_xlim([0, len(gene_value)-1])
                        axes[i][j].set_ylim([-1, 1])

            # remove ticks
            axes[i][j].tick_params(top='off', bottom='off', left='off', right='off',
                                   labelleft='off', labelright='off', labeltop='off', labelbottom='off')
            #if i <= 0 or j != 1:
            if i <= 0 or j != 0:
                # remove frames
                axes[i][j].spines['top'].set_visible(False)
                axes[i][j].spines['bottom'].set_visible(False)
                axes[i][j].spines['left'].set_visible(False)
                axes[i][j].spines['right'].set_visible(False)

    plt.subplots_adjust(wspace=0.03)
    #plt.savefig("GSEA.svg", format='svg', bbox_inches='tight')
    plt.savefig("GSEA" + "_" + note + ".png", dpi=500, bbox_inches='tight')
    #plt.show()
    plt.close()

        
"""
for j in range(ncols):
    if j ==1:
        axes[0][j].bar(range(len(gene_value)), sorted(gene_value.values(), reverse=True), color='tab:gray')
    axes[0][j].spines['top'].set_visible(False)
    axes[0][j].spines['bottom'].set_visible(False)
    axes[0][j].spines['left'].set_visible(False)
    axes[0][j].spines['right'].set_visible(False)
    axes[0][j].tick_params(top='off', bottom='off', left='off', right='off',
           labelleft='off', labelright='off', labeltop='off', labelbottom='off')
                         

for i in range(len(gs_list)):
    gs = gs_list[i]
    nes = gs['nes']
    name = ' '.join(gs['name'].split('_')[1:])

    # plot NES bar graph
    #if nes >= 0:
        #axes[i][0].barh([0],[-nes], height=1, color='tab:blue', alpha=0.5)
        #axes[i][0].annotate(str(nes), (0,-nes), ha='center', va='center')
        #axes[i][0].annotate(name, (0,0), ha='right', va='center', weight='bold', fontsize=10)
        #axes[i][0].set_xlim([-max_nes,0])
        #axes[i][0].set_ylim([-1,1])
    #else:
    #    axes[i][2].barh([0],[-nes], height=1, color='tab:red', alpha=0.5)
    #    axes[i][2].set_xlim([0, -min_nes])
    #    axes[i][2].set_ylim([-1,1])
    
    axes[i+1][0].annotate(name, (0,0), ha='right', va='center', weight='bold', fontsize=10)
    axes[i+1][0].set_xlim([-max_nes,0])
    axes[i+1][0].set_ylim([-1,1])

    axes[i+1][0].spines['top'].set_visible(False)
    axes[i+1][0].spines['bottom'].set_visible(False)
    axes[i+1][0].spines['left'].set_visible(False)
    axes[i+1][0].spines['right'].set_visible(False)

    axes[i+1][0].tick_params(top='off', bottom='off', left='off', right='off',
                       labelleft='off', labelright='off', labeltop='off', labelbottom='off')

    # put nes value
    
    #if nes >=0:
    #    barcolor='tab:blue'
    #else:
    #    barcolor='tab:red'
    
    #axes[i][2].barh([0],[nes], height=1, color=barcolor, alpha=0.5)
    #axes[i][2].set_ylim([-1,1])
    #axes[i][2].set_xlim([min_nes, max_nes])

    #axes[i][2].spines['top'].set_visible(False)
    #axes[i][2].spines['bottom'].set_visible(False)
    #axes[i][2].spines['left'].set_visible(False)
    #axes[i][2].spines['right'].set_visible(False)
        
    #axes[i][2].tick_params(top='off', bottom='off', left='off', right='off',
    #               labelleft='off', labelright='off', labeltop='off', labelbottom='off')


    # plot backgroun color
    binnum = 10
    binsize = len(gene_value)/10
    for k in range(binnum):
        color = float(k+0.5)/binnum
        axes[i][1].axvspan(k*binsize, (k+1)*binsize-1, color=bg_cmap(color), alpha=0.3, lw=0, zorder=1)
    axes[i][1].axhspan(0, 1, color='white', alpha=1, lw=0, zorder=1)

    # plot gene set ranks
    gene_info = gs['genes']
    pos_list = []
    color_list = []
    alpha_list = []
    for gene in gene_info:
        rank = gene_info[gene]['rank']
        value = gene_value[gene]
        pos_list.append(rank-1)
        #color = rescale(value, old_st=-lim_value, old_ed=lim_value, new_st=0, new_ed=1) #color by value
        color = rescale(rank, old_st=1, old_ed=len(gene_value), new_st=0, new_ed=1) #color by rank
        color_list.append(color)
        es = gene_info[gene]['es']
        if es < 0:
            alpha_list.append(-es)
        else:
            alpha_list.append(es)

    alpha_list = statis.rescale(alpha_list, old_st=0, old_ed=max(alpha_list), new_st=0.05, new_ed=1)
        
    for pos, color, alpha in zip(pos_list, color_list, alpha_list):
        #axes[i][1].axvline(x=pos, color='k', linewidth=0.25, alpha=alpha)
        #axes[i][1].axvline(x=pos, color=cmap(color), linewidth=0.25, alpha=alpha)
        axes[i+1][1].axvline(x=pos, color=cmap(alpha), linewidth=0.25, alpha=1)
        #axes[i][1].axvline(x=pos, color=cmap(alpha), linewidth=0.3, alpha=1)
        axes[i+1][1].set_xlim([0, len(gene_value)-1])
        axes[i+1][1].set_ylim([-1, 1])
        axes[i+1][1].tick_params(top='off', bottom='off', left='off', right='off',
                               labelleft='off', labelright='off', labeltop='off', labelbottom='off')

plt.subplots_adjust(wspace=0.03)
plt.savefig("GSEA.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()
"""
