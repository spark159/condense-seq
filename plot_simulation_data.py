import sys
import copy
import numpy as np
import matplotlib.pyplot as plt
import load_file
import statis
import sklearn
import sklearn.cluster


# convert contact matrix to score matrix
def contact_to_score (contact_matrix, scale=1.0):
    nrow, ncol = contact_matrix.shape
    score_matrix = np.zeros((nrow, ncol))
    for i in range(nrow):
        for j in range(ncol):
            prob = contact_matrix[i][j]
            score_matrix[i][j] = np.exp(-scale*(1.0-prob))
    return score_matrix

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


### parameters
path = "/home/spark159/../../storage/"
agent_pathname = {'sp':'spermine',
                  'HP1a':'HP1a'}

chr_choice = 'chr1'
region = (14, 64)
beadnum = 2000
beadsize = 25000

exp_list = [('H1', 'NCP', 'sp', 8),
            ('GM12878', 'NCP', 'sp', 8),
            ('H1', 'NCP', 'HP1a', 3)]

#exp_list = [('H1', 'NCP', 'sp', 8)]

# load simulation contact probabiltiy matrix
exp_cmatrix = {}
for exp in exp_list:
    cell, sample, agent, tnum  = exp

    subpath = '_'.join([cell,
                        'hg38',
                        'WT',
                        agent_pathname[agent],
                        'tp' + str(tnum),
                        chr_choice,
                        '%dM-%dM' % region,
                        str(beadnum) + 'beads',
                        str(beadsize) + 'bp']) + '/'

    
    score_fname = path + subpath + "condensability_scores.npy"
    exp_fname = path + subpath + "experiment_contact_map.npz"
    sim_fname = path + subpath + "simulation_contact_map.npz"
    analysis_fname = path + subpath + "analysis_results.npz"

    #scores = np.load(score_fname)
    #exp_data = np.load(exp_fname)
    sim_data = np.load(sim_fname)

    exp_cmatrix[exp] = sim_data['cmap']


# clustering based on contact probability
cluster_num = 5

for exp in exp_list:
    cell, sample, agent, tnum = exp
    
    contact_matrix = exp_cmatrix[exp]
    score_matrix = contact_to_score (contact_matrix)
    idx_cID, cID_idxs = Spectral_clustering (score_matrix, cluster_num=cluster_num)

    fig, axes = plt.subplots(nrows=2,
                             ncols=2,
                             figsize=(10,10),
                             gridspec_kw={'height_ratios':[0.1, 2],
                                          'width_ratios':[0.1, 2]})

    #axes[0,1].plot(range(len(idx_cID)), idx_cID)
    axes[0,1].imshow([idx_cID], aspect='auto')
    axes[0,1].set_xlim([0, len(idx_cID)])

    axes[1,0].imshow([[cID] for cID in idx_cID], aspect='auto')
    axes[1,0].set_ylim([0, len(idx_cID)])
    
    axes[1,1].imshow(np.log(contact_matrix), aspect='auto')
    axes[1,1].set_xlim([0, len(idx_cID)])
    axes[1,1].set_ylim([0, len(idx_cID)])
    #plt.subplots_adjust(wspace=0.02, hspace=0.02)
    plt.suptitle(exp)
    plt.show()
    plt.close()
    


    

    

    
    
