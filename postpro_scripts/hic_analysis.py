import sys
import math
import random
import copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import load_file
import scipy as sc
from scipy.linalg import eig
import matplotlib.gridspec as gridspec
import h5py
import cooler

def printname(name):
    print name

f = h5py.File("/home/spark159/../../media/spark159/sw/dataforcondense/4DNFIB59T7NN.mcool", 'r')

rough_data = f['resolutions/500000']

"""
target_ids = []
i = 0
while True:
    if rough_data['bins']['chrom'][i] != 0:
        break
    target_ids.append(i)
    i +=1

N = len(target_ids)
rough_matrix = np.zeros((N, N))

for i in range(len(rough_data['pixels']['bin1_id'][:])):
    id1 = rough_data['pixels']['bin1_id'][i]
    id2 = rough_data['pixels']['bin2_id'][i]
    weight1 = rough_data['bins']['weight'][id1]
    weight2 = rough_data['bins']['weight'][id2]
    if id1 not in target_ids:
        break
    if id2 not in target_ids:
        continue
    count = rough_data['pixels']['count'][i]
    rough_matrix[id1][id2] += count * weight1 * weight2
    rough_matrix[id2][id1] += count * weight2 * weight1

fig = plt.figure()
plt.imshow(np.log10(rough_matrix), cmap="Reds")
plt.show()
plt.close()



scale = 0.413 # total count ratio of test sample to control (UV data, sp9)
#scale = 0.154 # total count ratio of test sample to control (UV data, sp10)

def reaction_prob (score1, score2, metric, scale=scale):
    prob1 = 1.0 - scale*np.exp(-score1)
    if prob1 < 0:
        prob1 = 0.001
    if prob1 > 1:
        prob1 = 0.999
    #assert prob1 >=0 and prob1 <= 1
    prob2 = 1.0 - scale*np.exp(-score2)
    if prob2 < 0:
        prob2 = 0.001
    if prob2 > 1:
        prob2 = 0.999
    #assert prob2 >=0 and prob2 <= 1
    if metric == "product":
        prob = prob1*prob2
    elif metric == "GM":
        prob = np.sqrt(prob1*prob2)
    return prob

# load annotation file
path = "/home/spark159/../../media/spark159/sw/dataforcondense//hg19_chr1_171_everything_anot.cn"
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path)
ID_score = name_ID_value['work/condense_seq/sp9_hg19_chr1']
ID_AT = name_ID_value['ATcontent']

temp = {}
NCPst = 0
NCPnum = 1000
for i in range(NCPst, NCPst+NCPnum):
    temp[i] = ID_score[i]
ID_score = temp

IDs = sorted(ID_score.keys())

# binning the genome and create contact matrix
bin_size = 1000
bin_pair_probs = {}
for i in range(len(IDs)-1):
    for j in range(i+1, len(IDs)):
        ID1, ID2 = IDs[i], IDs[j]
        prob = reaction_prob (ID_score[ID1], ID_score[ID2], "product")
        pos1, pos2 = ID_pos[ID1], ID_pos[ID2]
        idx1, idx2 = pos1 / bin_size, pos2 / bin_size
        bin_pair = (idx1, idx2)
        if bin_pair not in bin_pair_probs:
            bin_pair_probs[bin_pair] = []
        bin_pair_probs[bin_pair].append(prob)

bin_num = ID_pos[IDs[-1]] / bin_size + 1
mean_matrix = np.zeros((bin_num, bin_num))
for bin_pair in bin_pair_probs:
    try:
        mean_prob = np.mean(bin_pair_probs[bin_pair])
    except:
        continue
    idx1, idx2 = bin_pair
    mean_matrix[idx1][idx2] = mean_prob
    mean_matrix[idx2][idx1] = mean_prob

fig = plt.figure()
plt.imshow(mean_matrix, cmap='Reds')
plt.colorbar()
plt.title("Mean contact matrix")
plt.show()
plt.close()


## normalize matrix
#values = np.sum(mean_matrix, axis=0)
#threshold = sorted(values)[int(bin_num*0.25)]
##threshold = 0.0

#excludes = []
##for i in range(len(values)):
##    if values[i] <= threshold:
##        excludes.append(i)

#cycle_num = 10
#norm_matrix = copy.deepcopy(mean_matrix)
#for k in range(cycle_num):
#    for i in range(bin_num):
#        norm_matrix[i,excludes] = 0.0
#        total = np.sum(norm_matrix[i,:])
#        if total > 0:
#            norm_matrix[i,:] = norm_matrix[i,:] / total
#    for j in range(bin_num):
#        norm_matrix[excludes,j] = 0.0
#        total = np.sum(norm_matrix[:,j])
#        if total > 0:
#            norm_matrix[:,j] = norm_matrix[:,j] / total

#norm_matrix = (norm_matrix + norm_matrix.T) / 2.0

#fig = plt.figure()
#plt.imshow(norm_matrix, cmap='Reds')
#plt.colorbar()
#plt.title("Normalized contact matrix")
#plt.show()
#plt.close()

norm_matrix = copy.deepcopy(mean_matrix)

# compute distance vs mean contact frequency
dbin_values = [ [] for i in range(bin_num) ]
for i in range(bin_num):
    for j in range(i, bin_num):
        dbin = (j-i)
        value = norm_matrix[i][j]
        dbin_values[dbin].append(value)

dist_list, mcontact_list = [], []
for i in range(len(dbin_values)):
    dist = i*bin_size
    mcontact = np.mean(dbin_values[i])
    dist_list.append(dist)
    mcontact_list.append(mcontact)

fig = plt.figure()
plt.plot(dist_list, mcontact_list)
plt.xscale("log")
plt.yscale("log")
plt.title("Distance VS mean contact frequency")
plt.xlabel("Genomic distance (bp)")
plt.ylabel("Mean contact frequency")
plt.show()
plt.close()

# compute observed/expected matrix
dist_matrix = np.zeros((bin_num, bin_num))
for i in range(bin_num):
    for j in range(i, bin_num):
        dbin = (j-i)
        mvalue = np.mean(dbin_values[dbin])
        dist_matrix[i][j] = mvalue
        dist_matrix[j][i] = mvalue

fig = plt.figure()
plt.imshow(dist_matrix, cmap='Reds')
plt.colorbar()
plt.title("expected contact matrix")
plt.show()
plt.close()

fold_matrix = norm_matrix / dist_matrix
fold_matrix[np.isnan(fold_matrix)] = 0.0

fig = plt.figure()
plt.imshow(fold_matrix)
plt.colorbar()
plt.title("Observed/expected contact matrix")
plt.show()
plt.close()

# compute correlation matrix
def get_corr(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = np.average(x)
    avg_y = np.average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff
    return diffprod / np.sqrt(xdiff2 * ydiff2)

corr_matrix = np.zeros((bin_num, bin_num))
for i in range(bin_num):
    for j in range(i, bin_num):
        X, Y = fold_matrix[i,:], fold_matrix[j,:]
        corr = get_corr(X, Y)
        if np.isnan(corr):
            corr = 0.0
        corr_matrix[i][j] = corr
        corr_matrix[j][i] = corr

fig = plt.figure()
plt.imshow(corr_matrix)
plt.colorbar()
plt.title("correlation matrix")
plt.show()
plt.close()

# assign A/B compartment
V,D = sc.linalg.eig(corr_matrix)
eigenvect = D[:,0]
gs = gridspec.GridSpec(2, 1, height_ratios=[5,1])
ax1 = plt.subplot(gs[0])
ax1.imshow(norm_matrix)
ax1.set_aspect("auto")
ax1.set_title("normalized matrix")

ax2 = plt.subplot(gs[1], sharex = ax1)
ax2.set_xlim([0, len(eigenvect)])
ax2.fill_between(range(0, len(eigenvect)), 0, eigenvect, color='blue')
ax2.fill_between(range(0, len(eigenvect)), 0, eigenvect, eigenvect<0, color='red')
plt.show()
"""
