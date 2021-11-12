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
from scipy import stats
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.cm as cm
from scipy import fft
import matplotlib as mpl
import random
import pickle

# get GC content of sequence
def GC_content(seq, percent=False):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    if percent:
        return (num/float(len(seq)))*100
    return (num/float(len(seq)))

# get reverse complementary of sequence
def rev_comp (seq):
    rev_seq=seq[::-1]; new_seq=''
    for nt in rev_seq:
        nt = nt.upper()
        if nt == 'A':
            new_seq += 'T'
        if nt == 'T':
            new_seq += 'A'
        if nt == 'C':
            new_seq += 'G'
        if nt == 'G':
            new_seq += 'C'
    return new_seq

# get all possible N-mers
def all_path (N, states='ATCG'):
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output

# count kmer-pair occurrence
def get_pair_dist_count (seq_list, klen=2, max_dist=sys.maxint, rev=True, norm=False):
    pair_dist_count = {}
    
    for seq in seq_list:
        dist_count = {k:0 for k in range(1, min(len(seq)-klen, max_dist)+1)}
        for i in range(len(seq)-klen):
            for j in range(i+1, min(len(seq)-klen+1, i+max_dist+1)):
                #dist = j-i
                dist = j - i - 1
                if dist <= 0:
                    continue
                pair_list = []

                kmer1 = seq[i:i+klen]
                kmer2 = seq[j:j+klen]
                pair = tuple(sorted([kmer1, kmer2]))
                pair_list.append(pair)

                if rev:
                    kmer1 = rev_comp(kmer1)
                    kmer2 = rev_comp(kmer2)
                    pair = tuple(sorted([kmer1, kmer2]))
                    pair_list.append(pair)

                for pair in pair_list:
                    if pair not in pair_dist_count:
                        pair_dist_count[pair] = copy.deepcopy(dist_count)
                    if norm:
                        pair_dist_count[pair][dist] +=1.0/len(seq_list)
                    else:
                        pair_dist_count[pair][dist] +=1

    return pair_dist_count

# pair helicity score
def get_hscore (corrf, inphases=[10, 20, 30], outphases=[5, 15, 25], offset=1):
    def phase_score (phases):
        score = 0.0
        for phase in phases:
            corrs = []
            for dist in range(phase-offset, phase+offset+1):
                corrs.append(corrf[dist])
            score += max(corrs)
        return score
    return phase_score(inphases) - phase_score(outphases)
            
                

# basic parameters
NCPlen = 147
clip = 10 # clip off the both ends of sequence 
max_dist = 48

# get mean kmer-pair occurrence for random sequences
path = ""
fname = "random_pair_dist_count"
try:
    random_pair_dist_count = pickle.load(open(fname + ".pickle", "rb"))
except:
    random.seed(0)
    random_num = 10000
    random_seq_list = []
    for i in range(random_num):
        seq = ""
        for k in range(NCPlen-2*clip):
            seq += random.choice('ATCG')
        random_seq_list.append(seq)
    random_pair_dist_count = get_pair_dist_count(random_seq_list, max_dist=48, norm=True)
    pickle.dump(random_pair_dist_count, open(fname + ".pickle", "wb"))
    del random_seq_list
print "random sequence counting is done"

# load annotation file
path = ''
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+"H1_NCP_sp_chr1_anot.cn")
ID_score2 = name_ID_value["work/2021_06_07_H1_sp_detail/H1-NCP-sp-8"]
ID_seq = name_ID_value['Sequence']

# get kmer-pair correlation function and helicity score for library
ID_pair_corrf = {}
pair_ID_hscore = {}
for ID, seq in ID_seq.items():
    seq = seq[clip:len(seq)-clip]
    pair_dist_count = get_pair_dist_count([seq], max_dist=48)
    pair_corrf = {}
    for pair in pair_dist_count:
        corrf = {}
        for dist in pair_dist_count[pair]:
            norm = random_pair_dist_count[pair][dist]
            corr = pair_dist_count[pair][dist] / norm
            corrf[dist] = corr
        pair_corrf[pair] = corrf    

        if pair not in pair_ID_hscore:
            pair_ID_hscore[pair] = {}
        assert ID not in pair_ID_hscore
        pair_ID_hscore[pair][ID] = get_hscore(corrf)
            
    ID_pair_corrf[ID] = pair_corrf
    hscore = get_pair_hscore(pair_corr)

del pair_dist_count
del random_pair_dist_count
print "get pair-correlation is done"

# check the correlation between hscore and condensability
pair_scorr = {}
for pair in pair_ID_hscore:
    ID_hscore = pair_ID_hscore[pair]
    X = [ID_hscore[ID] for ID in ID_pos]
    Y = [ID_score2[ID] for ID in ID_pos]
    pair_scorr[pair] = statis.get_spearman_corr(X, Y)


# plot kmer-pair correlation matrix
GCcontent_din = sorted([(GC_content(din), din) for din in all_path(2, 'ATCG')], cmp=tuple_cmp)
all_din = [din for GCcontent, din in GCcontent_din]

corr_img = np.zeros((len(all_din), len(all_din)))
for i in range(len(all_din)):
    for j in range(i, len(all_din)):
        din1, din2 = all_din[i], all_din[j]
        pair = tuple(sorted([din1, din2]))
        corr_img[i][j] = pair_scorr[pair]
        corr_img[j][i] = pair_scorr[pair]

fig = plt.figure()
plt.imshow(corr_img, cmap='Spectral', vmin=-1, vmax=1)
plt.xticks(range(len(all_din)), all_din)
plt.yticks(range(len(all_din)), all_din)
plt.title("Dinucleotide Helicity VS Condensability")
plt.colorbar()
plt.show()
plt.close()
