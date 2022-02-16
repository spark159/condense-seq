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

# read loop-seq data file
def read_loop_seq (fname):
    field_ID_value = {}
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            fields = ['score', 'amplt', 'phase', 'seq']
            First=False
            continue
        ID, seq = int(cols[0]), cols[1]
        values = [float(value) for value in cols[-3:]] + [seq]
        for field, value in zip(fields, values):
            if field not in field_ID_value:
                field_ID_value[field] = {}
            field_ID_value[field][ID] = value
    return field_ID_value

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
        seq = seq.upper()
        dist_count = {k:0 for k in range(klen, min(len(seq)-klen, max_dist)+1)}
        for i in range(len(seq)-klen):
            for j in range(i+1, min(len(seq)-klen+1, i+max_dist+1)):
                dist = j-i
                if dist < klen:
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


# load loop-seq data
path = ''
fname = 'loop_seq_random_library.txt'
field_ID_value = read_loop_seq(path+fname)
ID_score = field_ID_value['score']
ID_seq = field_ID_value['seq']
NCPlen = 100
clip = 25 # clip off the both ends of sequence 
max_dist = 48
note = '_loop'
vmin, vmax = -0.15, 0.15


# load condense-seq data 
#path = ''
#fname = "H1_NCP_sp_chr1_anot.cn"
#ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+fname)
#ID_score = name_ID_value["work/2021_06_07_H1_sp_detail/H1-NCP-sp-8"]
#ID_seq = name_ID_value['Sequence']
#NCPlen = 147
#clip = 10 # clip off the both ends of sequence 
#max_dist = 48
#note = '_condense'
#vmin, vmax = -0.05, 0.05


# get mean kmer-pair occurrence for random sequences
path = ""
fname = "random_pair_dist_count" + note
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
    random_pair_dist_count = get_pair_dist_count(random_seq_list, max_dist=max_dist, norm=True)
    pickle.dump(random_pair_dist_count, open(fname + ".pickle", "wb"))
    del random_seq_list
print "random sequence counting is done"


# get kmer-pair correlation function and helicity score for library
# save the correlation between hscore and loop/consensability score
path = ""
fname = "pair_scorr" + note
try:
    pair_scorr = pickle.load(open(fname + ".pickle", "rb"))
except:
    print "start reading sequence data"
    print "total %d sequences" % (len(ID_seq))
    count = 0 
    log_counter = -1
    #ID_pair_corrf = {}
    pair_ID_hscore = {}
    for ID, seq in ID_seq.items():
        # adjust sequence
        middle = len(seq)/2
        seq = seq[middle-NCPlen/2:middle+NCPlen/2+1]
        assert len(seq) == NCPlen
        seq = seq[clip:len(seq)-clip].upper()
        if set(seq) > set('ATCG'):
            continue
        # counting the sequences
        count +=1
        if int(np.log10(count)) > log_counter:
            print count
            log_counter +=1
    
        pair_dist_count = get_pair_dist_count([seq], max_dist=max_dist)
        pair_corrf = {}
        for pair in pair_dist_count:
            corrf = {}
            for dist in pair_dist_count[pair]:
                norm = random_pair_dist_count[pair][dist]
                corr = float(pair_dist_count[pair][dist]) / norm
                corrf[dist] = corr
            pair_corrf[pair] = corrf    

            if pair not in pair_ID_hscore:
                pair_ID_hscore[pair] = {}
            assert ID not in pair_ID_hscore
            pair_ID_hscore[pair][ID] = get_hscore(corrf)

        #ID_pair_corrf[ID] = pair_corrf
    print "get helicity score is done"

    # check the correlation between hscore and condensability
    pair_scorr = {}
    for pair in pair_ID_hscore:
        X, Y = [], []
        for ID in ID_seq:
            try:
                X.append(pair_ID_hscore[pair][ID])
            except:
                continue
            Y.append(ID_score[ID])
        pair_scorr[pair] = statis.get_spearman_corr(X, Y)

    del pair_ID_hscore

    pickle.dump(pair_scorr, open(fname + ".pickle", "wb"))
print "get helicity score VS condensabiltiy correlation is done"        


# plot kmer-pair correlation matrix
GCcontent_din = sorted([(GC_content(din), din) for din in all_path(2, 'ATCG')])
all_din = [din for GCcontent, din in GCcontent_din]

corr_img = np.zeros((len(all_din), len(all_din)))
for i in range(len(all_din)):
    for j in range(i, len(all_din)):
        din1, din2 = all_din[i], all_din[j]
        pair = tuple(sorted([din1, din2]))
        corr_img[i][j] = pair_scorr[pair]
        corr_img[j][i] = pair_scorr[pair]

fig = plt.figure(figsize=(3,3))
plt.imshow(corr_img, cmap='Spectral_r', vmin=vmin, vmax=vmax)
plt.xticks(range(len(all_din)), all_din, fontsize=8, rotation=48)
plt.yticks(range(len(all_din)), all_din, fontsize=8)
if note.endswith('loop'):
    title = "Dinucleotide Helicity VS %s" % ('Cyclizability')
elif note.endswith('condense'):
    title = "Dinucleotide Helicity VS %s" % ('Condensability')
plt.title(title, fontsize=8)
cbar = plt.colorbar(pad=0.05, ticks=[vmin, vmax], fraction=0.04)
cbar.ax.set_yticklabels([str(vmin), str(vmax)], fontsize=5)
cbar.ax.set_ylabel('Spearman correlation', rotation=-90, va="bottom", labelpad=-15, fontsize=8)
plt.savefig(title + '.svg', format='svg', bbox_inches='tight')
#plt.show()
plt.close()
