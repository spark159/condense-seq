#import load_file
#import graphics
#import statis
import sys
import copy
import pickle
#import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
import random
#import Interval_dict
from scipy import stats
from mpl_toolkits.axes_grid1 import ImageGrid
#from SliderClass import Slider

def rev_comp (seq):
    rev_seq=seq[::-1]; new_seq=''
    for nt in rev_seq.upper():
        if nt == 'A':
            new_seq += 'T'
        elif nt == 'T':
            new_seq += 'A'
        elif nt == 'C':
            new_seq += 'G'
        elif nt == 'G':
            new_seq += 'C'
        else:
            new_seq += nt
    return new_seq

def all_path (N, states='ATCG'):
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output

def get_kmer_probs (seq_list, weight_list, bscore_list, seq_len, knum, step, bnum, rev=False):
    # initialize
    pos_num = (seq_len - knum + 1) / step
    all_kmers = all_path(knum, 'ATCG')

    kmer_prob = {kmer:0.0 for kmer in all_kmers}
    kmer_bscore_prob = {}
    for kmer in all_kmers:
        for b in range(bnum):
            kmer_bscore_prob[(kmer, b)] = 0.0
            
    bscore_prob = {b:0.0 for b in range(bnum)}
    pos_kmer_prob = [copy.deepcopy(kmer_prob) for i in range(pos_num)]
    pos_kmer_bscore_prob = [copy.deepcopy(kmer_bscore_prob) for i in range(pos_num)]
                
    # counting
    for seq, weight, bscore in zip(seq_list, weight_list, bscore_list):
        seq = seq.upper()
        if rev:
            rev_seq = rev_comp(seq)
        for i in range(pos_num):
            kmer = seq[i*step:i*step+knum]
            try:
                pos_kmer_prob[i][kmer] += weight
                pos_kmer_bscore_prob[i][(kmer, bscore)] += weight
            except:
                continue
            if rev:
                rev_kmer = rev_seq[i*step:i*step+knum]
                try:
                    pos_kmer_prob[i][rev_kmer] += weight
                    pos_kmer_bscore_prob[i][(rev_kmer, bscore)] += weight
                except:
                    continue            
        bscore_prob[bscore] += weight

    # normalize
    total = float(sum(bscore_prob.values()))
    for bscore in bscore_prob:
        bscore_prob[bscore] = bscore_prob[bscore] / total

    for kmer_prob in pos_kmer_prob:
        total = float(sum(kmer_prob.values()))
        for kmer in kmer_prob:
            kmer_prob[kmer] = kmer_prob[kmer] / total

    for kmer_bscore_prob in pos_kmer_bscore_prob:
        total = float(sum(kmer_bscore_prob.values()))
        for kmer in kmer_bscore_prob:
            kmer_bscore_prob[kmer] = kmer_bscore_prob[kmer] / total

    return bscore_prob, pos_kmer_prob, pos_kmer_bscore_prob

def get_pos_MI (bscore_prob, pos_kmer_prob, pos_kmer_bscore_prob, seq_len, knum, step, bnum):
    def get_MI (kmer_bscore_prob, kmer_prob, bscore_prob):
        MI = 0.0
        for kmer, bscore in kmer_bscore_prob:
            joint_prob = kmer_bscore_prob[(kmer, bscore)]
            prob1, prob2 = kmer_prob[kmer], bscore_prob[bscore]
            if joint_prob <= 0 or prob1 <=0 or prob2 <= 0:
                continue
            MI += joint_prob*np.log2(joint_prob/(prob1*prob2))
        return MI

    pos_num = (seq_len - knum + 1) / step

    pos_MI = []
    for i in range(pos_num):
        kmer_bscore_prob = pos_kmer_bscore_prob[i]
        kmer_prob = pos_kmer_prob[i]
        MI = get_MI (kmer_bscore_prob, kmer_prob, bscore_prob)
        pos_MI.append(MI)

    return pos_MI


# read score file
def read_score (fname):
    First = True
    tname_ID_score = {} 
    for line in open(fname):
        cols = line.strip().split()
        if First:
            _, tname_list = cols[0], cols[1:]
            First = False
            continue
        
        ID, score_list = cols[0], cols[1:]
        ID = int(ID)

        for i in range(len(tname_list)):
            tname = tname_list[i]
            try:
                tname = int(tname)
            except:
                pass
            score = float(score_list[i])

            if tname not in tname_ID_score:
                tname_ID_score[tname] = {}
            tname_ID_score[tname][ID] = score
    return tname_ID_score

# read reference sequence
def read_ref (ref_fname):
    id_seq = {}
    for line in open(ref_fname):
        line = line.strip()
        if line.startswith('>'):
            id = int(line[4:])
            continue
        if line:
            assert id not in id_seq
            id_seq[id] = line
    return id_seq


tname_ID_score = read_score("YWlib_sp_score.txt")
ID_score = tname_ID_score[3]
ID_seq = read_ref('Blib.ref')


# binning the score
med = np.median(ID_score.values())
std = np.std(ID_score.values())/4
lines = [med-0.5*std-i*std for i in range(12)] + [med+0.5*std+i*std for i in range(12)]
lines = sorted(lines)
p_num = len(lines)+1
p_range = []
for i in range(p_num):
    if i == 0:
        st = -np.inf
        ed = lines[i]
    elif i == p_num-1:
        st = ed
        ed = np.inf
    else:
        st = ed
        ed = lines[i]
    p_range.append((st, ed))

p_IDs = [[] for i in range(p_num)]
for ID in ID_score:
    score = ID_score[ID]
    for i in range(p_num):
        st, ed = p_range[i]
        if score >= st and score < ed:
            break
    p_IDs[i].append(ID)

ID_bin = {}
for i in range(p_num):
    IDs = p_IDs[i]
    for ID in IDs:
        ID_bin[ID] = i

# pos-kmer x score mutual information analysis
bnum = p_num
knum = 3
step = 1
seq_len = 101

seq_list = []
weight_list = []
bscore_list = []
for ID in ID_seq.keys():
    seq = ID_seq[ID]
    middle = len(seq)/2
    seq = seq[middle-seq_len/2:middle+seq_len/2+1]
    seq_list.append(seq)
    weight_list.append(1.0)
    bscore_list.append(ID_bin[ID])
    

print "Data reading is done"

# for control
repeat = 20
sample_size = len(ID_seq)
pos_MI0 = []
for i in range(repeat):
    print i
    cbscore_list = [random.randrange(bnum) for k in range(sample_size)]
    bscore_prob, pos_kmer_prob, pos_kmer_bscore_prob = get_kmer_probs (seq_list, weight_list, cbscore_list, seq_len, knum, step, bnum, rev=True)

    pos_MI = get_pos_MI (bscore_prob, pos_kmer_prob, pos_kmer_bscore_prob, seq_len, knum, step, bnum)
    if i == 0:
        pos_MI0 = copy.deepcopy(pos_MI)
    else:
        for k in range(len(pos_MI)):
            pos_MI0[k] += pos_MI[k]

for k in range(len(pos_MI0)):
    pos_MI0[k] /= float(repeat)

print "MI calculation (control) is done"

# plot control MI profile
fig = plt.figure()
plt.plot(pos_MI0)
plt.show()
plt.close()


# for data
bscore_prob, pos_kmer_prob, pos_kmer_bscore_prob = get_kmer_probs (seq_list, weight_list, bscore_list, seq_len, knum, step, bnum, rev=True)
print "prob calculation is done"

pos_MI = get_pos_MI (bscore_prob, pos_kmer_prob, pos_kmer_bscore_prob, seq_len, knum, step, bnum)
print "MI calculation is done"

pos_MIfold = [pos_MI[pos]/pos_MI0[pos] for pos in range(len(pos_MI0))]


# plot enrichment MI profile
fig = plt.figure()
plt.plot(pos_MIfold)
plt.show()
plt.close()
