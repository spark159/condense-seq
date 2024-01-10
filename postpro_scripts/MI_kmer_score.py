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
import random
import Interval_dict
from scipy import stats
from mpl_toolkits.axes_grid1 import ImageGrid
from SliderClass import Slider

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

def get_kmer_probs (seq_list, weight_list, bscore_list, seq_len, knum, step, bnum, sym=False):
    # initialize
    pos_num = (seq_len - knum + 1) / step
    all_kmers = all_path(knum, 'ATCG')

    kmer_prob = {kmer:0.0 for kmer in all_kmers}
    kmer_bscore_prob = {}
    for kmer in all_kmers:
        for b in range(bnum):
            kmer_bscore_prob[(kmer, b)] = 0.0
            
    bscore_prob = {b:0.0 for b in range(bnum)}
    kmer_prob_list = [copy.deepcopy(kmer_prob) for i in range(pos_num)]
    kmer_bscore_prob_list = [copy.deepcopy(kmer_bscore_prob) for i in range(pos_num)]
                
    # counting
    for seq, weight, bscore in zip(seq_list, weight_list, bscore_list):
        seq = seq.upper()
        for i in range(pos_num):
            kmer = seq[i*step:i*step+knum]
            try:
                kmer_prob_list[i][kmer] += weight
                kmer_bscore_prob_list[i][(kmer, bscore)] += weight
            except:
                continue
        bscore_prob[bscore] += weight

    # normalize
    total = float(sum(bscore_prob.values()))
    for bscore in bscore_prob:
        bscore_prob[bscore] = bscore_prob[bscore] / total

    for kmer_prob in kmer_prob_list:
        total = float(sum(kmer_prob.values()))
        for kmer in kmer_prob:
            kmer_prob[kmer] = kmer_prob[kmer] / total

    for kmer_bscore_prob in kmer_bscore_prob_list:
        total = float(sum(kmer_bscore_prob.values()))
        for kmer in kmer_bscore_prob:
            kmer_bscore_prob[kmer] = kmer_bscore_prob[kmer] / total

    return bscore_prob, kmer_prob_list, kmer_bscore_prob_list

def get_MI_list (bscore_prob, kmer_prob_list, kmer_bscore_prob_list, seq_len, knum, step, bnum):
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

    MI_list = []
    for i in range(pos_num):
        kmer_bscore_prob = kmer_bscore_prob_list[i]
        kmer_prob = kmer_prob_list[i]
        MI = get_MI (kmer_bscore_prob, kmer_prob, bscore_prob)
        MI_list.append(MI)

    return MI_list

path = './data/'
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+"hg19_chr1_NCP_ics_anot.cn")
ID_score1 = name_ID_value["data/sp_spd_tests_detail/sp7"]
ID_score2 = name_ID_value["data/sp_spd_tests_detail/sp8"]
ID_seq = name_ID_value['Sequence']

# binning the score
med = np.median(ID_score1.values())
std = np.std(ID_score1.values())/4
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
for ID in ID_score1:
    score1 = ID_score1[ID]
    for i in range(p_num):
        st, ed = p_range[i]
        if score1 >= st and score1 < ed:
            break
    p_IDs[i].append(ID)

ID_bin = {}
for i in range(p_num):
    IDs = p_IDs[i]
    for ID in IDs:
        ID_bin[ID] = i

# kmer x score mutual information analysis
bnum = p_num
knum = 3
step = 1
seq_len = 147

seq_list = []
weight_list = []
bscore_list = []
for ID in ID_seq.keys():
    seq = ID_seq[ID]
    dyad = len(seq)/2
    for i in range(-1,2):
        subseq = seq[dyad+i-seq_len/2:dyad+i+seq_len/2+1]
        assert len(subseq) == seq_len
        seq_list.append(subseq)
        #seq_list.append(rev_comp(subseq))
        weight_list.append(1.0)
        bscore_list.append(ID_bin[ID])

print "Data reading is done"

bscore_prob, kmer_prob_list, kmer_bscore_prob_list = get_kmer_probs (seq_list, weight_list, bscore_list, seq_len, knum, step, bnum, sym=False)
print "prob calculation is done"

MI_list = get_MI_list (bscore_prob, kmer_prob_list, kmer_bscore_prob_list, seq_len, knum, step, bnum)
print "MI calculation is done"

fig = plt.figure()
plt.plot(MI_list)
plt.show()
plt.close()
