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

def get_kmer_probs (seq_list, weight_list, seq_len, knum, step, sym=False):
    pos_num = (seq_len - knum + 1) / step
    all_kmers = all_path(knum, 'ATCG')

    # initialize
    kmer_prob = {kmer:0.0 for kmer in all_kmers}
    kmer_pair_prob = {}
    for kmer1 in all_kmers:
        for kmer2 in all_kmers:
            kmer_pair_prob[(kmer1, kmer2)] = 0.0

    kmer_prob_list = [copy.deepcopy(kmer_prob) for i in range(pos_num)]
    kmer_pair_prob_matrix = {}
    for i in range(pos_num-1):
        kmer_pair_prob_matrix[i] = {}
        for j in range(i+1, pos_num):
            kmer_pair_prob_matrix[i][j] = copy.deepcopy(kmer_pair_prob)

    # counting
    for k in range(len(seq_list)):
        seq = seq_list[k].upper()
        if sym:
            rev_seq = rev_comp(seq)
        weight = weight_list[k]
        for i in range(pos_num):
            kmer1 = seq[i*step:i*step+knum]
            try:
                kmer_prob_list[i][kmer1] += weight
            except:
                pass
            if sym:
                kmer1_rev = rev_seq[i*step:i*step+knum]
                try:
                    kmer_prob_list[i][kmer1_rev] += weight
                except:
                    pass
            for j in range(i+1, pos_num):
                kmer2 = seq[j*step:j*step+knum]
                try:
                    kmer_pair_prob_matrix[i][j][(kmer1,kmer2)] += weight
                except:
                    pass
                if sym:
                    kmer2_rev = rev_seq[j*step:j*step+knum]
                    try:
                        kmer_pair_prob_matrix[i][j][(kmer1_rev,kmer2_rev)] += weight
                    except:
                        pass

    # normalize
    for kmer_prob in kmer_prob_list:
        total = float(sum(kmer_prob.values()))
        for kmer in kmer_prob:
            kmer_prob[kmer] = kmer_prob[kmer] / total

    for i in range(pos_num-1):
        for j in range(i+1, pos_num):
            kmer_pair_prob = kmer_pair_prob_matrix[i][j]
            total = float(sum(kmer_pair_prob.values()))
            for pair in kmer_pair_prob:
                kmer_pair_prob[pair] = kmer_pair_prob[pair] / total

    return kmer_prob_list, kmer_pair_prob_matrix

def get_kmer_MI_matrix (kmer_prob_list, kmer_pair_prob_matrix, seq_len, knum, step):
    def get_MI (kmer_pair_prob, kmer_prob1, kmer_prob2):
        MI = 0.0
        for kmer1, kmer2 in kmer_pair_prob:
            joint_prob = kmer_pair_prob[(kmer1, kmer2)]
            prob1, prob2 = kmer_prob1[kmer1], kmer_prob2[kmer2]
            if joint_prob <= 0 or prob1 <=0 or prob2 <= 0:
                continue
            MI += joint_prob*np.log2(joint_prob/(prob1*prob2))
        return MI
    
    pos_num = (seq_len - knum + 1) / step

    kmer_pair_MI_matrix = {}
    for i in range(pos_num-1):
        kmer_pair_MI_matrix[i] = {}
        for j in range(i+1, pos_num):
            kmer_prob1 = kmer_prob_list[i]
            kmer_prob2 = kmer_prob_list[j]
            kmer_pair_prob = kmer_pair_prob_matrix[i][j]
            MI = get_MI (kmer_pair_prob, kmer_prob1, kmer_prob2)
            kmer_pair_MI_matrix[i][j] = MI

    return kmer_pair_MI_matrix

def visualize_matrix (matrix):
    N = len(matrix)
    img = np.zeros((N,N))
    img[:] = np.nan
    for i in range(N-1):
        for j in range(i+1, N):
            img[i][j] = matrix[i][j]

    fig = plt.figure()
    plt.imshow(img)
    plt.colorbar()
    plt.show()
    plt.close()

    return img
            

path = './data/'
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+"hg19_chr1_NCP_ics_anot.cn")
ID_score1 = name_ID_value["data/sp_spd_tests_detail/sp7"]
ID_score2 = name_ID_value["data/sp_spd_tests_detail/sp8"]
ID_seq = name_ID_value['Sequence']


# pairwise kmer mutual information analysis
knum = 2
step = 2
seq_len = 147


seq_list = []
weight_list1, weight_list2 = [], [] 
for ID in ID_seq.keys():
    seq = ID_seq[ID]
    dyad = len(seq)/2
    for i in range(-1,2):
    #for i in [0]:
        subseq = seq[dyad+i-seq_len/2:dyad+i+seq_len/2+1]
        assert len(subseq) == seq_len
        seq_list.append(subseq)
        #seq_list.append(rev_comp(subseq))
        weight_list1.append(1.0)
        weight = np.exp(-ID_score1[ID])
        weight_list2.append(weight)

"""
def AG_freq (NCP_seq_list):
    nucleosome_dna_len = len(NCP_seq_list[0])
    Afreq=np.zeros(nucleosome_dna_len - 1); Gfreq=np.zeros(nucleosome_dna_len - 1)
    for seq in NCP_seq_list:
        for i in range(len(seq)-1):
            dint = seq[i:i+2]
            if dint in ['AA','AT','TA','TT']:
                Afreq[i] += 1.0
            elif dint in ['CC','CG','GC','GG']:
                Gfreq[i] += 1.0
    return Afreq / len(NCP_seq_list), Gfreq / len(NCP_seq_list)

Afreq1, Gfreq1 = AG_freq(seq_list)

fig, ax1 = plt.subplots()
ax1.plot(Afreq1, 'r')
ax1.set_ylabel('AA/AT/TA/TT freqeuncy', color='r')
ax1.tick_params('y', colors='r')
ax2 = ax1.twinx()
ax2.plot(Gfreq1, 'b')
ax2.set_ylabel('CC/CG/GC/GG freqeuncy', color='b')
ax2.tick_params('y', colors='b')
mid = 147/2
line_list = [mid - i*20 for i in range(4)] + [mid + i*20 for i in range(1, 4)]
xlabel_list = ['SHL' + str(-2*i) for i in range(4)] + ['SHL' + str(2*i) for i in range(1, 4)]

for line in line_list:
    plt.axvline(x=line, color='k', linestyle='--', alpha=0.5)

ax1.set_xticks(line_list)
ax1.set_xticklabels(xlabel_list)
plt.title("Dinucleotide frequency")
plt.savefig("DinFreq.png", bbox_inches='tight')
plt.show()
plt.close()

sys.exit(1)
"""

"""
seq_list = []
weight_list1, weight_list2 = [], []
for i in range(1000):
    seq = ""
    for j in range(seq_len):
        seq += random.choice(list('ATCG'))
    seq_list.append(seq)
    weight_list1.append(1.0)
    weight_list2.append(1.0)
"""

"""
path = "/home/spark159/script/slide-seq/"
with open(path+"plusonelib_new:corrected_0.pickle", "rb") as f:
    key_slider = pickle.load(f)

seq_list = []
weight_list1, weight_list2 = [], []
for key in key_slider:
    seq = key_slider[key].seq
    seq = seq[40:len(seq)-40]
    assert len(seq) == seq_len
    seq_list.append(seq)
    weight_list1.append(1.0)
    weight_list2.append(1.0)
"""

"""
seq_list = []
weight_list1, weight_list2 = [], []
for i in range(1000):
    seq = ""
    for j in range(10):
        seq += 'AA'
        for u in range(5):
            seq += random.choice(list('ATCG'))
    seq_len = len(seq)
    seq_list.append(seq)
    weight_list1.append(1.0)
    weight_list2.append(1.0)
"""

print "Data reading is done"

kmer_prob_list1, kmer_pair_prob_matrix1 = get_kmer_probs (seq_list, weight_list1, seq_len, knum, step, sym=False)
print "prob calculation is done"

MI_matrix1 = get_kmer_MI_matrix (kmer_prob_list1, kmer_pair_prob_matrix1, seq_len, knum, step)
print "MI matrix is done"

kmer_prob_list2, kmer_pair_prob_matrix2 = get_kmer_probs (seq_list, weight_list2, seq_len, knum, step, sym=False)
print "prob calculation is done"

MI_matrix2 = get_kmer_MI_matrix (kmer_prob_list2, kmer_pair_prob_matrix2, seq_len, knum, step)
print "MI matrix is done"

# plot MI matrix
N = len(MI_matrix1)

img = np.zeros((N, N))
img[:] = np.nan
for i in range(N-1):
    for j in range(i+1, N):
        img[i][j] = MI_matrix1[i][j]
        img[j][i] = MI_matrix2[i][j]
        
fig = plt.figure()
plt.imshow(img)
plt.colorbar()
plt.title('Pairwise ' + str(knum) + '-mer MI matrix')
plt.show()
plt.close()
