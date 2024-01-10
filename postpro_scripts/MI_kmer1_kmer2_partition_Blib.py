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

def get_kmer_probs (seq_list, weight_list, seq_len, knum, step, sym=False):
    pos_num = (seq_len - knum + 1) / step
    all_kmers = all_path(knum, 'ATCG')

    # initialize
    kmer_prob = {kmer:0.0 for kmer in all_kmers}
    pair_prob = {}
    for kmer1 in all_kmers:
        for kmer2 in all_kmers:
            pair_prob[(kmer1, kmer2)] = 0.0

    pos_kmer_prob = [copy.deepcopy(kmer_prob) for i in range(pos_num)]
    pos1_pos2_pair_prob = {}
    for i in range(pos_num-1):
        pos1_pos2_pair_prob[i] = {}
        for j in range(i+1, pos_num):
            pos1_pos2_pair_prob[i][j] = copy.deepcopy(pair_prob)

    # counting
    for k in range(len(seq_list)):
        seq = seq_list[k].upper()
        if sym:
            rev_seq = rev_comp(seq)
        weight = weight_list[k]
        for i in range(pos_num):
            kmer1 = seq[i*step:i*step+knum]
            try:
                pos_kmer_prob[i][kmer1] += weight
            except:
                pass
            if sym:
                kmer1_rev = rev_seq[i*step:i*step+knum]
                try:
                    pos_kmer_prob[i][kmer1_rev] += weight
                except:
                    pass
            for j in range(i+1, pos_num):
                kmer2 = seq[j*step:j*step+knum]
                try:
                    pos1_pos2_pair_prob[i][j][(kmer1,kmer2)] += weight
                except:
                    pass
                if sym:
                    kmer2_rev = rev_seq[j*step:j*step+knum]
                    try:
                        pos1_pos2_pair_prob[i][j][(kmer1_rev,kmer2_rev)] += weight
                    except:
                        pass

    # normalize
    for kmer_prob in pos_kmer_prob:
        total = float(sum(kmer_prob.values()))
        for kmer in kmer_prob:
            kmer_prob[kmer] = kmer_prob[kmer] / total

    for i in range(pos_num-1):
        for j in range(i+1, pos_num):
            pair_prob = pos1_pos2_pair_prob[i][j]
            total = float(sum(pair_prob.values()))
            for pair in pair_prob:
                pair_prob[pair] = pair_prob[pair] / total

    return pos_kmer_prob, pos1_pos2_pair_prob

def get_kmer_MI_matrix (pos_kmer_prob, pos1_pos2_pair_prob, seq_len, knum, step):
    def get_MI (pair_prob, kmer_prob1, kmer_prob2):
        MI = 0.0
        for kmer1, kmer2 in pair_prob:
            joint_prob = pair_prob[(kmer1, kmer2)]
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
            kmer_prob1 = pos_kmer_prob[i]
            kmer_prob2 = pos_kmer_prob[j]
            pair_prob = pos1_pos2_pair_prob[i][j]
            MI = get_MI (pair_prob, kmer_prob1, kmer_prob2)
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



# partition by score
def quantile (ID_score, num):
    def value_cmp(a, b):
        if a[1] <= b[1]:
            return -1
        else:
            return 1

    assert len(ID_score) >= num
    IDscore = [(ID, ID_score[ID]) for ID in ID_score.keys()]
    IDscore = sorted(IDscore, cmp=value_cmp)

    size, remain = len(IDscore) / num, len(IDscore) % num
    size_list = [size]*num

    for i in range(remain):
        size_list[i] +=1

    assert sum(size_list) == len(IDscore)

    q_IDs = []
    for i in range(num):
        size = size_list[i]
        st = i*size
        ed = st + size
        q_IDs.append([IDscore[j][0] for j in range(st,ed)])
    return q_IDs

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

seq_len = 101
clip = 0 # clip off the both ends of sequence 
knum = 3
step = knum
p_num = 5
note = '_YW'

# slice out the sequence
for ID in ID_seq:
    seq = ID_seq[ID]
    middle = len(seq)/2
    seq = seq[middle-seq_len/2:middle+seq_len/2+1]
    ID_seq[ID] = seq[clip:len(seq)-clip]
    #ID_seq[ID] = ''.join([random.choice('ATCG') for i in range(seq_len)])

# MI for all data
repeat = 10
sample_size = len(ID_seq) / p_num
MI_matrix0 = {}
for i in range(repeat):
    print i
    seq_list = random.sample(ID_seq.values(), sample_size) 
    weight_list = [1]*len(seq_list)

    pos_kmer_prob, pos1_pos2_pair_prob = get_kmer_probs (seq_list, weight_list, seq_len, knum, step, sym=True)

    MI_matrix = get_kmer_MI_matrix (pos_kmer_prob, pos1_pos2_pair_prob, seq_len, knum, step)

    for pos1 in MI_matrix:
        for pos2 in MI_matrix[pos1]:
            MI = MI_matrix[pos1][pos2]
            if pos1 not in MI_matrix0:
                MI_matrix0[pos1] = {}
            if pos2 not in MI_matrix0[pos1]:
                MI_matrix0[pos1][pos2] = 0.0
            MI_matrix0[pos1][pos2] += MI / float(repeat)

print "MI matrix0 is done"
print

# plot MI matrix
N = len(MI_matrix0)

img = np.zeros((N, N))
img[:] = np.nan
for i in range(N-1):
    for j in range(i+1, N):
        img[i][j] = MI_matrix0[i][j]
        img[j][i] = MI_matrix0[i][j]
        
fig = plt.figure()
plt.imshow(img, cmap='jet')
plt.colorbar()
plt.title('Pairwise ' + str(knum) + '-mer MI matrix (control)')
plt.show()
plt.close()


# quantile the group
p_IDs = quantile(ID_score, p_num)

p_MImatrix = []
for p in range(p_num):
    IDs = p_IDs[p]
    
    seq_list = [ID_seq[ID] for ID in IDs]
    weight_list = [1]*len(IDs)

    pos_kmer_prob, pos1_pos2_pair_prob = get_kmer_probs (seq_list, weight_list, seq_len, knum, step, sym=True)
    print "prob calculation is done"

    MI_matrix = get_kmer_MI_matrix (pos_kmer_prob, pos1_pos2_pair_prob, seq_len, knum, step)
    print "MI matrix is done"

    p_MImatrix.append(copy.deepcopy(MI_matrix))


# plot MI fold enrichment
for p in range(p_num):
    MI_matrix = p_MImatrix[p]
    img = np.zeros((N, N))
    img[:] = np.nan
    for i in range(N-1):
        for j in range(i+1, N):
            img[i][j] = MI_matrix[i][j] / MI_matrix0[i][j]
            img[j][i] = MI_matrix[i][j] / MI_matrix0[i][j]

    fig = plt.figure()
    plt.imshow(img, cmap='seismic', vmin=0.2, vmax=1.8)
    plt.colorbar()
    plt.title('Pairwise ' + str(knum) + '-mer MI matrix (enrich) ' + str(p+1))
    plt.show()
    plt.close()


"""
# pairwise kmer mutual information analysis
seq_list = []
weight_list1, weight_list2 = [], [] 
for ID in ID_seq:
    seq = ID_seq[ID]
    middle = len(seq)/2
    seq = seq[middle-seq_len/2:middle+seq_len/2+1]
    #ID_seq[ID] = seq[clip:len(seq)-clip]
    seq_list.append(seq)
    weight_list1.append(1.0)
    weight_list2.append(np.exp(-ID_score[ID]))

print "Data reading is done"

pos_kmer_prob1, pos1_pos2_pair_prob1 = get_kmer_probs (seq_list, weight_list1, seq_len, knum, step, sym=True)
print "prob calculation is done"

MI_matrix1 = get_kmer_MI_matrix (pos_kmer_prob1, pos1_pos2_pair_prob1, seq_len, knum, step)
print "MI matrix is done"

pos_kmer_prob2, pos1_pos2_pair_prob2 = get_kmer_probs (seq_list, weight_list2, seq_len, knum, step, sym=True)
print "prob calculation is done"

MI_matrix2 = get_kmer_MI_matrix (pos_kmer_prob2, pos1_pos2_pair_prob2, seq_len, knum, step)
print "MI matrix is done"

# plot MI matrix
N = len(MI_matrix1)

img = np.zeros((N, N))
img[:] = np.nan
for i in range(N-1):
    for j in range(i+1, N):
        img[i][j] = MI_matrix1[i][j]
        img[j][i] = MI_matrix1[i][j]
        
fig = plt.figure()
plt.imshow(img)
plt.colorbar()
plt.title('Pairwise ' + str(knum) + '-mer MI matrix')
plt.show()
plt.close()

# plot MI fold enrichment
img = np.zeros((N, N))
img[:] = np.nan
for i in range(N-1):
    for j in range(i+1, N):
        img[i][j] = MI_matrix2[i][j] / MI_matrix1[i][j]
        img[j][i] = MI_matrix2[i][j] / MI_matrix1[i][j]
        
fig = plt.figure()
plt.imshow(img)
plt.colorbar()
plt.title('Pairwise ' + str(knum) + '-mer MI matrix (enrichment)')
plt.show()
plt.close()


"""
