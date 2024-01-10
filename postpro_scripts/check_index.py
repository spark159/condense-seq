import copy
import numpy as np
import matplotlib.pyplot as plt

def Hamming_dist(seq1, seq2):
    assert len(seq1) == len(seq2)
    dist = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist +=1
    return dist

def get_pos_fract (seq_list):
    seqlen = len(seq_list[0])
    pos_nt_fract = [{'A':0.0, 'T':0.0, 'C':0.0, 'G':0.0} for i in range(seqlen)] 
    for seq in seq_list:
        assert len(seq) == seqlen
        for i in range(seqlen):
            nt = seq[i]
            nt_fract = pos_nt_fract[i]
            nt_fract[nt] += 1.0/len(seq_list)
    return pos_nt_fract

def read_table(fname):
    sample_idxpair = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        if First:
            First = False
            continue
        sample, index1, index2 = cols[1], cols[5], cols[7]
        assert sample not in sample_idxpair
        sample_idxpair[sample] = (index1, index2)
    return sample_idxpair
sample_idxpair = read_table('2022.08.30 GM_NCP_sp_deepseq(retry)_H1_HP1a_deepseq.csv')

def get_hamming_matrix (sample_idxpair):
    sample_list = sorted(sample_idxpair.keys())
    nanmatrix = np.zeros((len(sample_list), len(sample_list))); nanmatrix[:] = np.nan
    Hamming_matrix1 = copy.deepcopy(nanmatrix)
    Hamming_matrix2 = copy.deepcopy(nanmatrix)
    for i in range(len(sample_list)-1):
        for j in range(i+1, len(sample_list)):
            sample1, sample2 = sample_list[i], sample_list[j]
            idx_pair1, idx_pair2 = sample_idxpair[sample1], sample_idxpair[sample2]
            Hamming_matrix1[i][j] = Hamming_dist(idx_pair1[0], idx_pair2[0])
            Hamming_matrix1[j][i] = copy.deepcopy(Hamming_matrix1[i][j])
            Hamming_matrix2[i][j] = Hamming_dist(idx_pair1[1], idx_pair2[1])
            Hamming_matrix2[j][i] = copy.deepcopy(Hamming_matrix2[i][j])
    return sample_list, Hamming_matrix1, Hamming_matrix2
sample_list, Hamming_matrix1, Hamming_matrix2 = get_hamming_matrix(sample_idxpair)

def plot_hamming_matrix (sample_list, Hamming_matrix, title=None, note=""):
    fig = plt.figure()
    plt.imshow(Hamming_matrix)
    for i in range(len(sample_list)):
        for j in range(len(sample_list)):
            if i == j:
                continue
            dist = int(Hamming_matrix[i][j])
            if dist < 6:
                color = 'white'
            else:
                color = 'black'
            plt.text(j, i, str(dist), ha='center', va='center', weight='bold', size=12, color=color)

    plt.xticks(range(len(sample_list)), sample_list, rotation=45, ha='right', rotation_mode='anchor')
    plt.yticks(range(len(sample_list)), sample_list)
    plt.colorbar()
    plt.title(title)
    plt.savefig("hamming_matrix_" + note + ".png", bbox_inches='tight')
    #plt.show()
    plt.close()

def plot_dual_hamming_matrix (sample_list, Hamming_matrix1, Hamming_matrix2, title=None, note=""):
    dual_matrix = np.zeros((len(sample_list), len(sample_list)))
    dual_matrix[:] = np.nan
    for i in range(len(sample_list)):
        for j in range(len(sample_list)):
            if i == j:
                continue
            elif i < j:
                dist = Hamming_matrix1[i][j]
            else:
                dist = Hamming_matrix2[i][j]
            dual_matrix[i][j] = copy.deepcopy(dist)

    fig = plt.figure()
    plt.imshow(dual_matrix)
    for i in range(len(sample_list)):
        for j in range(len(sample_list)):
            if i == j:
                continue
            dist = int(dual_matrix[i][j])
            if dist < 6:
                color = 'white'
            else:
                color = 'black'
            plt.text(j, i, str(dist), ha='center', va='center', weight='bold', size=12, color=color)
            
    plt.xticks(range(len(sample_list)), sample_list, rotation=45, ha='right', rotation_mode='anchor')
    plt.yticks(range(len(sample_list)), sample_list)
    plt.colorbar()
    plt.title(title)
    plt.savefig("dual_hamming_matrix_" + note + ".png", bbox_inches='tight')
    #plt.show()
    plt.close()


plot_hamming_matrix (sample_list, Hamming_matrix1, title='i7 index Hamming distance', note='i7')
plot_hamming_matrix (sample_list, Hamming_matrix2, title='i5 index Hamming distance', note='i5')
plot_dual_hamming_matrix (sample_list,
                          Hamming_matrix1,
                          Hamming_matrix2,
                          title='i7&i5 index Hamming distance')

