#import load_file
#import graphics
from scipy import stats
import sys
import copy
#import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
#import Interval_dict
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.cm as cm

def GC_content(seq, percent=False):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    if percent:
        return (num/float(len(seq)))*100
    return (num/float(len(seq)))

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

def all_path (N, states='ATCG'):
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output

def poly_score (seq, nts='AT', pos=False):
    num = []
    num_pos = {}
    i = 0
    while i < len(seq):
        if seq[i] in nts:
            nt = seq[i]
            count = 1
            j = i + 1
            while j < len(seq):
                if seq[j] != nt:
                    break
                count +=1
                j +=1
            num.append(count)
            if count not in num_pos:
                num_pos[count] = []
            num_pos[count].append(i)
            i = j
        else:
            i +=1
    if pos:
        return num_pos
    if len(num) == 0:
        return 0
    return max(num)

def get_dincount(seq, din=None, both=False):
    if din:
        count = 0
        for i in range(len(seq)-1):
            if seq[i:i+2].upper() == din.upper():
                count +=1
        return count
    din_count = {}
    seq = seq.upper()
    for i in range(len(seq)-1):
        din = seq[i:i+2]
        if 'N' in din:
            continue
        if din not in din_count:
            din_count[din] = 0
        din_count[din] += 1
        if both:
            din = rev_comp(din)
            if din not in din_count:
                din_count[din] = 0
            din_count[din] += 1
    return din_count

def get_kmer_count(seq, klen, both=False):
    kmer_count = {}
    seq = seq.upper()
    seq_list = seq.split('N')
    for seq in seq_list:
        if len(seq) < klen:
            continue
        if both:
            rev_seq = rev_comp(seq)
        for i in range(len(seq)-klen-1):
            kmer = seq[i:i+klen]
            if kmer not in kmer_count:
                kmer_count[kmer] = 0
            kmer_count[kmer] += 1
            if both:
                kmer = rev_seq[i:i+klen]
                if kmer not in kmer_count:
                    kmer_count[kmer] = 0
                kmer_count[kmer] += 1
    return kmer_count

def get_zscore(p_data):
    linear_data = []
    p_num = []
    for data in p_data:
        p_num.append(len(data))
        for value in data:
            linear_data.append(value)
    linear_data = stats.zscore(linear_data)
    new_data = []
    st = 0
    for num in p_num:
        new_data.append(linear_data[st:st+num])
        st += num
    assert st == len(linear_data)
    return new_data

def stat_Markov(seq_list, NCPlen, order):
    ntdic = {}
    for nt in all_path(order+1, 'ATCG'):
        ntdic[nt] = 0.0

    sample_num = len(seq_list)

    freq = [ copy.deepcopy(ntdic) for i in range(NCPlen - order) ]
    for seq in seq_list:
        for i in range(len(seq) - order):
            nt = seq[i:i+1+order]
            freq[i][nt] += 1.0 / sample_num

    mean, std = [], []
    for ntdic in freq:
        mean.append(np.mean(ntdic.values()))
        std.append(np.std(ntdic.values()))

    stdz_freq = []
    for i in range(len(freq)):
        ntdic = freq[i]
        temp = {}
        for nt, value in ntdic.items():
            temp[nt] = (value - mean[i]) / std[i]
        stdz_freq.append(temp)

    return freq, sample_num, mean, std, stdz_freq

def stat_Kmer(seq_list, NCPlen, knum, bnum):
    seqlen = NCPlen / bnum
    assert seqlen >= knum

    extra = NCPlen % bnum
    boundoff = extra / 2
    centeroff = extra % 2

    ntdic = {}
    for nt in all_path(knum, 'ATCG'):
        ntdic[nt] = 0.0
    freq = [ copy.deepcopy(ntdic) for i in range(bnum) ]

    sample_num = len(seq_list)*(seqlen-knum+1)

    for seq in seq_list:
        for k in range(bnum):
            if k < bnum/2:
                st = boundoff + k*seqlen
            if k >= bnum/2:
                st = boundoff + centeroff + k*seqlen
            bseq = seq[st:st+seqlen]                
            for i in range(seqlen - knum + 1):
                nt = bseq[i:i+knum]
                freq[k][nt] += 1.0 / sample_num

    mean, std = [], []
    for ntdic in freq:
        mean.append(np.mean(ntdic.values()))
        std.append(np.std(ntdic.values()))

    stdz_freq = []
    for i in range(len(freq)):
        ntdic = freq[i]
        temp = {}
        for nt, value in ntdic.items():
            temp[nt] = (value - mean[i]) / std[i]
        stdz_freq.append(temp)

    return freq, sample_num, mean, std, stdz_freq


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
tname_ID_score = read_score("YWlib_sp_score.txt")
ID_score1 = tname_ID_score[3]
ID_score2 = tname_ID_score[5]

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
ID_seq = read_ref('Blib.ref')


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

p_num = 5
p_IDs = quantile(ID_score1, p_num)

# parameters
#poly_st, poly_ed = 1, 7 # homopolymer length range
clip = 0
NCPlen = 101


# clip off sequence
for ID in ID_seq:
    seq = ID_seq[ID]
    middle = len(seq)/2
    seq = seq[middle-NCPlen/2:middle+NCPlen/2+1]
    assert len(seq) == NCPlen
    ID_seq[ID] = seq[clip:len(seq)-clip].upper()


def get_loc_lens (seq_list, nts):
    loc_lens = {}
    for seq in seq_list:
        len_pos = poly_score(seq, nts=nts, pos=True)
        for length, pos_list in len_pos.items():
            for pos in pos_list:
                for i in range(pos, pos+length):
                    if i not in loc_lens:
                        loc_lens[i] = []
                    loc_lens[i].append(length)
    return loc_lens

nts_list = ['AT','GC']

# Check poly(dA:dT) and poly(dG:dC) tract length and locations (all data)
control_nts_loc_meanlen = {}
for nts in nts_list:
    loc_lens = get_loc_lens(ID_seq.values(), nts=nts)
    loc_meanlen = {loc:np.mean(loc_lens[loc]) for loc in loc_lens}
    control_nts_loc_meanlen[nts] = loc_meanlen



# Check poly(dA:dT) and poly(dG:dC) tract length and locations (each quantile)
nts_p_loc_meanlen = {}
for i in range(p_num):
    IDs = p_IDs[i]
    seq_list = [ID_seq[ID] for ID in IDs]

    for nts in nts_list:
        loc_lens = get_loc_lens (seq_list, nts=nts)
        loc_meanlen = {loc:np.mean(loc_lens[loc]) for loc in loc_lens}
        if nts not in nts_p_loc_meanlen:
            nts_p_loc_meanlen[nts] = []
        nts_p_loc_meanlen[nts].append(loc_meanlen)


    
# plot the fold change
cmap_list = ['OrRd', 'GnBu']
color_list = np.linspace(0.2, 1, num=p_num)
for i in range(len(nts_list)):
    nts = nts_list[i]
    control_loc_meanlen = control_nts_loc_meanlen[nts]
    assert len(control_loc_meanlen) == NCPlen - 2*clip
    p_loc_meanlen = nts_p_loc_meanlen[nts]
    cmap = cm.get_cmap(cmap_list[i])
    fig = plt.figure(figsize=(2.6, 2))
    for p in range(p_num):
        loc_meanlen = p_loc_meanlen[p]
        assert len(loc_meanlen) == NCPlen - 2*clip
        fold_list = []
        for j in range(NCPlen - 2*clip):
            fold = float(loc_meanlen[j]) / control_loc_meanlen[j]
            fold_list.append(fold)
        
        plt.plot(fold_list, color=cmap(color_list[p]), lw=1.5, alpha=1)
        #plt.ylabel("Mean tract length (bp)", fontsize=8)
        plt.ylabel("Fold enrichment", fontsize=8)
        plt.xlabel("Position (bp)", fontsize=8)
        line_list = [101/2 - k*20 for k in range(3)] + [101/2 + k*20 for k in range(1, 3)]
        for line in line_list:
            plt.axvline(x=line, color='k', linestyle='--', alpha=0.25)
        plt.xticks([101/2 + 10*k for k in range(-5, 6, 2)], [str(10*k) for k in range(-5, 6, 2)], fontsize=5)
        plt.gca().tick_params('x', labelsize=6)
        plt.gca().tick_params('y', labelsize=6)

    plt.title('Poly(d%s:d%s)' % tuple(nts), fontsize=10)
    plt.savefig('polyloc_' + nts + '.svg', format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()
