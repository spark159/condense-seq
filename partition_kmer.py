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


#path = './data/'
path = ''
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+"H1_NCP_sp_chr1_anot.cn")
ID_score1 = name_ID_value["work/2021_06_07_H1_sp_detail/H1-NCP-sp-8"]
#ID_score2 = name_ID_value["work/2021_06_07_H1_sp_detail/H1-NCP-sp-8"]
ID_seq = name_ID_value['Sequence']
ID_AT = name_ID_value['ATcontent']
ID_CG = name_ID_value['CNumber(CpG)']
ID_me = name_ID_value['meCNumber(CpG)']

# Partition by score
med = np.median(ID_score1.values())
std = np.std(ID_score1.values())
lines = [med-0.5*std-i*std for i in range(3)] + [med+0.5*std+i*std for i in range(3)]
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
fig = plt.figure()
plt.hist(ID_score1.values(), bins=1000)
for line in lines:
    plt.axvline(x=line, color='k', linestyle='--')
num_rom = {1:'I', 2:'II', 3:'III', 4:'IV', 5:'V', 6:'VI', 7:'VII'}
for i in range(p_num):
    st, ed = p_range[i]
    if i == 0:
        x = np.mean([-3, ed])
    elif i == p_num - 1:
        x = np.mean([st, 3])
    else:
        x = np.mean([st, ed])
    plt.text(x, 10000, num_rom[i+1], fontsize=20, va='center', ha='center')
plt.xlim([-3,3])
plt.title("Chromosome1")
plt.xlabel("Condensability (A.U.)")
plt.ylabel("Nucleosome Counts")
#plt.savefig("partition_hist.png")
#plt.show()
plt.close()

p_IDs = [[] for i in range(p_num)]
for ID in ID_score1:
    score1 = ID_score1[ID]
    for i in range(p_num):
        st, ed = p_range[i]
        if score1 >= st and score1 < ed:
            break
    p_IDs[i].append(ID)

# parameters
klen = 4
clip = 10
NCPlen = 147

# Count kmers
all_kmers = [kmer for GC, kmer in sorted([(GC_content(kmer), kmer) for kmer in all_path(klen, 'ATCG')])]
kmer_p_sigs = {}
for i in range(p_num):
    IDs = p_IDs[i]
    for ID in IDs:
        seq = ID_seq[ID]
        middle = len(seq)/2
        seq = seq[middle-NCPlen/2:middle+NCPlen/2+1]
        assert len(seq) == NCPlen
        seq = seq[clip:len(seq)-clip].upper()
        kmer_count = get_kmer_count(seq, klen=klen, both=True) # survey both strands
        for kmer in all_kmers:
            try:
                count = kmer_count[kmer]
            except:
                count = 0
            if kmer not in kmer_p_sigs:
                kmer_p_sigs[kmer] = [[] for k in range(p_num)]
            kmer_p_sigs[kmer][i].append(count)

print "Sequence feature reading is done"

    
#sys.exit(1)

# calculate mean z-scores over partitions
group_num = klen
group_size = len(all_kmers) / group_num
img_list = [[] for i in range(group_num)]
ylabels_list = [[] for i in range(group_num)]
for i in range(len(all_kmers)):
    kmer = all_kmers[i]
    g_idx = i / group_size
    img_list[g_idx].append([np.mean(zscores) for zscores in get_zscore(kmer_p_sigs[kmer])])
    ylabels_list[g_idx].append(kmer)


# plot the z-scores over partitions
fig, axes = plt.subplots(nrows=1, ncols=group_num, figsize=(5.8,8.4))
for i in range(group_num):
    img = img_list[i]
    ylabels = ylabels_list[i]
    axes[i].imshow(img, vmin=-1.5, vmax=1.5, cmap='coolwarm', aspect='auto')
    axes[i].set_xticks(range(p_num))
    axes[i].set_xticklabels(range(1, p_num+1), fontsize=8)
    axes[i].xaxis.tick_top()
    axes[i].set_yticks(range(len(img)))    
    axes[i].set_yticklabels(ylabels, horizontalalignment='right', fontname='monospace', fontsize=9)
plt.tight_layout()
plt.savefig('partition_kmer.svg', format='svg', bbox_inches='tight')
#plt.show()
plt.close()
