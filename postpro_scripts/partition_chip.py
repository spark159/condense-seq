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

def get_density (ID_CNum, ID_meNum):
    ID_mefrac = {}
    for ID in ID_CNum:
        CNum = ID_CNum[ID]
        if CNum <= 0:
            continue
        meNum = ID_meNum[ID]
        mefrac = float(meNum) / (CNum)
        ID_mefrac[ID] = mefrac
    return ID_mefrac

    

# load files
#path = './data/'
path = ''
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+"H1_NCP_sp_chr1_extended_anot.cn")
ID_score1 = name_ID_value["work/2021_06_07_H1_sp_detail/H1-NCP-sp-8"]
#ID_score2 = name_ID_value["work/2021_06_07_H1_sp_detail/H1-NCP-sp-8"]
ID_seq = name_ID_value['Sequence']
ID_AT = name_ID_value['ATcontent']
ID_meCpGfrac = get_density(name_ID_value['CNumber(CpG)'], name_ID_value['meCNumber(CpG)'])
ID_meCHGfrac = get_density(name_ID_value['CNumber(CHG)'], name_ID_value['meCNumber(CHG)'])
ID_meCHHfrac = get_density(name_ID_value['CNumber(CHH)'], name_ID_value['meCNumber(CHH)'])

# parameters
group_num = 1
vmin, vmax = -1.5, 1.5


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
fig = plt.figure(figsize=(2.4, 2))
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
    plt.text(x, 10000, num_rom[i+1], fontsize=12, va='center', ha='center')
plt.xlim([-3,3])
plt.title("Chromosome1", fontsize=8)
plt.xlabel("Condensability (A.U.)", fontsize=8)
plt.ylabel("Nucleosome Counts", fontsize=8)
plt.gca().tick_params(axis='both', which='major', labelsize=5)
plt.gca().tick_params(axis='both', which='minor', labelsize=5)
#plt.savefig("partition_hist.png")
#plt.savefig('partition.svg', format='svg', bbox_inches='tight')
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

# Reading BS for each partitions
bs_p_sigs = {}
for i in range(p_num):
    IDs = p_IDs[i]
    for ID in IDs:
        for bs, ID_mefrac in zip(['CpG', 'CHG', 'CHH'], [ID_meCpGfrac, ID_meCHGfrac, ID_meCHHfrac]):
            if bs not in bs_p_sigs:
                bs_p_sigs[bs] = [[] for k in range(p_num)]
            try:
                mefrac = ID_mefrac[ID]
            except:
                continue
                #mefrac = 0.0
            bs_p_sigs[bs][i].append(mefrac)


# Reading chip seq data for each partitions
chip_p_sigs = {}
for i in range(p_num):
    IDs = p_IDs[i]
    for ID in IDs:
        for name in name_ID_value:
            if not name.startswith('H'):
                continue
            if name not in chip_p_sigs:
                chip_p_sigs[name] = [[] for k in range(p_num)]
            ID_value = name_ID_value[name]
            chip_p_sigs[name][i].append(ID_value[ID])

print "Epigenetic marks reading is done"

    
# calculate z-scores
names = sorted(bs_p_sigs.keys()) + sorted(chip_p_sigs.keys())
group_size = len(names)/group_num
img_list = [[] for i in range(group_num)]
ylabels_list = [[] for i in range(group_num)]

for i in range(len(names)):
    name = names[i]
    if i < len(bs_p_sigs):
        name_p_sigs = bs_p_sigs
    else:
        name_p_sigs = chip_p_sigs
    g_idx = i / group_size
    img_list[g_idx].append([np.mean(zscores) for zscores in get_zscore(name_p_sigs[name])])
    if name.startswith('C'):
        name += ' density'
    ylabels_list[g_idx].append(name)


# plot the z-scores over partitions
#fig, axes = plt.subplots(nrows=1, ncols=group_num, figsize=(5.8,4))
fig, axes = plt.subplots(nrows=1, ncols=group_num, figsize=(2.5,8))
for i in range(group_num):
    img = img_list[i]
    ylabels = ylabels_list[i]
    if group_num > 1:
        ax = axes[i]
    else:
        ax = plt.gca()
    im = ax.imshow(img, vmin=vmin, vmax=vmax, cmap='coolwarm', aspect='auto')
    ax.set_xticks(range(p_num))
    ax.set_xticklabels(range(1, p_num+1), fontsize=8)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(img)))    
    ax.set_yticklabels(ylabels, horizontalalignment='right', fontname='monospace', fontsize=9)
plt.tight_layout()
plt.savefig('partition_chip.svg', format='svg', bbox_inches='tight')
#plt.show()
plt.close()


# make color bar
fig = plt.figure(figsize=(1.2,1))
plt.subplot(1,2,1)
cbar = plt.colorbar(im, cax=plt.gca(), ticks=[vmin, vmax])
cbar.ax.set_yticklabels([vmin, vmax], fontsize=8)
cbar.ax.set_ylabel('Z-score', rotation=-90, va="bottom", fontsize=8)
plt.tight_layout()
plt.savefig('partition_chip_cbar.svg', format='svg', bbox_inches='tight')
plt.close()
