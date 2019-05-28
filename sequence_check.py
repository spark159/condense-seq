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
import LinModel

def all_path (N, states='ATCG'):
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output

def get_dincount(seq, din=None):
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
    return din_count
        
def Amer_len(seq, pos=False):
    num = []
    num_pos = {}
    i = 0
    while i < len(seq):
        if seq[i] in 'AT':
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

ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("data/hg19_chr1_171_icdseq_anot.cn")
print "reading is done"

ID_score1 = name_ID_value['work/condense_seq/sp9_hg19_chr1']
ID_score2 = name_ID_value['work/condense_seq/sp10_hg19_chr1']
ID_AT = name_ID_value['ATcontent']
ID_CpG = name_ID_value['CpGNumber']
ID_me = name_ID_value['meGCNumber']
ID_seq = name_ID_value['Sequence']

new_ID_score2 = statis.neutralize_score_by_target(ID_score2, ID_AT)

ID_din_count = {}
for ID in ID_seq:
    seq = ID_seq[ID]
    din_count = get_dincount(seq)
    ID_din_count[ID] = din_count

    
# dinucleotide counting
din_count_scores = {}
for ID in ID_din_count:
    score = new_ID_score2[ID]
    #score = ID_score2[ID]
    din_count = ID_din_count[ID]
    for din in din_count:
        count = din_count[din]
        if din not in din_count_scores:
            din_count_scores[din] = {}
        if count not in din_count_scores[din]:
            din_count_scores[din][count] = []
        din_count_scores[din][count].append(score)

for din in din_count_scores:
    count_scores = din_count_scores[din]
    X = sorted(count_scores.keys())
    Y, Z = [], []
    for count in X:
        mean = np.mean(count_scores[count])
        err = np.std(count_scores[count]) / np.sqrt(len(count_scores[count]))
        Y.append(mean)
        Z.append(err)
    fig = plt.figure()
    plt.errorbar(X, Y, yerr=Z, fmt='.', color='skyblue')
    plt.plot(X, Y, '.', color='black', zorder=10)
    plt.title(din)
    plt.xlabel("Dinucleotide number")
    plt.ylabel("Condensability")
    plt.show()
    plt.close()

    
# Partition by condensability
group_IDs = statis.quantile (new_ID_score2, 1000)

din_group_counts = {}
for i in range(len(group_IDs)):
    IDs = group_IDs[i]
    for ID in IDs:
        din_count = ID_din_count[ID]
        for din in din_count:
            if din not in din_group_counts:
                din_group_counts[din] = [ [] for i in range(len(group_IDs)) ]
            din_group_counts[din][i].append(din_count[din])

din_group_mean = {}
din_group_err = {}
for din in din_group_counts:
    group_counts = din_group_counts[din]
    for i in range(len(group_counts)):
        counts = group_counts[i]
        if not counts:
            continue
        mean = np.mean(counts)
        err = np.std(counts) / np.sqrt(len(counts))
        if din not in din_group_mean:
            din_group_mean[din] = [ np.nan for i in range(len(group_IDs)) ]
        din_group_mean[din][i] = mean
        if din not in din_group_err:
            din_group_err[din] = [ np.nan for i in range(len(group_IDs)) ]
        din_group_err[din][i] = err

X = range(1, len(group_IDs)+1)
for din in din_group_mean:
    Y = din_group_mean[din]
    Z = din_group_err[din]
    fig = plt.figure()
    plt.errorbar(X, Y, yerr=Z, fmt='.', color='skyblue')
    plt.plot(X, Y, '.', color='black', zorder=10)
    plt.title(din)
    plt.xlabel("Partition by Condensability")
    plt.ylabel("Dinucleotide number")
    plt.show()
    plt.close()


"""
# Partition by din count

for ID in ID_AT:
    ID_AT[ID] = ID_AT[ID]*100

ID_Alen = {}
for ID in ID_seq:
    seq = ID_seq[ID]
    num_pos = Amer_len(seq, pos=True)
    score = 0.0
    for num, pos in num_pos.items():
        if num < 3:
            continue
        score += num*len(pos)
    ID_Alen[ID] = score

new_ID_score2 = statis.neutralize_score_by_target(ID_score2, ID_AT)

frac = [(4**i) for i in range(1,11)]
#frac = [(2**i) for i in range(1,11)]
frac = frac[::-1]

all_din = all_path(2)

for din in all_din:
    ID_din = {}
    for ID in ID_seq:
        seq = ID_seq[ID]
        count = get_dincount(seq, din=din)
        ID_din[ID] = count

    graphics.PartitionMeanplot(ID_AT, ID_score2, ID_din, frac, note="sp10_" + din)
    graphics.PartitionScatterplot(ID_AT, ID_score2, ID_din, frac, note="sp10_" + din)
    graphics.PartitionBoxplot(new_ID_score2, ID_din, frac, note="sp10_" + din)


graphics.PartitionMeanplot(ID_AT, ID_score2, ID_CpG, frac, note="sp10_CpG")
graphics.PartitionScatterplot(ID_AT, ID_score2, ID_CpG, frac, note="sp10_CpG")
graphics.PartitionBoxplot(new_ID_score2, ID_CpG, frac, note="sp10_CpG")

graphics.PartitionMeanplot(ID_AT, ID_score2, ID_me, frac, note="sp10_me")
graphics.PartitionScatterplot(ID_AT, ID_score2, ID_me, frac, note="sp10_me")
graphics.PartitionBoxplot(new_ID_score2, ID_me, frac, note="sp10_me")

graphics.PartitionMeanplot(ID_AT, ID_score2, ID_Alen, frac, note="sp10_Alen")
graphics.PartitionScatterplot(ID_AT, ID_score2, ID_Alen, frac, note="sp10_Alen")
graphics.PartitionBoxplot(new_ID_score2, ID_Alen, frac, note="sp10_Alen")


X, Y, Z = [], [], []
for ID in ID_score2:
    score2 = ID_score2[ID]
    AT = ID_AT[ID]
    Alen = ID_Alen[ID]
    X.append(AT)
    Y.append(Alen)
    Z.append(score2)

fig = plt.figure()
plt.scatter(X, Y, c=Z)
plt.show()
plt.close()
"""

"""
# use LinModel module
seq_list, score_list = [], []
count_list = []

for ID in ID_seq:
    seq = ID_seq[ID].upper()
    if 'N' in seq:
        continue
    seq_list.append(seq)
    score = new_ID_score2[ID]
    #score = ID_score2[ID]
    score_list.append(score)
    count_list.append(0)

m = LinModel.SeqLinearModel(seq_list, score_list, count_list)
#group_freq = m.spectrum(MM_orders=[1], Kmer_k_b=[5,1], PolyA_b=1, GC_b=1, Harmonic=False, gnum=10, sym=False, norm=False)
group_freq = m.spectrum(MM_orders=False, Kmer_k_b=[2,1], PolyA_b=False, GC_b=False, Harmonic=False, gnum=1000, sym=False, norm=False)


fig = plt.figure()
for i in range(len(group_freq)):
    freq = group_freq[i]
    din_value = freq['Kmer0']
    if i == 0:
        keys = din_value.keys()
    temp = [din_value[key] for key in keys]
    alpha = float(i)/len(group_freq)
    plt.plot(temp, 'kx-', label="group" + str(i), alpha=alpha)
    plt.xticks(range(len(keys)), keys)

#plt.legend()
plt.show()
plt.close()

for key in keys:
    Y = []
    for i in range(len(group_freq)):
        freq = group_freq[i]
        din_value = freq['Kmer0']
        value = din_value[key]
        Y.append(value)
    fig = plt.figure()
    plt.plot(range(len(Y)), Y, 'x-')
    plt.title(key)
    plt.show()
    plt.close()
"""
