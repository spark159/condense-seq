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
from scipy import stats

def GC_content(seq):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    return (num/float(len(seq)))*100


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

def get_kmer_count(seq, knum=None, kmer=None):
    seq = seq.upper()
    if kmer:
        kmer, knum = kmer.upper(), len(knum)
        count = 0
        for i in range(len(seq)-knum+1):
            if seq[i:i+knum] == kmer:
                count +=1
        return count
    kmer_count = {}
    for i in range(len(seq)-knum+1):
        kmer = seq[i:i+knum]
        if "N" in kmer:
            continue
        if kmer not in kmer_count:
            kmer_count[kmer] = 0
        kmer_count[kmer] += 1
    return kmer_count
        
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


path = './data/'
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+"hg19_chr1_NCP_ics_anot.cn")
ID_score1 = name_ID_value['data/sp_spd_tests_detail/sp7']
ID_seq = name_ID_value['Sequence']
ID_AT = name_ID_value['ATcontent']
ID_CpG = name_ID_value['CpGNumber']
ID_me = name_ID_value['meGCNumber']

knum = 2
all_kmers = all_path(knum, "ATCG")

#ID_kmer_count = {}
kmer_counts, kmer_scores = {}, {}
for ID in ID_seq:
    seq = ID_seq[ID]
    kmer_count = get_kmer_count(seq, knum=knum)
    #ID_kmer_count[ID] = kmer_count
    score = ID_score1[ID]
    for kmer, count in kmer_count.items():
        if kmer not in kmer_counts:
            kmer_counts[kmer] = []
        if kmer not in kmer_scores:
            kmer_scores[kmer] = []
        kmer_counts[kmer].append(count)
        kmer_scores[kmer].append(score)

kmer_corr = {}
kmer_pvalue = {}
for kmer in all_kmers:
    X, Y = kmer_counts[kmer], kmer_scores[kmer]
    corr, pvalue = stats.spearmanr(X, Y)
    #corr = statis.get_corr(X, Y)
    kmer_corr[kmer] = corr
    kmer_pvalue[kmer] = pvalue

ATkmer = []
X, Y = [], []
for kmer, corr in kmer_corr.items():
    ATcontent = 100 - GC_content(kmer) 
    #X.append(ATcontent)
    #Y.append(corr)
    X.append(corr)
    pvalue = kmer_pvalue[kmer]
    Y.append(pvalue)
    #Y.append(-np.log(pvalue))
    ATkmer.append([ATcontent, kmer])

fig = plt.figure()
plt.scatter(X, Y, s=5)
plt.show()
plt.close()

ATkmer = sorted(ATkmer)
data = []
xlabels = []
for AT, kmer in ATkmer:
    data.append(kmer_corr[kmer])
    xlabels.append(kmer)

fig = plt.figure()
plt.bar(range(len(data)), data, width=0.5, color='g')
plt.xticks(range(len(data)), xlabels, rotation=0)
plt.ylabel("Pearson correlation")
plt.axhline(y=0, color='k', linestyle='--')
#plt.title("Conditional correlation with Condensability")
plt.tight_layout()
#plt.savefig("bar_cdcorr.png",bbox_inches='tight')
plt.show()
plt.close()

sys.exit(1)

fig = plt.figure()
plt.plot(X, Y, 'b.', alpha=0.2)
#plt.scatter(X, Y, c=C, cmap='jet', s=5, alpha=0.2)
plt.xlabel("AT content (%)")
plt.ylabel("Correlation with condensability")
plt.title("5-mers VS condensability")
plt.show()
plt.close()
    

"""
ID_din_count = {}
for ID in ID_seq:
    seq = ID_seq[ID]
    din_count = get_dincount(seq)
    ID_din_count[ID] = din_count

    
# dinucleotide counting
din_count_scores = {}
for ID in ID_din_count:
    score = ID_score1[ID]
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
    #plt.show()
    plt.close()

  
# Partition by condensability
group_IDs = statis.quantile (ID_score1, 7)

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

# use LinModel module
seq_list, score_list = [], []
count_list = []

for ID in ID_seq:
    seq = ID_seq[ID].upper()
    dyad = len(seq)/2
    if 'N' in seq:
        continue
    score = ID_score1[ID]
    for i in range(-1,2):
        seq_list.append(seq[dyad+i-147/2:dyad+i+147/2+1])
        score_list.append(score)
        count_list.append(0)

m = LinModel.SeqLinearModel(seq_list, score_list, count_list)
group_freq = m.spectrum(MM_orders=[1], Kmer_k_b=False, PolyA_b=False, GC_b=False, Harmonic=False, gnum=10, sym=False, norm=False)
#group_freq = m.spectrum(MM_orders=False, Kmer_k_b=[2,1], PolyA_b=False, GC_b=False, Harmonic=False, gnum=1000, sym=False, norm=False)

def get_ATGC_sig (freq):
    freq_MM1 = freq['MM1']
    length = len(freq_MM1)
    nts = all_path(2, 'ATCG')
    AT_sig, GC_sig = np.zeros(length), np.zeros(length)
    for nt in nts:
        row = [ freq_MM1[i][nt] for i in range(length)]
        if nt in ['AA', 'AT', 'TA', 'TT']:
            AT_sig += np.asarray(row)
        if nt in ['GG', 'GC', 'CG', 'CC']:
            GC_sig += np.asarray(row)
    for i in range(length):
        AT_sig[i] = float(AT_sig[i]) / sum(freq_MM1[i].values())
        GC_sig[i] = float(GC_sig[i]) / sum(freq_MM1[i].values())
    return AT_sig, GC_sig


fig = plt.figure()
for i in range(len(group_freq)):
    freq = group_freq[i]
    AT_sig, GC_sig = get_ATGC_sig(freq)
    #fig = plt.figure()
    plt.plot(AT_sig, label='AT-rich ' + str(i))
    #plt.plot(AT_sig, 'hotpink', label='AT-rich ' + str(i))
    #plt.plot(GC_sig, 'lightskyblue', label='GC-rich')
plt.xticks([147/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
plt.xlabel("Super Helical Location")
plt.ylabel("Relative frequency")
plt.legend()
#plt.ylim([0.22, 0.28])
#plt.savefig('ATGCperiod_' + "slide" + '.png')
plt.show()
plt.close()



"""

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
