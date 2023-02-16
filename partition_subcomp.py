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

# read bin score file
def read_bin_score_new (fname, bin_size):
    name_chr_binID_score = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        if First:
            names = cols[4:]
            First = False
            continue
        _, chr, st, ed = cols[:4]
        st, ed = int(st), int(ed)
        scores = [float(score) for score in cols[4:]]
        binID = st / int(bin_size)
        for name, score in zip(names, scores):
            if name not in name_chr_binID_score:
                name_chr_binID_score[name] = {}
            if chr not in name_chr_binID_score[name]:
                name_chr_binID_score[name][chr] = {}
            name_chr_binID_score[name][chr][binID] = score
    return name_chr_binID_score

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

### parameters
path = "/home/spark159/../../storage/"

#exp_list = [('H1', 'NCP', 'sp', 8),
#            ('H1', 'NCP', 'spd', 6),
#            ('H1', 'NCP', 'CoH', 5),
#            ('H1', 'NCP', 'PEG', 6),
#            ('H1', 'NCP', 'Ca', 5),
#            ('H1', 'NCP', 'Mg', 5),
#            ('H1', 'NCP', 'HP1a', 3),
#            ('H1', 'NCP', 'HP1bSUV', 4),
#            ('H1', 'NCP', 'LKH', 3),
#            ('H1', 'NCP', 'Ki67', 4),
#            ('H1', 'NCP', 'FUS', 5)]

exp_list = [('H1', 'NCP', 'sp', 8),
            ('H1', 'NCP', 'spd', 6),
            ('H1', 'NCP', 'CoH', 5),
            ('H1', 'NCP', 'PEG', 6),
            ('H1', 'NCP', 'Ca', 5),
            ('H1', 'NCP', 'HP1a', 3),
            ('H1', 'NCP', 'LKH', 3),
            ('H1', 'NCP', 'Ki67', 4)]

exp_list = [('H1', 'NCP', 'sp', 8),
            ('H1', 'NCP', 'HP1a', 3),
            ('H1', 'NCP', 'LKH', 3),
            ('H1', 'NCP', 'Ki67', 4)]

exp_list = [('H1', 'NCP', 'sp', 8)]

# set domain names
domain_names = {'NSpeckle':['SON'],
                 'Trnx':['POLR2A', 'POLR2AphosphoS5', 'H3K9ac', 'H3K4me3', 'H3K27ac'],
                 'Polycomb':['CBX8', 'EZH2', 'RNF2', 'SUZ12'],
                 'Hetero':['H3K9me3', 'CBX5'],
                 'Nucleolus':['Nucleolar'],
                 'Lamin':['LaminB1'],
                 'Compartment':['eigen'],
                 'other':['ATcontent']}

domains = ['NSpeckle', 'Trnx', 'Polycomb', 'Hetero', 'Nucleolus', 'Lamin']

names = []
for domain in domains:
    names += domain_names[domain]

# binsize
bin_size = 10000

# other parameters
#dtype = 'zscore'
dtype = 'score'


# get enrichment zscore for each data
exp_name_p_mzscore = {}
for exp in exp_list:
    cell, sample, agent, tnum = exp
    print exp
    
    # load data
    fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', 'anot']) + '.txt'
    score_field = "%s-%s-%s-%d" % (cell, sample, agent, tnum)
    
    field_ID_value = load_file.read_tabular_file(fname, mode='col')
    ID_chr = field_ID_value['Chromosome']
    ID_score = field_ID_value[score_field]


    # Partition domain by value
    med = np.median(ID_score.values())
    std = np.std(ID_score.values())
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
    plt.hist(ID_score.values(), bins=1000)

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
        plt.text(x, 0 , num_rom[i+1], fontsize=12, va='center', ha='center')

    plt.xlim([med-4*std, med+4*std])
    #plt.title("Chromosome1", fontsize=8)
    #plt.xlabel("Condensability (A.U.)", fontsize=8)
    plt.ylabel("Bin Counts", fontsize=8)
    plt.gca().tick_params(axis='both', which='major', labelsize=5)
    plt.gca().tick_params(axis='both', which='minor', labelsize=5)
    plt.savefig('partition_%s.svg' % (agent), format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

    p_IDs = [[] for i in range(p_num)]
    for ID in ID_score:
        score = ID_score[ID]
        for i in range(p_num):
            st, ed = p_range[i]
            if score >= st and score < ed:
                break
        p_IDs[i].append(ID)


    # get mean zscore of each paritition for each data
    for name in names:
        p_values = [[] for i in range(p_num)]
        for i in range(p_num):
            for ID in p_IDs[i]:
                p_values[i].append(field_ID_value[name][ID])
                
        p_zscores = get_zscore(p_values)
        p_mzscore = [np.mean(zscores) for zscores in p_zscores]

        if exp not in exp_name_p_mzscore:
            exp_name_p_mzscore[exp] = {}

        exp_name_p_mzscore[exp][name] = p_mzscore

# plot mean zscore for each partitions
for exp in exp_list:
    cell, sample, agent, tnum = exp
    name_p_mzscore = exp_name_p_mzscore[exp]
    
    for domain in domains:
        img = []
        for name in domain_names[domain]:
            ylabels = []
            img.append(name_p_mzscore[name])
            ylabels.append(name)

        height = 0.1*len(domain_names[domain])
        width = 0.1*p_num

        fig = plt.figure()
        im = plt.imshow(img, cmap='bwr', vmin=-1.5, vmax=1.5)
        plt.yticks(range(len(img)), ylabels)
        plt.savefig("p_subcompt_%s_%s.svg" % (agent, name), format='svg', bbox_inches='tight')
        #plt.show()
        plt.close()


# plot colorbar only
fig = plt.figure(figsize=(1.2,1))
plt.subplot(1,2,1)
cbar = plt.colorbar(im, cax=plt.gca(), ticks=[-1.5, 1.5])
cbar.ax.set_yticklabels(['-1.5', '1.5'], fontsize=8)
cbar.ax.set_ylabel('Z-score', rotation=-90, va="bottom", fontsize=8)
plt.tight_layout()
plt.savefig('partition_cbar.svg', format='svg', bbox_inches='tight')
plt.close()
