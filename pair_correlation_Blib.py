#import load_file
#import graphics
#import statis
import sys
import copy
#import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
#import Interval_dict
from scipy import stats
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.cm as cm
from scipy import fft
import matplotlib as mpl
import random
import pickle

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


# read loop-seq data file
def read_loop_seq (fname):
    field_ID_value = {}
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            fields = ['score', 'amplt', 'phase', 'seq']
            First=False
            continue
        ID, seq = int(cols[0]), cols[1]
        values = [float(value) for value in cols[-3:]] + [seq]
        for field, value in zip(fields, values):
            if field not in field_ID_value:
                field_ID_value[field] = {}
            field_ID_value[field][ID] = value
    return field_ID_value

# get GC content of sequence
def GC_content(seq, percent=False):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    if percent:
        return (num/float(len(seq)))*100
    return (num/float(len(seq)))

# get reverse complementary of sequence
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

# get all possible N-mers
def all_path (N, states='ATCG'):
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output

# count kmer-pair occurrence of sequence
def get_pair_dist_count (seq, klen, max_dist=sys.maxint):
    pair_dist_count = {}
    dist_count = {k:0 for k in range(klen, min(len(seq)-klen, max_dist)+1)}
    seq = seq.upper()
    for i in range(len(seq)-klen):
        for j in range(i+1, min(len(seq)-klen+1, i+max_dist+1)):
            dist = j-i
            if dist < klen:
                continue
            pair_list = []

            kmer1 = seq[i:i+klen]
            kmer2 = seq[j:j+klen]
            pair = tuple(sorted([kmer1, kmer2]))
            pair_list.append(pair)

            #if rev:
            #    kmer1 = rev_comp(kmer1)
            #    kmer2 = rev_comp(kmer2)
            #    pair = tuple(sorted([kmer1, kmer2]))
            #    pair_list.append(pair)

            for pair in pair_list:
                if pair not in pair_dist_count:
                    pair_dist_count[pair] = copy.deepcopy(dist_count)
                pair_dist_count[pair][dist] +=1

    ## make probability mass function
    #if norm:
    #    for pair in pair_dist_count:
    #        total = float(sum(pair_dist_count[pair].values()))
    #        for dist in pair_dist_count[pair]:
    #            pair_dist_count[pair][dist] /= total

    return pair_dist_count

# get weighted averaged kmer-pair occurence of sequence pool
def get_mean_pair_dist_count (ID_seq, ID_weight, klen, max_dist=sys.maxint, rev=True, norm=False):
    total_weight = sum(ID_weight.values())
    mean_pair_dist_count = {}
    for ID in ID_seq:
        seq = ID_seq[ID]
        prob = float(ID_weight[ID]) / total_weight
        pair_dist_count = get_pair_dist_count(seq, klen, max_dist=max_dist, rev=rev, norm=norm)
        for pair in pair_dist_count:
            for dist, count in pair_dist_count[pair].items():
                if pair not in mean_pair_dist_count:
                    mean_pair_dist_count[pair] = {}
                if dist not in mean_pair_dist_count[pair]:
                    mean_pair_dist_count[pair][dist] = 0.0
                mean_pair_dist_count[pair][dist] += prob*count
    return mean_pair_dist_count

# combine all kmer-pair occurent of sequence pool
def combine_pair_dist_count (seq_list, klen, max_dist=sys.maxint, rev=True, norm=False):
    total_pair_dist_count = {}

    if rev:
        total_seq_list = []
        for seq in seq_list:
            total_seq_list.append(seq)
            total_seq_list.append(rev_comp(seq))
    else:
        total_seq_list = seq_list
        
    for seq in total_seq_list:
        pair_dist_count = get_pair_dist_count(seq, klen, max_dist=max_dist)
        for pair in pair_dist_count:
            for dist, count in pair_dist_count[pair].items():
                if pair not in total_pair_dist_count:
                    total_pair_dist_count[pair] = {}
                if dist not in total_pair_dist_count[pair]:
                    total_pair_dist_count[pair][dist] = 0
                total_pair_dist_count[pair][dist] += count

    # make probability mass function
    if norm:
        for pair in total_pair_dist_count:
            total = float(sum(total_pair_dist_count[pair].values()))
            for dist in total_pair_dist_count[pair]:
                total_pair_dist_count[pair][dist] /= total

    return total_pair_dist_count
        


# load loop-seq data
#path = ''
#fname = 'loop_seq_random_library.txt'
#field_ID_value = read_loop_seq(path+fname)
#ID_score = field_ID_value['score']
#ID_seq = field_ID_value['seq']
#NCPlen = 100
#clip = 25 # clip off the both ends of sequence 
#max_dist = 48
#note = '_loop'
#vmin, vmax = -0.15, 0.15


# load condense-seq data 
#path = ''
#fname = "H1_NCP_sp_chr1_anot.cn"
#ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+fname)
#ID_score = name_ID_value["work/2021_06_07_H1_sp_detail/H1-NCP-sp-8"]
#ID_seq = name_ID_value['Sequence']
#NCPlen = 147
#clip = 10 # clip off the both ends of sequence 
#max_dist = 48
#note = '_condense'
#vmin, vmax = -0.05, 0.05

# load YW library condense-seq data
tname_ID_score = read_score("YWlib_sp_score.txt")
ID_score = tname_ID_score[3]
ID_seq = read_ref('Blib.ref')

NCPlen = 101
clip = 0 # clip off the both ends of sequence 
klen = 1
max_dist = 10
note = '_YW'

p_num = 5


# slice out sequence
for ID in ID_seq:
    seq = ID_seq[ID]
    middle = len(seq)/2
    seq = seq[middle-NCPlen/2:middle+NCPlen/2+1]
    ID_seq[ID] = seq[clip:len(seq)-clip]


# get mean kmer-pair occurrence with uniform sampling
repeat = 1
sample_size = len(ID_seq)
control_pair_dist_count = {}
for i in range(repeat):
    seq_list = random.sample(ID_seq.values(), sample_size) 
    pair_dist_count = combine_pair_dist_count (seq_list, klen, max_dist=max_dist, norm=True)

    for pair in pair_dist_count:
        for dist, count in pair_dist_count[pair].items():
            if pair not in control_pair_dist_count:
                control_pair_dist_count[pair] = {}
            if dist not in control_pair_dist_count[pair]:
                control_pair_dist_count[pair][dist] = 0.0
            control_pair_dist_count[pair][dist] += count / float(repeat)

print "Done for control"

# quantile the group
p_IDs = quantile(ID_score, p_num)

pair_p_corrf = {}
for i in range(p_num):
    IDs = p_IDs[i]
    seq_list = [ID_seq[ID] for ID in IDs]
    
    # get mean kmer-pair occurrence with weighted sampling
    pair_dist_count = combine_pair_dist_count (seq_list, klen, max_dist=max_dist, norm=True)
    print "Done for data"

    # compute pair-correlation function
    pair_corrf = {}
    for pair in pair_dist_count:
        for dist in pair_dist_count[pair]:
            count = pair_dist_count[pair][dist]
            norm = control_pair_dist_count[pair][dist]
            if pair not in pair_corrf:
                pair_corrf[pair] = {}
            pair_corrf[pair][dist] = float(count)/norm
    print "pair corrf is computed"

    for pair in pair_corrf:
        corrf = pair_corrf[pair]
        if pair not in pair_p_corrf:
            pair_p_corrf[pair] = []
        pair_p_corrf[pair].append(corrf)
    
    # plot pair-correlation function
    fig = plt.figure()
    for pair in pair_corrf:
        X, Y = [], []
        for dist, corr in sorted(pair_corrf[pair].items()):
            X.append(dist)
            Y.append(corr)
        plt.plot(X, Y, '.-', label=pair)
    plt.legend()
    #plt.show()
    plt.close()

color_list = np.linspace(0.3, 1, num=p_num)
for pair in pair_p_corrf:
    pair_seq = ''.join(list(pair))
    rev_pair_seq = rev_comp(pair_seq)

    if pair_seq == rev_pair_seq:
        palindrome = True
    else:
        palindrome = False

    AT = 1.0 - GC_content(pair_seq)
    if AT > 0.5:
        if not palindrome:
            cmapname = 'Reds'
        else:
            cmapname = 'Purples'
    elif AT < 0.5:
        if not palindrome:
            cmapname = 'Blues'
        else:
            cmapname = 'Greens'
    else:
        cmapname = 'Greys'
    cmap = cm.get_cmap(cmapname)

    if palindrome:
        title = '%s-%s' % (pair_seq[0], pair_seq[1])
    else:
        title = '%s-%s / %s-%s' % (pair_seq[0], pair_seq[1], rev_pair_seq[0], rev_pair_seq[1])

    fig = plt.figure(figsize=(1.4, 1.2))
    for i in range(p_num):
        corrf = pair_p_corrf[pair][i]
        X, Y = [], []
        for dist, corr in sorted(corrf.items()):
            X.append(dist)
            Y.append(corr)
        plt.plot(X, Y, '.-', color=cmap(color_list[i]), lw=1.5)
    #plt.axhline(y=1, color='k', linestyle='--', alpha=0.5)    
    plt.title(title, fontsize=10)
    plt.ylim([1-0.08, 1+0.08])
    plt.xticks(X, [str(x) for x in X])
    plt.gca().tick_params('x', labelsize=6)
    plt.gca().tick_params('y', labelsize=6)
    #plt.xlabel("Separation (bp)", fontsize=8)
    #plt.ylabel("Pair corrf", fontsize=8)
    plt.savefig('paircorr_' + pair[0] + pair[1] + '.svg', format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()
    

        


    
"""
# get mean kmer-pair occurrence with uniform sampling
ID_weight = {ID:1.0 for ID in ID_score}
#control_pair_dist_count = get_mean_pair_dist_count(ID_seq, ID_weight, klen=klen, max_dist=max_dist, norm=True, rev=False)
print "Done for control"

# get mean kmer-pair occurrence with weighted sampling
ID_weight = {ID:1.0 for ID in ID_score}
#ID_weight = {ID:np.exp(ID_score[ID]) for ID in ID_score}
pair_dist_count = get_mean_pair_dist_count(ID_seq, ID_weight, klen=klen, max_dist=max_dist, norm=True, rev=False)
print "Done for data"

# compute pair-correlation function
pair_corrf = {}
for pair in pair_dist_count:
    for dist in pair_dist_count[pair]:
        count = pair_dist_count[pair][dist]
        norm = control_pair_dist_count[pair][dist]
        if pair not in pair_corrf:
            pair_corrf[pair] = {}
        pair_corrf[pair][dist] = float(count)/norm
print "pair corrf is computed"

# plot pair-correlation function
fig = plt.figure()
for pair in pair_corrf:
    X, Y = [], []
    for dist, corr in sorted(pair_corrf[pair].items()):
        X.append(dist)
        Y.append(corr)
    plt.plot(X, Y, '.-', label=pair)
plt.legend()
plt.show()
plt.close()



# group by AT content
AT_IDs = {}
for ID in ID_seq:
    seq = ID_seq[ID]
    AT = 1.0 - GC_content(seq)
    if AT not in AT_IDs:
        AT_IDs[AT] = []
    AT_IDs[AT].append(ID)

# compute correlation fucntion for each group
for AT, IDs in AT_IDs.items():
    if len(IDs) < 10:
        continue

    print 'AT content %f' % (AT)
    
    sID_seq = {ID:ID_seq[ID] for ID in IDs}

    # get mean kmer-pair occurrence with uniform sampling
    ID_weight = {ID:1.0 for ID in IDs}
    control_pair_dist_count = get_mean_pair_dist_count(sID_seq, ID_weight, klen=klen, max_dist=max_dist, norm=True)
    print "Done for control"

    # get mean kmer-pair occurrence with weighted sampling
    ID_weight = {ID:np.exp(ID_score[ID]) for ID in IDs}
    pair_dist_count = get_mean_pair_dist_count(sID_seq, ID_weight, klen=klen, max_dist=max_dist, norm=True)
    print "Done for data"

    # compute pair-correlation function
    pair_corrf = {}
    for pair in pair_dist_count:
        for dist in pair_dist_count[pair]:
            count = pair_dist_count[pair][dist]
            norm = control_pair_dist_count[pair][dist]
            if pair not in pair_corrf:
                pair_corrf[pair] = {}
            pair_corrf[pair][dist] = float(count)/norm
    print "pair corrf is computed"

    # plot pair-correlation function
    fig = plt.figure()
    for pair in pair_corrf:
        X, Y = [], []
        for dist, corr in sorted(pair_corrf[pair].items()):
            X.append(dist)
            Y.append(corr)
        plt.plot(X, Y, '.-', label=pair)
    plt.legend()
    plt.show()
    plt.close()


path = ""
fname = "random_pair_dist_count" + note
try:
    random_pair_dist_count = pickle.load(open(fname + ".pickle", "rb"))
except:
    random.seed(0)
    random_num = 10000
    random_seq_list = []
    for i in range(random_num):
        seq = ""
        for k in range(NCPlen-2*clip):
            seq += random.choice('ATCG')
        random_seq_list.append(seq)
    random_pair_dist_count = get_pair_dist_count(random_seq_list, max_dist=max_dist, norm=True)
    pickle.dump(random_pair_dist_count, open(fname + ".pickle", "wb"))
    del random_seq_list
print "random sequence counting is done"


# get kmer-pair correlation function and helicity score for library
# save the correlation between hscore and loop/consensability score
path = ""
fname = "pair_corrf" + note
try:
    pair_corrf = pickle.load(open(fname + ".pickle", "rb"))
except:
    print "start reading sequence data"
    print "total %d sequences" % (len(ID_seq))
    count = 0 
    log_counter = -1
    pair_corrf = {}
    #pair_ID_hscore = {}
    for ID, seq in ID_seq.items():
        # adjust sequence
        middle = len(seq)/2
        seq = seq[middle-NCPlen/2:middle+NCPlen/2+1]
        assert len(seq) == NCPlen
        seq = seq[clip:len(seq)-clip].upper()
        if set(seq) > set('ATCG'):
            continue
        # counting the sequences
        count +=1
        if int(np.log10(count)) > log_counter:
            print count
            log_counter +=1
    
        pair_dist_count = get_pair_dist_count([seq], max_dist=max_dist)

        # compute pair-correaltion function
        score = ID_score[ID]
        prob = 1.0 - 2**(-score)
        
        for pair in pair_dist_count:
            for dist in pair_dist_count[pair]:
                #norm = random_pair_dist_count[pair][dist]
                count = pair_dist_count[pair][dist]
                if pair not in pair_corrf:
                    pair_corrf[pair] = {}
                if dist not in pair_corrf[pair]:
                    pair_corrf[pair][dist] = [0.0 for k in range(klen, min(len(seq)-klen, max_dist)+1)]
                pair_corrf[pair][dist] += prob*count/norm

    pickle.dump(pair_corrf, open(fname + ".pickle", "wb"))
print "get pair correlation function is done"        


# plot kmer-pair correlation matrix
GCcontent_din = sorted([(GC_content(din), din) for din in all_path(2, 'ATCG')])
all_din = [din for GCcontent, din in GCcontent_din]

corr_img = np.zeros((len(all_din), len(all_din)))
for i in range(len(all_din)):
    for j in range(i, len(all_din)):
        din1, din2 = all_din[i], all_din[j]
        pair = tuple(sorted([din1, din2]))
        corr_img[i][j] = pair_scorr[pair]
        corr_img[j][i] = pair_scorr[pair]

fig = plt.figure(figsize=(3,3))
plt.imshow(corr_img, cmap='Spectral_r', vmin=vmin, vmax=vmax)
plt.xticks(range(len(all_din)), all_din, fontsize=8, rotation=48)
plt.yticks(range(len(all_din)), all_din, fontsize=8)
if note.endswith('loop'):
    title = "Dinucleotide Helicity VS %s" % ('Cyclizability')
elif note.endswith('condense'):
    title = "Dinucleotide Helicity VS %s" % ('Condensability')
plt.title(title, fontsize=8)
cbar = plt.colorbar(pad=0.05, ticks=[vmin, vmax], fraction=0.04)
cbar.ax.set_yticklabels([str(vmin), str(vmax)], fontsize=5)
cbar.ax.set_ylabel('Spearman correlation', rotation=-90, va="bottom", labelpad=-15, fontsize=8)
plt.savefig(title + '.svg', format='svg', bbox_inches='tight')
#plt.show()
plt.close()
"""
