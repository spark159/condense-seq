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
import statis

def rescale (value_list, old_st, old_ed, new_st, new_ed):
    output = []
    for value in value_list:
        assert value >= old_st and value <= old_ed
        new_value = new_st + (new_ed - new_st)*float(value-old_st)/(old_ed-old_st)
        output.append(new_value)
    return output

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

# parameters
clip = 0
NCPlen = 101

# group by AT content
AT_IDs = {}
for ID in ID_seq:
    seq = ID_seq[ID]
    middle = len(seq)/2
    seq = seq[middle-NCPlen/2:middle+NCPlen/2+1]
    AT = 1.0 - GC_content(seq)
    if AT not in AT_IDs:
        AT_IDs[AT] = []
    AT_IDs[AT].append(ID)

# check the correaltion between max poly (dA:dT) length and condensability
def mean_length (num_pos):
    total = 0.0
    mean = 0.0
    for num, pos in num_pos.items():
        mean += num*len(pos)
        total += len(pos)
    return mean/total

AT_pAcorr = {}
AT_pGcorr = {}
AT_weight = {}
for AT, IDs in AT_IDs.items():
    AT_weight[AT] = len(IDs)
    if len(IDs) < 5:
        AT_pAcorr[AT] = np.nan
        AT_pGcorr[AT] = np.nan
        continue
    #X1 = [poly_score(ID_seq[ID], nts='AT') for ID in IDs]
    #X2 = [poly_score(ID_seq[ID], nts='GC') for ID in IDs]
    X1 = [mean_length(poly_score(ID_seq[ID], nts='AT', pos=True)) for ID in IDs]
    X2 = [mean_length(poly_score(ID_seq[ID], nts='GC', pos=True)) for ID in IDs]

    Y = [ID_score1[ID] for ID in IDs]

    fig = plt.figure()
    plt.plot(X1, Y, '.')
    plt.show()
    plt.close()
    
    pAcorr = statis.get_spearman_corr(X1, Y)
    pGcorr = statis.get_spearman_corr(X2, Y)
    AT_pAcorr[AT] = pAcorr
    AT_pGcorr[AT] = pGcorr
    
# plot AT vs correlation
wmin, wmax = min(AT_weight.values()), max(AT_weight.values())
smin, smax = 10, 100
X1, Y1, S1 = [], [], []
X2, Y2, S2 = [], [], []
for AT in AT_weight:
    weight = AT_weight[AT]
    size = smin + (smax-smin)*float(weight-wmin)/(wmax-wmin)
    if not np.isnan(AT_pAcorr[AT]):
        X1.append(AT*100)
        Y1.append(AT_pAcorr[AT])
        S1.append(size)
    if not np.isnan(AT_pGcorr[AT]):
        X2.append(AT*100)
        Y2.append(AT_pGcorr[AT])
        S2.append(size)

fig = plt.figure(figsize=(7, 3))
plt.scatter(X1, Y1, s=S1, color='tab:red', edgecolor='k', label='poly(dA:dT)')
plt.scatter(X2, Y2, s=S2, color='tab:blue', edgecolor='k', label='poly(dG:dC)')
plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
plt.xlabel("AT content (%)", fontsize=8)
plt.ylabel("Spearman correlation", fontsize=8)
plt.title("Homopolymer length vs Condensability", fontsize=10)
plt.gca().tick_params('x', labelsize=6)
plt.gca().tick_params('y', labelsize=8)
plt.legend(fontsize=8)
plt.savefig('stratification.svg', format='svg', bbox_inches='tight')
plt.close()
