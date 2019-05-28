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

#ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("data/hg19_chr1_171_everything_anot.cn")
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("data/hg19_chr1_171_icdseq_anot.cn")
ID_score1 = name_ID_value['work/condense_seq/sp9_hg19_chr1']
ID_score2 = name_ID_value['work/condense_seq/sp10_hg19_chr1']
ID_AT = name_ID_value['ATcontent']
ID_seq = name_ID_value['Sequence']

# define sticky and fluffy nucleosomes
new_ID_score2 = statis.neutralize_score_by_target(ID_score2, ID_AT)
fluffy_IDs, middle_IDs, sticky_IDs = statis.quantile (new_ID_score2, 3, [5, 90, 5])
total_IDs = ID_pos.keys()

for ID in ID_AT:
    ID_AT[ID] = 100*ID_AT[ID]

graphics.Scatter_plot (ID_AT, ID_score2, ylim=[-3, 3.5], note='raw')
graphics.Scatter_plot (ID_AT, new_ID_score2, ylim=[-5, 5.5], ylabel='Condensability - AT dependence (A.U.)', note='neutral')

X1, Y1 = [], []
for ID in sticky_IDs:
    X1.append(ID_AT[ID])
    Y1.append(new_ID_score2[ID])
X2, Y2 = [], []
for ID in middle_IDs:
    X2.append(ID_AT[ID])
    Y2.append(new_ID_score2[ID])
X3, Y3 = [], []
for ID in fluffy_IDs:
    X3.append(ID_AT[ID])
    Y3.append(new_ID_score2[ID])

fig = plt.figure()
plt.plot(X1, Y1, '.', label='Sticky nucleosomes', alpha=0.01)
plt.plot(X2, Y2, 'k.', alpha=0.01)
plt.plot(X3, Y3, '.', label='Fluffy nucleosomes', alpha=0.01)
leg = plt.legend(loc='best', numpoints=1, prop={'size': 10})
for lh in leg.legendHandles:
    lh._legmarker.set_markersize(10)
    lh._legmarker.set_alpha(1)
plt.xlabel('AT content (%)')
plt.ylabel("Condensability - AT dependence (A.U.)")
#plt.title("Chromosome 1")
plt.xlim([0, 100])
plt.ylim([-5, 5.5])
plt.savefig("StickyFluffy" + "_scatter.png")
#plt.show()
plt.close()

# compare sequence feature
ID_din_count = {}
for ID in ID_seq:
    seq = ID_seq[ID]
    din_count = get_dincount(seq)
    ID_din_count[ID] = din_count

ID_TA = {}
for ID in ID_seq:
    seq = ID_seq[ID]
    count = get_dincount(seq, din="TA")
    ID_TA[ID] = count

ID_polyAT, ID_polyGC = {}, {}
for ID in ID_seq:
    seq = ID_seq[ID]
    num_pos = poly_score(seq, nts='AT', pos=True)
    score = 0.0
    for num, pos in num_pos.items():
        if num < 3:
            continue
        score += len(pos)*(num**2)
    ID_polyAT[ID] = score
    num_pos = poly_score(seq, nts='GC', pos=True)
    score = 0.0
    for num, pos in num_pos.items():
        if num < 3:
            continue
        score += len(pos)*(num**2)
    ID_polyGC[ID] = score


count1 = [ID_polyAT[ID] for ID in sticky_IDs]
count2 = [ID_polyAT[ID] for ID in total_IDs]
count3 = [ID_polyAT[ID] for ID in fluffy_IDs]

fig = plt.figure()
data = [count1, count2, count3]
plt.boxplot(data, 0, "")
plt.xticks(range(1, len(data)+1), ['Sticky', 'Total', 'Fluffy'])
plt.title("Poly-A score")
plt.savefig("StickyFluffy_box_" + "polyA" + ".png")
#plt.show()
plt.close()


all_din = all_path(2)
for din in all_din:
    count1, count2, count3 = [], [], []

    for ID in sticky_IDs:
        try:
            count = ID_din_count[ID][din]
        except:
            count = 0
        count1.append(count)

    for ID in total_IDs:
        try:
            count = ID_din_count[ID][din]
        except:
            count = 0
        count2.append(count)

    for ID in fluffy_IDs:
        try:
            count = ID_din_count[ID][din]
        except:
            count = 0
        count3.append(count)

    fig = plt.figure()
    data = [count1, count2, count3]
    plt.boxplot(data, 0, "")
    plt.xticks(range(1, len(data)+1), ['Sticky', 'Total', 'Fluffy'])
    plt.title(din + " Dinucleotide count")
    plt.savefig("StickyFluffy_box_" + din + ".png")
    #plt.show()
    plt.close()


# compare all epigenetic marks
for name in name_ID_value:
    if name == "Sequence":
        continue
    ID_value = name_ID_value[name]
    count1, count2, count3 = [], [], []

    for ID in sticky_IDs:
        count = ID_value[ID]
        count1.append(count)

    for ID in total_IDs:
        count = ID_value[ID]
        count2.append(count)

    for ID in fluffy_IDs:
        count = ID_value[ID]
        count3.append(count)

    fig = plt.figure()
    data = [count1, count2, count3]
    plt.boxplot(data, 0, "")
    plt.xticks(range(1, len(data)+1), ['Sticky', 'Total', 'Fluffy'])
    plt.title(name)
    plt.savefig("StickyFluffy_box_" + name.split('/')[-1] + ".png")
    #plt.show()
    plt.close()


# gene domain fraction difference
sticky_pos = []
for ID in sticky_IDs:
    sticky_pos.append(ID_pos[ID])
total_pos = []
for ID in total_IDs:
    total_pos.append(ID_pos[ID])
fluffy_pos = []
for ID in fluffy_IDs:
    fluffy_pos.append(ID_pos[ID])

genome_size = load_file.read_genome_size("data/hg19.fa")
ID_field_values, field_ID_values = load_file.read_GTF ("data/Homo_sapiens.GRCh37.87.gtf", "chr1", mode="both")

TSS_range = 500
TTS_range = 500
Prom_range = 5000
dID_interval = {}
for ID in ID_field_values:
    field_values = ID_field_values[ID]
    if field_values['geneType'] != "protein_coding":
        continue
    TSS = field_values['TSS']
    TTS = field_values['TTS']
    TSS_st = max(TSS - TSS_range, 0)
    TSS_ed = min(TSS + TSS_range, genome_size['chr1'])
    TTS_st = max(TTS - TTS_range, 0) 
    TTS_ed = min(TTS + TTS_range, genome_size['chr1'])
    strand = field_values['strand']
    if strand == '+':
       Body_st = TSS_ed
       Body_ed = TTS_st
       Prom_st = max(TSS - Prom_range, 0)
       Prom_ed = TSS_st
    if strand == '-':
       Body_st = TTS_ed
       Body_ed = TSS_st
       Prom_st = TSS_ed
       Prom_ed = min(TSS + Prom_range, genome_size['chr1'])
    if TSS_ed - TSS_st > 0:
        dID_interval["TSS:" + ID] = [TSS_st, TSS_ed]
    if TTS_ed - TTS_st > 0:
        dID_interval["TTS:" + ID] = [TTS_st, TTS_ed]
    if Body_ed - Body_st > 0:
        dID_interval["Body:" + ID] = [Body_st, Body_ed]
    if Prom_ed - Prom_st > 0:
        dID_interval["Prom:" + ID] = [Prom_st, Prom_ed]
dinterval_dict = Interval_dict.double_hash(dID_interval, 100000, genome_size['chr1'])

ID_dtype = {}

for ID in ID_pos:
    pos = ID_pos[ID]
    dIDs = dinterval_dict.find(pos)
    assert ID not in ID_dtype
    ID_dtype[ID] = []
    if len(dIDs) <= 0:
        ID_dtype[ID].append("Intergenic")
        continue
    for dID in dIDs:
        domain = dID.split(':')[0]
        if domain not in ID_dtype[ID]:
            ID_dtype[ID].append(domain)
                
domain_count1 = {}
domain_count2 = {}
domain_count3 = {}

for ID in sticky_IDs:
    domains = ID_dtype[ID]
    for domain in domains:
        if domain not in domain_count1:
            domain_count1[domain] = 0
        domain_count1[domain] += 1

for ID in total_IDs:
    domains = ID_dtype[ID]
    for domain in domains:
        if domain not in domain_count2:
            domain_count2[domain] = 0
        domain_count2[domain] += 1

for ID in fluffy_IDs:
    domains = ID_dtype[ID]
    for domain in domains:
        if domain not in domain_count3:
            domain_count3[domain] = 0
        domain_count3[domain] += 1

X = ['Prom', 'TSS', 'Body', 'TTS', 'Intergenic']

Y1, Y2, Y3 = [], [], []
for name in X:
    try:
        Y1.append(domain_count1[name])
    except:
        Y1.append(0)
    try:
        Y2.append(domain_count2[name])
    except:
        Y2.append(0)
    try:
        Y3.append(domain_count3[name])
    except:
        Y3.append(0)


Y1 = np.asarray([ float(value)/sum(Y1) for value in Y1 ])
Y2 = np.asarray([ float(value)/sum(Y2) for value in Y2 ])
Y3 = np.asarray([ float(value)/sum(Y3) for value in Y3 ])

diff1 = (Y1-Y2)/Y2
diff2 = (Y3-Y2)/Y2

fig = plt.figure()
plt.bar([i-0.1 for i in range(len(X))], diff1, width=0.2, label='Sticky nucleosomes')
#plt.bar([i for i in range(len(X))], Y2-Y2, width=0.2, label='total')
plt.bar([i+0.1 for i in range(len(X))], diff2, width=0.2, label='Fluffy nucleosomes')
plt.axhline(y=0, linestyle='--', color='k', linewidth=1)
plt.xticks(range(len(X)), X)
plt.ylabel("Fraction difference")
plt.title("Protein coding genes")
plt.legend()
plt.savefig("StickyFluffy_genecoding.png")
#plt.show()
plt.close()


# ChromHHM domain fraction difference
def read_chromHMM(fname, chr_target, change=False):
    state_intervals = {}
    for line in open(fname):
        cols = line.strip().split()
        chr, st, ed, state = cols
        if chr != chr_target:
            continue
        st, ed = int(st), int(ed)
        if change:
            state = change[state]
        if state not in state_intervals:
            state_intervals[state] = []
        state_intervals[state].append((st,ed))
    return state_intervals

#state_intervals = read_chromHMM("data/E098_15_coreMarks_mnemonics.bed", chr_target='chr1')
#name_dict = {'E1':'Polycomb', 'E2':'me2Het', 'E3':'Quies', 'E4':'me3Het', 'E5':'Tx', 'E6':'Enh', 'E7':'TssA'}
#state_intervals = read_chromHMM("data/38-Per_7_segments.bed", chr_target='chr1', change=name_dict)
name_dict = {'E1':'TssBiv', 'E2':'TssA', 'E3':'EnhA', 'E4':'TxWk', 'E5':'Tx', 'E6':'me3Het', 'E7':'Quies', 'E8':'me2Het', 'E9':'PcWk', 'E10':'Pc'}
state_intervals = read_chromHMM("data/38-Per_10_segments.bed", chr_target='chr1', change=name_dict)

dID_interval = {}
for state in state_intervals:
    intervals = state_intervals[state]
    for i in range(len(intervals)):
        dID = state + ':' + str(i)
        assert dID not in dID_interval
        dID_interval[dID] = intervals[i]

dinterval_dict = Interval_dict.double_hash(dID_interval, 10000, 250000000)

ID_dtype = {}

for ID in ID_pos:
    pos = ID_pos[ID]
    dIDs = dinterval_dict.find(pos)
    assert ID not in ID_dtype
    ID_dtype[ID] = []
    if len(dIDs) <= 0:
        ID_dtype[ID].append("None")
        continue
    for dID in dIDs:
        domain = dID.split(':')[0]
        if domain not in ID_dtype[ID]:
            ID_dtype[ID].append(domain)
                
domain_count1 = {}
domain_count2 = {}
domain_count3 = {}

for ID in sticky_IDs:
    domains = ID_dtype[ID]
    for domain in domains:
        if domain not in domain_count1:
            domain_count1[domain] = 0
        domain_count1[domain] += 1

for ID in total_IDs:
    domains = ID_dtype[ID]
    for domain in domains:
        if domain not in domain_count2:
            domain_count2[domain] = 0
        domain_count2[domain] += 1

for ID in fluffy_IDs:
    domains = ID_dtype[ID]
    for domain in domains:
        if domain not in domain_count3:
            domain_count3[domain] = 0
        domain_count3[domain] += 1

#X = ['1_TssA', '2_TssAFlnk', '3_TxFlnk', '4_Tx', '5_TxWk', '6_EnhG', '7_Enh', '8_ZNF/Rpts', '9_Het', '10_TssBiv', '11_BivFlnk', '12_EnhBiv', '13_ReprPC', '14_ReprPCWk', '15_Quies']
#X = ['TssA', 'Enh', 'Tx', 'me2Het', 'me3Het', 'Polycomb', 'Quies']
X = ['TssA', 'EnhA', 'Tx', 'TxWk', 'TssBiv', 'me2Het', 'me3Het', 'PcWk', 'Pc', 'Quies']


Y1, Y2, Y3 = [], [], []
for name in X:
    try:
        Y1.append(domain_count1[name])
    except:
        Y1.append(0)
    try:
        Y2.append(domain_count2[name])
    except:
        Y2.append(0)
    try:
        Y3.append(domain_count3[name])
    except:
        Y3.append(0)

Y1 = np.asarray([ 100*float(value)/sum(Y1) for value in Y1 ])
Y2 = np.asarray([ 100*float(value)/sum(Y2) for value in Y2 ])
Y3 = np.asarray([ 100*float(value)/sum(Y3) for value in Y3 ])

diff1 = (Y1-Y2)/Y2
diff2 = (Y3-Y2)/Y2
        
fig = plt.figure()
plt.bar([i-0.1 for i in range(len(X))], diff1, width=0.2, label='Sticky nucleosomes')
#plt.bar([i for i in range(len(X))], Y2-Y2, width=0.2, label='total')
plt.bar([i+0.1 for i in range(len(X))], diff2, width=0.2, label='Fluffy nucleosomes')
plt.axhline(y=0, linestyle='--', color='k', linewidth=1)
plt.xticks(range(len(X)), X, rotation=75)
plt.ylabel("Fraction difference")
plt.title("ChromHMM")
plt.legend()
plt.savefig("StickyFluffy" + "_chromHMM.png", bbox_inches='tight')
#plt.show()
plt.close()
