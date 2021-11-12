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
plt.savefig("partition_hist.png")
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


# Compare sequence feature
for ID in ID_AT:
    ID_AT[ID] = ID_AT[ID]*100

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

p_ATs = [[] for i in range(p_num)]
din_p_sigs = {}
poly_p_sigs = {'polyA/T':[[] for i in range(p_num)], 'polyG/C':[[] for i in range(p_num)]}
for i in range(p_num):
    IDs = p_IDs[i]
    for ID in IDs:
        p_ATs[i].append(ID_AT[ID])
        poly_p_sigs['polyA/T'][i].append(ID_polyAT[ID])
        poly_p_sigs['polyG/C'][i].append(ID_polyGC[ID])
        seq = ID_seq[ID]
        din_count = get_dincount(seq, both=True) # survey both strands
        for din in din_count:
            if din not in din_p_sigs:
                din_p_sigs[din] = [[] for k in range(p_num)]
            din_p_sigs[din][i].append(din_count[din])

NCPlen = 147
MMfreq_list = []
for i in range(p_num):
    seq_list = []
    IDs = p_IDs[i]
    for ID in IDs:
        seq = ID_seq[ID].upper()
        if 'N' in seq:
            continue
        dyad = len(seq)/2
        for i in range(-1,2):
            seq_list.append(seq[dyad+i-NCPlen/2:dyad+i+NCPlen/2+1])
    freq, sample_num, mean, std, stdz_freq = stat_Markov(seq_list, NCPlen, 1)
    MMfreq_list.append(freq)
    

print "Sequence feature reading is done"


fig = plt.figure()
plt.xlabel('Partitions')
plt.ylabel('%')
plt.title('AT content')
plt.boxplot(p_ATs, 0, "")
plt.savefig("ATcontent"+"_pbox.png")
#plt.show()
plt.close()

for din in din_p_sigs:
    p_sigs = din_p_sigs[din]
    fig = plt.figure()
    plt.xlabel('Partitions')
    plt.ylabel('Counts')
    plt.title(din[0] + 'p' + din[1])
    plt.boxplot(p_sigs, 0, "")
    plt.savefig(din+"_pbox.png")
    #plt.show()
    plt.close()

for poly in poly_p_sigs:
    p_sigs = poly_p_sigs[poly]
    fig = plt.figure()
    plt.xlabel('Partitions')
    plt.ylabel('Score (A.U.)')
    plt.title(poly)
    plt.boxplot(p_sigs, 0, "")
    plt.savefig('_'.join(poly.split('/'))+"_pbox.png")
    #plt.show()
    plt.close()


def get_ATGC_sig (freq_MM1):
    #freq_MM1 = freq['MM1']
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

AT_sig_list, GC_sig_list = [], []
for i in range(p_num):
    freq = MMfreq_list[i]
    AT_sig, GC_sig = get_ATGC_sig(freq)
    AT_sig_list.append(AT_sig)
    GC_sig_list.append(GC_sig)

color_list = np.linspace(0.2, 1, num=p_num)
cmap1, cmap2 = cm.get_cmap("OrRd"), cm.get_cmap("GnBu")
cmap = cm.get_cmap("rainbow")
fig = plt.figure()
#plt.subplot(1,2,1)
ax1 = plt.gca()
ax2 = ax1.twinx()
for i in range(p_num):
    AT_sig, GC_sig = AT_sig_list[i], GC_sig_list[i]
    ax1.plot(AT_sig, color=cmap1(color_list[i]))
    ax2.plot(GC_sig, color=cmap2(color_list[i]))
ax1.set_ylabel('AA/AT/TA/TT freqeuncy', color='r')
ax2.set_ylabel('CC/CG/GC/GG freqeuncy', color='b')
ax1.set_xlabel("Super Helical Location")
ax1.tick_params('y', colors='r')
ax2.tick_params('y', colors='b')
line_list = [147/2 - i*20 for i in range(4)] + [147/2 + i*20 for i in range(1, 4)]
for line in line_list:
    plt.axvline(x=line, color='k', linestyle='--', alpha=0.25)
plt.xticks([147/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
#plt.ylabel("Relative frequency")
#plt.legend()
#plt.ylim([0.22, 0.28])
plt.savefig('ATGCperiod_' + "slide" + '.png')
plt.title("Dinucleotide periodicity")
#plt.show()
plt.close()

fig = plt.figure()
#plt.subplot(1,2,2)
img = [ [color, np.nan] for color in color_list ]
plt.imshow(img, cmap='OrRd', aspect=3)
img = [ [np.nan, color] for color in color_list ]
plt.imshow(img, cmap='GnBu', aspect=3)
plt.yticks(range(p_num), range(1, p_num+1))
ax = plt.gca()
ax.yaxis.tick_right()
#plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
plt.xticks([])
#plt.show()


    
#sys.exit(1)

# calculate z-scores
img1 = []
ylabels1 = []

p_zATs = get_zscore(p_ATs)
img1.append([np.mean(zATs) for zATs in p_zATs])
ylabels1.append('AT content')
    
din_p_zscores = {}
#din_list = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
din_list = ['AT', 'TA', 'CG', 'GC', 'AA', 'AC', 'GA', 'AG', 'CA', 'CC']
name_list = ['ApT', 'TpA', 'CpG', 'GpC', 'ApA or TpT', 'ApC or GpT', 'GpA or TpC', 'ApG or CpT', 'CpA or TpG', 'CpC or GpG']
for din, name in zip(din_list, name_list):
    din_p_zscores[din] = get_zscore(din_p_sigs[din])
    img1.append([np.mean(zscores) for zscores in din_p_zscores[din]])
    ylabels1.append(name)

poly_p_zscores = {}
for poly in sorted(poly_p_sigs):
    poly_p_zscores[poly] = get_zscore(poly_p_sigs[poly])
    img1.append([np.mean(zscores) for zscores in poly_p_zscores[poly]])
    ylabels1.append(poly)

# Compare epigenetic marks
ID_mefrac = {}
for ID in ID_CG:
    CG = ID_CG[ID]
    if CG <= 0:
        continue
    me = ID_me[ID]
    mefrac = float(me) / (CG)
    ID_mefrac[ID] = mefrac
    

me_p_sigs = {'CpG Number':[[] for k in range(p_num)], 'meGC Number':[[] for k in range(p_num)], 'meCpG Density':[[] for k in range(p_num)]}
chip_p_sigs = {}
for i in range(p_num):
    IDs = p_IDs[i]
    for ID in IDs:
        me_p_sigs['CpG Number'][i].append(ID_CG[ID])
        me_p_sigs['meGC Number'][i].append(ID_me[ID])
        try:
            me_p_sigs['meCpG Density'][i].append(ID_mefrac[ID])
        except:
            continue
        for name in name_ID_value:
            if not name.startswith('H'):
                continue
            if name not in chip_p_sigs:
                chip_p_sigs[name] = [[] for k in range(p_num)]
            ID_value = name_ID_value[name]
            chip_p_sigs[name][i].append(ID_value[ID])

print "Epigenetic marks reading is done"

for me in me_p_sigs:
    p_sigs = me_p_sigs[me]
    fig = plt.figure()
    plt.xlabel('Partitions')
    if me == 'meCpG Density':
        plt.ylabel('Fraction')
    else:
        plt.ylabel('Counts')
    plt.title(me)
    plt.boxplot(p_sigs, 0, "")
    plt.savefig(me+"_pbox.png")
    #plt.show()
    plt.close()

for chip in chip_p_sigs:
    p_sigs = chip_p_sigs[chip]
    fig = plt.figure()
    plt.xlabel('Partitions')
    plt.ylabel('Strength (A.U.)')
    plt.title(chip)
    plt.boxplot(p_sigs, 0, "")
    plt.savefig(chip+"_pbox.png")
    #plt.show()
    plt.close()

# calculate z-scores
img2 = []
ylabels2 = []

me_p_zscores = {}
me_names = ['CpG Number', 'meGC Number', 'meCpG Density']
for me in me_names:
    me_p_zscores[me] = get_zscore(me_p_sigs[me])
    img2.append([np.mean(zscores) for zscores in me_p_zscores[me]])
    ylabels2.append(me)

chip_p_zscores = {}
#chip_names = ['k27ac', 'k9ac', 'k4me3', 'k36me3_2', 'k9me2_2', 'k9me3_2', 'k27me3a_2']
#chip_names = ['k27ac', 'k9ac', 'k4me3', 'k36me3', 'k9me2', 'k9me3', 'k27me3a'
chip_names = ['H3k27ac', 'H3K9ac', 'H3K4me3', 'H2AFZ', 'H3K36me3', 'H3K9me3', 'H3K27me3']
for chip in chip_names:
    chip_p_zscores[chip] = get_zscore(chip_p_sigs[chip])
    img2.append([np.mean(zscores) for zscores in chip_p_zscores[chip]])
    ylabels2.append(chip)


# plot the z-scores of all features
img_list = [img1, img2]
ylabel_list = [ylabels1, ylabels2]
title_list = ['Sequence features', 'Epigenetic marks']

fig, (ax1, ax2, cax) = plt.subplots(ncols=3, gridspec_kw={"width_ratios":[1,1, 0.05]})
#fig, (ax1, ax2) = plt.subplots(ncols=2)
for ax, img, ylabels, title in zip([ax1, ax2], img_list, ylabel_list, title_list):
    im = ax.imshow(img, vmin=-1.5, vmax=1.5, cmap='coolwarm', aspect='auto')
    ax.set_xticks(range(p_num))
    ax.set_yticks(range(len(img)))
    ax.set_xticklabels(range(1, p_num+1))
    ax.set_yticklabels(ylabels)
    ax.set_title(title)
    plt.tight_layout()

cbar = fig.colorbar(im, cax=cax)
#cbar.ax.tick_params(labelsize=15) 
#cbar.ax.set_ylabel('Z-score', rotation=-90, va="bottom", fontsize=15)
plt.show()
plt.close()


fig = plt.figure()
plt.subplot(1,3,1)
plt.imshow(img1, vmin=-1.5, vmax=1.5, cmap='coolwarm', aspect='auto')
plt.xticks(range(p_num), range(1, p_num+1))
plt.yticks(range(len(img1)), ylabels1)
plt.title("Sequences")
plt.tight_layout() 
plt.subplot(1,3,2)
im = plt.imshow(img2, vmin=-1.5, vmax=1.5, cmap='coolwarm', aspect='auto')
plt.xticks(range(p_num), range(1, p_num+1))
plt.yticks(range(len(img2)), ylabels2)
plt.title("Epigentic marks")
plt.tight_layout()
plt.subplot(1,3,3)
cbar = plt.colorbar(im, cax = plt.gca())
cbar.ax.set_ylabel('Z-score', rotation=-90, va="bottom")
plt.tight_layout()
#plt.colorbar()
plt.show()
plt.close()
