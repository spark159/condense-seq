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
from scipy import fft
import matplotlib as mpl

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

def stat_Markov(seq_list, NCPlen, order, sym=False):
    ntdic = {}
    for nt in all_path(order+1, 'ATCG'):
        ntdic[nt] = 0.0

    sample_num = len(seq_list)

    if sym:
        sample_num = 2*sample_num

    freq = [ copy.deepcopy(ntdic) for i in range(NCPlen - order) ]
    for seq in seq_list:
        if sym:
            rev_seq = rev_comp(seq)
        for i in range(len(seq) - order):
            nt = seq[i:i+1+order]
            freq[i][nt] += 1.0 / sample_num
            if sym:
                nt = rev_seq[i:i+1+order]
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
#path = ''
path = "/home/spark159/../../storage/"

#fname = 'H1_NCP_sp_chr1_extended_anot.cn'
#field = "work/2021_06_07_H1_sp_detail/H1-NCP-sp-8"

#fname = 'H1_NCP_HP1a_chr1_anot.cn'
#field = "/home/spark159/scratch4-tha4/sangwoo/2022_09_08_GM_sp_H1_HP1a_deep/H1-new-NCP-HP1a-3"

fname = 'H1_NCP_LKH_chr1_anot.cn'
field = "/home/spark159/scratch/2022_12_13_H1_LKH_deep/H1-NCP-LKH-3"

#fname = 'H1_NCP_Ki67_chr1_anot.cn'
#field = "/home/spark159/scratch/2022_10_28_H1_Ki67_deep/H1-NCP-Ki67-4"



ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+fname)
ID_score1 = name_ID_value[field]
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
    freq, sample_num, mean, std, stdz_freq = stat_Markov(seq_list, NCPlen, 1, sym=True)
    MMfreq_list.append(freq)
    

print "Sequence feature reading is done"

# plot AT vs GC periodicity

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

fig = plt.figure(figsize=(2.6, 2))
#plt.subplot(1,2,1)
ax1 = plt.gca()
ax2 = ax1.twinx()
for i in range(p_num):
    AT_sig, GC_sig = AT_sig_list[i], GC_sig_list[i]
    ax1.plot(AT_sig, color=cmap1(color_list[i]), lw=2)
    ax2.plot(GC_sig, color=cmap2(color_list[i]), lw=2)
ax1.set_ylabel('AA/AT/TA/TT freqeuncy', color='r', fontsize=8)
ax2.set_ylabel('CC/CG/GC/GG freqeuncy', color='b', fontsize=8)
ax1.set_xlabel("Super Helical Location", fontsize=8)
ax1.tick_params('y', colors='r', labelsize=8)
ax2.tick_params('y', colors='b', labelsize=8)
line_list = [147/2 - i*20 for i in range(4)] + [147/2 + i*20 for i in range(1, 4)]
for line in line_list:
    plt.axvline(x=line, color='k', linestyle='--', alpha=0.25)
plt.xticks([147/2 + 10*i for i in range(-7, 8, 2)], [str(10*i) for i in range(-7, 8, 2)], fontsize=5)
#plt.ylabel("Relative frequency")
#plt.legend()
#plt.ylim([0.22, 0.28])
#plt.savefig('ATGCperiod_' + "slide" + '.png')
plt.savefig('ATGCperiod_' + "slide" + '.svg', format='svg', bbox_inches='tight')
plt.title("Dinucleotide periodicity", fontsize=8)
#plt.show()
plt.close()

fig = plt.figure(figsize=(1,1))
plt.subplot(1,2,1)
img = [ [color, np.nan] for color in color_list ]
plt.imshow(img, cmap='OrRd', aspect=3)
img = [ [np.nan, color] for color in color_list ]
plt.imshow(img, cmap='GnBu', aspect=3)
plt.yticks(range(p_num), range(1, p_num+1), fontsize=8)
ax = plt.gca()
ax.yaxis.tick_right()
#plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
plt.xticks([])
#plt.savefig("cbar.png", bbox_inches='tight')
plt.savefig('cbar.svg', format='svg', bbox_inches='tight')
#plt.show()
plt.close()


# correlation analysis
def correlate (sig1, sig2, max_dist=sys.maxint, circular=False, clip=10):
    sig1 = sig1[clip:len(sig1)-clip]
    sig2 = sig2[clip:len(sig2)-clip]
    sig1 = np.asarray(sig1) - np.mean(sig1)
    sig2 = np.asarray(sig2) - np.mean(sig2)
    dist_products = {}
    for i in range(len(sig1)):
        for j in range(i, min(len(sig1), i+max_dist+1)):
            dist = j - i
            product = sig1[i]*sig2[j]
            if dist not in dist_products:
                dist_products[dist] = []
            dist_products[dist].append(product)
    corr_sig = [0.0]*(max(dist_products.keys())+1)
    for dist, products in dist_products.items():
        corr_sig[dist] = np.mean(products)
    return corr_sig

def FFT (sig, clip=10):
    sig = sig[clip:len(sig)-clip]
    N = len(sig)
    sig_ft = fft(sig)[1:N/2]
    periods = [float(N)/k for k in range(1, N/2)]
    amplts = np.abs(sig_ft)/float(N)
    phases = np.arctan2(sig_ft.imag, sig_ft.real) / np.pi # unit of pi
    shifts = np.asarray([(phases[k-1] * N) / (2*np.pi*k) for k in range(1, N/2)]) /np.pi # unit of pi
    return periods, amplts, phases

pair_corr_list = []
for i in range(p_num):
    pair_corr = {}
    AT_sig, GC_sig = AT_sig_list[i], GC_sig_list[i]
    pair_corr[('AT','AT')] = correlate(AT_sig, AT_sig, max_dist=49)
    pair_corr[('GC','GC')] = correlate(GC_sig, GC_sig, max_dist=49)
    pair_corr[('AT','GC')] = correlate(AT_sig, GC_sig, max_dist=49)
    pair_corr_list.append(pair_corr)
    #fig = plt.figure()
    #plt.plot(pair_corr[('AT','AT')])
    #plt.plot(pair_corr[('AT','GC')])
    #plt.title(str(i+1))
    #plt.show()
    plt.close()

pair_FFT_list = []
for i in range(p_num):
    pair_FFT = {}
    pair_corr = pair_corr_list[i]
    for pair in pair_corr:
        corr_sig = pair_corr[pair]
        periods, amplts, phases = FFT(corr_sig)
        #idx = sorted([(amplts[i], i) for i in range(len(amplts))], reverse=True)[0][1]
        idx = 4
        print periods[idx], amplts[idx], phases[idx]
        pair_FFT[pair] = (periods[idx], amplts[idx], abs(phases[idx]))
        #fig = plt.figure()
        #plt.title(str(i+1))
        #plt.plot(periods, phases)
        #plt.plot(periods, amplts)
        #plt.show()
        #plt.close()
    pair_FFT_list.append(pair_FFT)

def select_best (periods, amplts, min_period, max_period):
    amplt_idx = []
    for i in range(len(periods)):
        period = periods[i]
        amplt = amplts[i]
        if period < min_period:
            continue
        elif period > max_period:
            continue
        amplt_idx.append((amplt, i))
    return sorted(amplt_idx, reverse=True)[0][1]

AT_FFT_list, GC_FFT_list = [], []
for i in range(p_num):
    AT_sig, GC_sig = AT_sig_list[i], GC_sig_list[i]
    AT_sig = np.asarray(AT_sig) - np.mean(AT_sig)
    GC_sig = np.asarray(GC_sig) - np.mean(GC_sig)
    periods1, amplts1, phases1 = FFT(AT_sig)
    idx1 = select_best(periods1, amplts1, 8, 12)
    print periods1[idx1], amplts1[idx1], phases1[idx1]
    AT_FFT_list.append((periods1[idx1], amplts1[idx1], phases1[idx1]))
    periods2, amplts2, phases2 = FFT(GC_sig)
    idx2 = select_best(periods2, amplts2, 8, 12)
    print periods2[idx2], amplts2[idx2], phases2[idx2]
    GC_FFT_list.append((periods2[idx2], amplts2[idx2], phases2[idx2]))
    #fig = plt.figure()
    #plt.title(str(i+1))
    #plt.plot(periods1, amplts1)
    #plt.plot(periods2, amplts2)
    #plt.show()
    #plt.close()

AT_amplt_list, AT_phase_list, AT_mean_list = [], [], []
GC_amplt_list, GC_phase_list, GC_mean_list = [], [], []
for i in range(p_num):
    period1, amplt1, phase1 = AT_FFT_list[i]
    period2, amplt2, phase2 = GC_FFT_list[i]
    AT_amplt_list.append(amplt1)
    AT_phase_list.append(phase1*np.pi)
    AT_mean_list.append(np.mean(AT_sig_list[i]))
    GC_amplt_list.append(amplt2)
    GC_phase_list.append(phase2*np.pi)
    GC_mean_list.append(np.mean(GC_sig_list[i]))

def rescale (value, old_st, old_ed, new_st, new_ed):
    new_value = new_st + (new_ed - new_st)*float(value-old_st)/(old_ed-old_st)
    return new_value

min_mean, max_mean = np.min(AT_mean_list+GC_mean_list), np.max(AT_mean_list+GC_mean_list)
color_list = np.linspace(0.2, 1, num=p_num)
cmap = cm.get_cmap("jet")
cmap1, cmap2 = cm.get_cmap("OrRd"), cm.get_cmap("GnBu")

fig = plt.figure(figsize=(2.2, 2.2))
for i in range(p_num):
    plt.polar(AT_phase_list[i], np.log(AT_amplt_list[i]), '.', markersize=16, color=cmap1(color_list[i]), alpha=0.5)
    plt.text(AT_phase_list[i], np.log(AT_amplt_list[i]), '%d' % (int(i+1)), horizontalalignment='center', verticalalignment='center', color='k', size=7)
    plt.polar(GC_phase_list[i], np.log(GC_amplt_list[i]), '.', markersize=16, color=cmap2(color_list[i]), alpha=0.5)
    plt.text(GC_phase_list[i], np.log(GC_amplt_list[i]), '%d' % (int(i+1)), horizontalalignment='center', verticalalignment='center', color='k', size=7)

plt.text(1.17*np.pi, -6.1, 'AT-rich', horizontalalignment='center', verticalalignment='center', color='r', size=8)
plt.text(1.84*np.pi, -6.75, 'GC-rich', horizontalalignment='center', verticalalignment='center', color='b', size=8)
ax = plt.gca()
ax.set_rlabel_position(135)
rtick_list = [-7, -6.25, -6.5, -6.75, -6]
rlabel_list = ['$10^{' + str(tick) + '}$' for tick in rtick_list]
rlabel_list[0] = ''
rlabel_list[1] = ''
rlabel_list[3] = ''
ax.set_rticks(rtick_list) 
ax.set_yticklabels(rlabel_list, fontsize=5)
ax.tick_params(axis='both', which='major', labelsize=5, pad=-5)
ax.tick_params(axis='both', which='minor', labelsize=5, pad=-5)
#plt.savefig("ATGC_amplt_phase.png", bbox_inches='tight') 
plt.savefig("ATGC_amplt_phase.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()



fig = plt.figure(figsize=(2.6, 2.6))
for i in range(p_num):
    color1 = rescale(AT_mean_list[i], min_mean, max_mean, 0, 1)
    #plt.polar(AT_phase_list[i], np.log(AT_amplt_list[i]), '.', markersize=8, color=cmap(color1), alpha=0.5)
    plt.polar(AT_phase_list[i], np.log(AT_amplt_list[i]), '.', markersize=12, color=cmap1(color_list[i]), alpha=0.5)
    if i == 0:
        alignment = "top"
    else:
        alignment = "bottom"
    plt.text(AT_phase_list[i], np.log(AT_amplt_list[i]), '%d' % (int(i+1)), horizontalalignment='center', verticalalignment=alignment, color='k', size=5)
    color2 = rescale(GC_mean_list[i], min_mean, max_mean, 0, 1)
    #plt.polar(GC_phase_list[i], np.log(GC_amplt_list[i]), '.', markersize=8, color=cmap(color2), alpha=0.5)
    plt.polar(GC_phase_list[i], np.log(GC_amplt_list[i]), '.', markersize=12, color=cmap2(color_list[i]), alpha=0.5)
    if i == 3:
        alignment = "top"
    else:
        alignment = "bottom"
    plt.text(GC_phase_list[i], np.log(GC_amplt_list[i]), '%d' % (int(i+1)), horizontalalignment='center', verticalalignment=alignment, color='k', size=5)
plt.text(1.17*np.pi, -6.1, 'AT-rich', horizontalalignment='center', verticalalignment='center', color='r', size=8)
plt.text(1.84*np.pi, -6.75, 'GC-rich', horizontalalignment='center', verticalalignment='center', color='b', size=8)
ax = plt.gca()
ax.set_rlabel_position(135)
rtick_list = [-7, -6.25, -6.5, -6.75, -6]
rlabel_list = ['$10^{' + str(tick) + '}$' for tick in rtick_list]
rlabel_list[0] = ''
rlabel_list[1] = ''
rlabel_list[3] = ''
ax.set_rticks(rtick_list) 
ax.set_yticklabels(rlabel_list, fontsize=5)
ax.tick_params(axis='both', which='major', labelsize=5, pad=-5)
ax.tick_params(axis='both', which='minor', labelsize=5, pad=-5)
#ax.set_xlabel("Phase")
#ax.set_ylabel("Amplitude")
cmap = mpl.cm.jet
norm = mpl.colors.Normalize(vmin=min_mean, vmax=max_mean)
Map = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
Map.set_array([])
#cbar = plt.colorbar(Map, fraction=0.03, ticks=[min_mean, max_mean])
#cbar.ax.set_yticklabels([str(round(min_mean,2)), str(round(max_mean,2))], fontsize=5)
#cbar.ax.set_ylabel('Mean frequency', rotation=-90, va="bottom", labelpad=-2, fontsize=8)
plt.tight_layout()
#plt.savefig("ATGC_amplt_phase.png", bbox_inches='tight')
#plt.savefig("ATGC_amplt_phase.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()

# check all dinucleotides
din_sig_list = []
for i in range(p_num):
    pos_din_freq = MMfreq_list[i]
    din_sig = {}
    for i in range(len(pos_din_freq)):
        din_freq = pos_din_freq[i]
        for din, freq in din_freq.items():
            if din not in din_sig:
                din_sig[din] = [0.0]*len(pos_din_freq)
            din_sig[din][i] += freq
    din_sig_list.append(din_sig)

for i in range(p_num):
    din_sig = din_sig_list[i]
    fig = plt.figure()
    for din in din_sig:
        plt.plot(din_sig[din], label=din)
    plt.title(str(i+1))
    plt.legend()
    #plt.show()
    plt.close()
        

din_FFT_list = []
din_amplt_list, din_phase_list = [], []
for i in range(p_num):
    din_FFT = {}
    din_amplt, din_phase = {}, {}
    din_sig = din_sig_list[i]
    for din in din_sig:
        sig = np.asarray(din_sig[din])
        sig = sig - np.mean(sig)
        periods, amplts, phases = FFT(sig)
        idx = select_best(periods, amplts, 9, 11)
        print i, din, periods[idx], amplts[idx], phases[idx]
        din_FFT[din] = (periods[idx], amplts[idx], phases[idx])
        din_amplt[din] = amplts[idx]
        din_phase[din] = phases[idx]
    din_FFT_list.append(din_FFT_list)
    din_amplt_list.append(din_amplt)
    din_phase_list.append(din_phase)

for i in range(p_num):
    din_amplt = din_amplt_list[i]
    din_phase = din_phase_list[i]
    fig = plt.figure()
    for din in din_amplt:
        GC = GC_content(din)
        if GC < 0.5:
            color = 'r'
        elif GC > 0.5:
            color = 'b'
        else:
            color = 'g'
        plt.polar(din_phase[din]*np.pi, din_amplt[din], color+'o')
        plt.text(din_phase[din]*np.pi, din_amplt[din], din, horizontalalignment='center', verticalalignment='bottom', color=color)
    plt.title(str(i+1))
    #plt.show()
    plt.close()

din_amplt_trace, din_phase_trace = {}, {}
all_din = all_path(2)
for din in all_din:
    if din not in din_amplt_trace:
        din_amplt_trace[din] = []
    if din not in din_phase_trace:
        din_phase_trace[din] = []
    for i in range(0, p_num):
        amplt = din_amplt_list[i][din]
        phase = din_phase_list[i][din]
        din_amplt_trace[din].append(amplt)
        din_phase_trace[din].append(phase)

#temp = []
fig = plt.figure()
for din in din_amplt_trace:
    #if rev_comp(din) in temp:
    #    continue
    #temp.append(din)
    GC = GC_content(din)
    if GC < 0.5:
        #color = 'tab:red'
        color = 'red'
    elif GC > 0.5:
        #color = 'tab:blue'
        color = 'blue'
    else:
        color = 'green'
        #color = 'black'
    if color == 'green':
        continue
    amplt_trace = np.log(np.asarray(din_amplt_trace[din]))
    phase_trace = np.asarray(din_phase_trace[din])*np.pi
    p = plt.polar(phase_trace, amplt_trace, 'o--', markersize=12, alpha=0.5)
    if din == 'AT':
        valignment = 'top'
        halignment = 'center'
        theta = phase_trace[0] + 0.05
    else:
        valignment = 'bottom'
        halignment = 'left'
        theta = phase_trace[0]
    plt.text(theta, amplt_trace[0]+0.01, din, horizontalalignment=halignment, verticalalignment=valignment, color=color, size=15)
    #plt.text(phase_trace[0], amplt_trace[0]+0.01, din, horizontalalignment='center', verticalalignment='bottom', color=p[0].get_color(), size=15)
    for i in range(len(phase_trace)):
        plt.text(phase_trace[i], amplt_trace[i], str(i+1), horizontalalignment='center', verticalalignment='center', color='k', size=12)
#plt.legend()
ax = plt.gca()
ax.set_rlabel_position(270)
rtick_list = [-11, -10, -9, -8, -7]
rlabel_list = ['$10^{' + str(tick) + '}$' for tick in rtick_list]
rlabel_list[0] = ''
rlabel_list[2] = ''
rlabel_list[4] = ''
ax.set_rticks(rtick_list) 
ax.set_yticklabels(rlabel_list)
plt.title("Dinucleotide peridiocity change")
plt.tight_layout()
#plt.show()
plt.close() 
   
    
