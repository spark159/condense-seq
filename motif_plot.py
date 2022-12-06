import numpy as np
import matplotlib.pyplot as plt
import statis
import matplotlib.cm as cm

def tuple_cmp(a, b):
    if a[0] < b[0]:
        return -1
    elif a[0] > b[0]:
        return 1
    else:
        return 0

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

def read_motif (fname, data_choice=None):
    name_ID_info = {}
    name_ID_seq = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.strip().split()
        if First:
            First = False
            continue
        name, ID, chr, pos, strand, weight, seq = cols

        if data_choice != None and name != data_choice:
            continue
        
        if name not in name_ID_info:
            name_ID_info[name] = {}
        assert ID not in name_ID_info[name]
        name_ID_info[name][ID] = (chr, int(pos), strand, int(weight))

        if name not in name_ID_seq:
            name_ID_seq[name] = {}
        assert ID not in name_ID_seq[name] 
        name_ID_seq[name][ID] = seq.upper()

    if data_choice != None:
        return name_ID_info[data_choice], name_ID_seq[data_choice]
        
    return name_ID_info, name_ID_seq


def AG_freq (NCP_seq_list):
    nucleosome_dna_len = len(NCP_seq_list[0])
    Afreq=np.zeros(nucleosome_dna_len - 1); Gfreq=np.zeros(nucleosome_dna_len - 1)
    for seq in NCP_seq_list:
        for i in range(len(seq)-1):
            dint = seq[i:i+2]
            if dint in ['AA','AT','TA','TT']:
                Afreq[i] += 1.0
            elif dint in ['CC','CG','GC','GG']:
                Gfreq[i] += 1.0
    return Afreq / len(NCP_seq_list), Gfreq / len(NCP_seq_list)

#names = ['work/2021_06_07_H1_sp_detail/H1-NCP-sp-0', 'work/2021_06_07_H1_sp_detail/H1-NCP-sp-4', 'work/2021_06_07_H1_sp_detail/H1-NCP-sp-8']

#names = ["/home/spark159/scratch4-tha4/sangwoo/2022_09_08_GM_sp_H1_HP1a_deep/GM-NCP-sp-0", "/home/spark159/scratch4-tha4/sangwoo/2022_09_08_GM_sp_H1_HP1a_deep/GM-NCP-sp-4", "/home/spark159/scratch4-tha4/sangwoo/2022_09_08_GM_sp_H1_HP1a_deep/GM-NCP-sp-8"]

#names = ["/home/spark159/scratch4-tha4/sangwoo/2022_09_08_GM_sp_H1_HP1a_deep/H1-DNA-HP1a-0", "/home/spark159/scratch4-tha4/sangwoo/2022_09_08_GM_sp_H1_HP1a_deep/H1-DNA-HP1a-3"]


# mouse CD8 T cell in detail
#path = "/home/spark159/../../media/spark159/sw/"
#names = ["/home/spark159/scratch4-tha4/sangwoo/MouseCD8Tcell_detail/mCD8T-WT-NCP-sp-0",
#         "/home/spark159/scratch4-tha4/sangwoo/MouseCD8Tcell_detail/mCD8T-WT-NCP-sp-4",
#         "/home/spark159/scratch4-tha4/sangwoo/MouseCD8Tcell_detail/mCD8T-WT-NCP-sp-8"]
#fname = "mCD8T_WT-NCP_sp_chr1_motif.txt"

#path = "/home/spark159/../../media/spark159/sw/"
#names = ["/home/spark159/scratch4-tha4/sangwoo/MouseCD8Tcell_detail/mCD8T-inht-NCP-sp-0",
#         "/home/spark159/scratch4-tha4/sangwoo/MouseCD8Tcell_detail/mCD8T-inht-NCP-sp-4",
#         "/home/spark159/scratch4-tha4/sangwoo/MouseCD8Tcell_detail/mCD8T-inht-NCP-sp-8"]
#fname = "mCD8T_inht-NCP_sp_chr1_motif.txt"

path = "/home/spark159/../../media/spark159/sw/"
names = ["/home/spark159/scratch/2022_11_08_mCD8T_KO_detail/mCD8T-KO-NCP-sp-0",
         "/home/spark159/scratch/2022_11_08_mCD8T_KO_detail/mCD8T-KO-NCP-sp-4",
         "/home/spark159/scratch/2022_11_08_mCD8T_KO_detail/mCD8T-KO-NCP-sp-8"]
fname = "mCD8T_KO-NCP_sp_chr1_motif.txt"

name_Afreq, name_Gfreq = {}, {}
for name in names:
    ID_info, ID_seq = read_motif(path + fname, data_choice=name)
    print name, "reading done"
    
    weight_seq = []
    for ID in ID_info:
        info = ID_info[ID]
        seq = ID_seq[ID]
        weight = info[-1]
        weight_seq.append([weight, seq])

    weight_seq = sorted(weight_seq, cmp=tuple_cmp, reverse=True)

    # select top 90 %
    weight_seq = weight_seq[:int(len(weight_seq)*0.9)]

    seq_list = []
    for weight, seq in weight_seq:
        seq_list.append(seq[:147])
        seq_list.append(seq[1:148])
        seq_list.append(seq[2:149])
        seq_list.append(rev_cmp(seq)[:147])
        seq_list.append(rev_cmp(seq)[1:148])
        seq_list.append(rev_cmp(seq)[2:149])
    del weight_seq
    
    print name, "pooling sequences done"
    
    Afreq, Gfreq = AG_freq(seq_list)
    del seq_list

    name_Afreq[name] = Afreq
    name_Gfreq[name] = Gfreq

    print name, "frequency calculation done"


# plot dinucleotide frequency
color_list = np.linspace(0.3, 1, num=len(names))
cmap1, cmap2 = cm.get_cmap("OrRd"), cm.get_cmap("GnBu")

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

for i in range(len(names)):
    name = names[i]
    Afreq, Gfreq = name_Afreq[name], name_Gfreq[name]
    ax1.plot(Afreq, color=cmap1(color_list[i]), label = name.split('/')[-1])
    ax2.plot(Gfreq, color=cmap2(color_list[i]), label = name.split('/')[-1])

ax1.set_ylabel('AA/AT/TA/TT freqeuncy', color='r')
ax1.tick_params('y', colors='r')
ax2.set_ylabel('CC/CG/GC/GG freqeuncy', color='b')
ax2.tick_params('y', colors='b')

mid = 147/2
line_list = [mid - i*20 for i in range(4)] + [mid + i*20 for i in range(1, 4)]
xlabel_list = ['SHL' + str(-2*i) for i in range(4)] + ['SHL' + str(2*i) for i in range(1, 4)]

for line in line_list:
    plt.axvline(x=line, color='k', linestyle='--', alpha=0.5)

ax1.set_xticks(line_list)
ax1.set_xticklabels(xlabel_list)
ax1.legend(loc='upper left')
ax2.legend(loc='lower right')
plt.title("Dinucleotide frequency")
plt.savefig("DinFreq.png", bbox_inches='tight')
#plt.show()
plt.close()
