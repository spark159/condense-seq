import numpy as np
import matplotlib.pyplot as plt
import statis

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

def read_motif (fname, data_choice):
    ID_info = {}
    ID_seq = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.strip().split()
        if First:
            First = False
            continue
        sample, ID, chr, pos, strand, weight, seq = cols
        if sample != data_choice:
            continue
        assert ID not in ID_info
        assert ID not in ID_seq
        ID_info[ID] = (sample, chr, int(pos), strand, int(weight))
        ID_seq[ID] = seq.upper()
    return ID_info, ID_seq

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


ID_info1, ID_seq1 = read_motif("output_motif.txt", data_choice="data/sp_spd_tests_detail/sp1")
print "reading done"


weight_seq = []
for ID in ID_info1:
    info = ID_info1[ID]
    seq = ID_seq1[ID]
    weight = info[-1]
    weight_seq.append([weight, seq])

weight_seq = sorted(weight_seq, cmp=tuple_cmp, reverse=True)

seq_list1 = []
for ID in ID_seq1:
    seq = ID_seq1[ID]
    seq_list1.append(seq[:147])
    seq_list1.append(seq[1:148])
    seq_list1.append(seq[2:149])
    seq_list1.append(rev_cmp(seq)[:147])
    seq_list1.append(rev_cmp(seq)[1:148])
    seq_list1.append(rev_cmp(seq)[2:149])

print "combine done"
Afreq1, Gfreq1 = AG_freq(seq_list1)

fig, ax1 = plt.subplots()
ax1.plot(Afreq1, 'r')
ax1.set_ylabel('AA/AT/TA/TT freqeuncy', color='r')
ax1.tick_params('y', colors='r')
ax2 = ax1.twinx()
ax2.plot(Gfreq1, 'b')
ax2.set_ylabel('CC/CG/GC/GG freqeuncy', color='b')
ax2.tick_params('y', colors='b')
mid = 147/2
line_list = [mid - i*20 for i in range(4)] + [mid + i*20 for i in range(1, 4)]
xlabel_list = ['SHL' + str(-2*i) for i in range(4)] + ['SHL' + str(2*i) for i in range(1, 4)]

for line in line_list:
    plt.axvline(x=line, color='k', linestyle='--', alpha=0.5)

ax1.set_xticks(line_list)
ax1.set_xticklabels(xlabel_list)
plt.title("Dinucleotide frequency")
plt.savefig("DinFreq.png", bbox_inches='tight')
plt.show()
plt.close()
