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

def read(fname, data_choice='None'):
    ID_profile = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        if First:
            total = [0.0] * len(cols[6:])
            First = False
            continue
        sample, ID = cols[:2]
        if sample != data_choice:
            continue
        profile = []
        for i in range(len(cols[6:])):
            value = float(cols[6:][i])
            profile.append(value)
        #print ID
        #assert ID not in ID_profile
        profile = np.asarray(profile)#/sum(profile)
        ID_profile[ID] = profile
        for i in range(len(profile)):
            total[i] += profile[i]
    
    return np.asarray([total[i]/len(ID_profile) for i in range(len(total))]), ID_profile

#total_list = [4027035549.0, 4106574224.0, 3867424201.0]
#mean_profile = read("data/sp1_chr1_profile.txt")
#mean_profile = read("output_profile.txt")
#mean_profile1, ID_profile1 = read("score_TSS_profile.txt", data_choice="work/condense_seq/sp1_hg19_chr1")#/total_list[-1]
#mean_profile2, ID_profile2 = read("score_ETS_profile.txt", data_choice="work/condense_seq/sp9_hg19_chr1")#/total_list[0]
#mean_profile3, ID_profile3 = read("score_ETS_profile.txt", data_choice="work/condense_seq/sp10_hg19_chr1")#/total_list[1]

mean_profile1, ID_profile1 = read("test_profile.txt", data_choice="work/condense_seq/sp10_hg19_chr1")
mean_profile2, ID_profile2 = read("test_profile.txt", data_choice="meGCNumber")

mean_profile1 = statis.moving_average(mean_profile1, 100)
mean_profile2 = statis.moving_average(mean_profile2, 100)

fig = plt.figure()
X = [ i - 1000 for i in range(len(mean_profile1))]
#plt.plot(X, mean_profile1/sum(mean_profile1), label='sp1')
plt.plot(X, mean_profile1, label='Condensability')
#plt.plot(X, mean_profile2/sum(mean_profile2), label='Condensibility')
plt.plot(X, mean_profile2, label='me G orC')
plt.legend()
plt.xlabel("Distance from TSS (bp)")
plt.ylabel("Condensibility (A.U.)")
#plt.ylabel("Condensibility (A.U.)")
#plt.savefig("condensibility_TTS_profile.png")
#plt.savefig("condensibility_profile.png")
#xplt.savefig("fake.norm.png")
plt.show()

"""
SumProfile = []
for ID, profile in ID_profile1.items():
    total = sum(profile[500:1500])
    SumProfile.append([total, ID])
SumProfile = sorted(SumProfile, cmp=tuple_cmp, reverse=True)

img = []
for total, ID in SumProfile:
    img.append(np.log2(ID_profile1[ID]))

fig = plt.figure()
plt.imshow(img, interpolation='none', aspect='auto', cmap='jet')
plt.show()


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


#ID_info, ID_seq = read_motif("output_motif.txt")
#ID_info, ID_seq = read_motif("GS_motif.txt")
ID_info1, ID_seq1 = read_motif("data/hg19_chr1_peak_motif.txt", data_choice="work/condense_seq/sp1_hg19")
ID_info2, ID_seq2 = read_motif("data/hg19_chr1_peak_motif.txt", data_choice="work/condense_seq/sp9_hg19")
ID_info3, ID_seq3 = read_motif("data/hg19_chr1_peak_motif.txt", data_choice="work/condense_seq/sp10_hg19")
print "done"


weight_seq = []
for ID in ID_info1:
    info = ID_info[ID]
    seq = ID_seq[ID]
    weight = info[-1]
    weight_seq.append([weight, seq])

weight_seq = sorted(weight_seq, cmp=tuple_cmp, reverse=True)


seq_list1 = []
seq_list2 = []
seq_list3 = []
for ID in ID_seq1:
    seq = ID_seq1[ID]
    seq_list1.append(seq[:147])
    seq_list1.append(seq[1:148])
    seq_list1.append(seq[2:149])
    seq_list1.append(rev_cmp(seq)[:147])
    seq_list1.append(rev_cmp(seq)[1:148])
    seq_list1.append(rev_cmp(seq)[2:149])
for ID in ID_seq2:
    seq = ID_seq2[ID]
    seq_list2.append(seq[:147])
    seq_list2.append(seq[1:148])
    seq_list2.append(seq[2:149])
    seq_list2.append(rev_cmp(seq)[:147])
    seq_list2.append(rev_cmp(seq)[1:148])
    seq_list2.append(rev_cmp(seq)[2:149])
for ID in ID_seq3:
    seq = ID_seq3[ID]
    seq_list3.append(seq[:147])
    seq_list3.append(seq[1:148])
    seq_list3.append(seq[2:149])
    seq_list3.append(rev_cmp(seq)[:147])
    seq_list3.append(rev_cmp(seq)[1:148])
    seq_list3.append(rev_cmp(seq)[2:149])

print "combine done"
Afreq1, Gfreq1 = AG_freq(seq_list1)
Afreq2, Gfreq2 = AG_freq(seq_list2)
Afreq3, Gfreq3 = AG_freq(seq_list3)

fig = plt.figure()
plt.plot(Afreq1, label='sp1')
plt.plot(Afreq2, label='sp9')
plt.plot(Afreq3, label='sp10')
#mid = len(seq_list[0])/2
#line_list = [mid - i*10 for i in range(7)] + [mid + i*10 for i in range(1, 7)]
#for line in line_list:
#    plt.axvline(x=line, color='k', linestyle='--')
#plt.plot(Gfreq, label='GC')
plt.ylabel('AA/AT/TA/TT frequency')
plt.legend()
plt.show()
"""
