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

ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("/home/spark159/../../media/spark159/sw/sp_spd_tests_detail/hg19_chr1_NCP_ics_anot.cn")
ID_score1 = name_ID_value["data/sp_spd_tests_detail/sp7"]
ID_score2 = name_ID_value["data/sp_spd_tests_detail/sp8"]
ID_seq = name_ID_value['Sequence']
ID_AT = name_ID_value['ATcontent']
ID_CG = name_ID_value['CpGNumber']
ID_me = name_ID_value['meGCNumber']

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
plt.show()
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
        din_count = get_dincount(seq)
        for din in din_count:
            if din not in din_p_sigs:
                din_p_sigs[din] = [[] for k in range(p_num)]
            din_p_sigs[din][i].append(din_count[din])

print "Sequence feature reading is done"

fig = plt.figure()
plt.xlabel('Partitions')
plt.ylabel('%')
plt.title('AT content')
plt.boxplot(p_ATs, 0, "")
plt.savefig("ATcontent"+"_pbox.png")
plt.show()
plt.close()

for din in din_p_sigs:
    p_sigs = din_p_sigs[din]
    fig = plt.figure()
    plt.xlabel('Partitions')
    plt.ylabel('Counts')
    plt.title(din[0] + 'p' + din[1])
    plt.boxplot(p_sigs, 0, "")
    plt.savefig(din+"_pbox.png")
    plt.show()
    plt.close()

for poly in poly_p_sigs:
    p_sigs = poly_p_sigs[poly]
    fig = plt.figure()
    plt.xlabel('Partitions')
    plt.ylabel('Score (A.U.)')
    plt.title(poly)
    plt.boxplot(p_sigs, 0, "")
    plt.savefig('_'.join(poly.split('/'))+"_pbox.png")
    plt.show()
    plt.close()


# Compare epigenetic marks
ID_mefrac = {}
for ID in ID_CG:
    CG = ID_CG[ID]
    if CG <= 0:
        continue
    me = ID_me[ID]
    mefrac = float(me) / (2*CG)
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
            if not name.startswith('k'):
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
    if me == 'meCpG Denssity':
        plt.ylabel('Fraction')
    else:
        plt.ylabel('Counts')
    plt.title(me)
    plt.boxplot(p_sigs, 0, "")
    plt.savefig(me+"_pbox.png")
    plt.show()
    plt.close()

for chip in chip_p_sigs:
    p_sigs = chip_p_sigs[chip]
    fig = plt.figure()
    plt.xlabel('Partitions')
    plt.ylabel('Strength (A.U.)')
    plt.title(chip)
    plt.boxplot(p_sigs, 0, "")
    plt.savefig(chip+"_pbox.png")
    plt.show()
    plt.close()
