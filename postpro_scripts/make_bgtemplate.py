import sys
import copy
import random


def GC_content(seq):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    return (num/float(len(seq)))*100

def poly_score (seq, nts='ATCG', pos=False):
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
    #return max(num)
    score = 0
    for count in num:
        if count > 1:
            score +=count
    return score

def hamming_dist (seq1, seq2):
    assert len(seq1) == len(seq2)
    dist = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist +=1
    return dist

def all_path(N, states='ATCG'):
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output

def random_DNA (N, GC=False):
    while True:
        output = ""
        for i in range(N):
            output += random.choice(list('ATCG'))
        if GC == False or GC_content(output) == GC:
            break
    return output

def read_table (fname):
    BC_Histone, Histone_BC = {}, {}
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        cols = line.strip().split('\t')
        BC = cols[-2]
        Histone = cols[0]
        if BC == '-':
            continue
        assert BC not in BC_Histone
        assert Histone not in Histone_BC
        BC_Histone[BC] = Histone
        Histone_BC[Histone] = BC
    return BC_Histone, Histone_BC

def read_PTMtable (fname):
    BC_Histone, Histone_BC = {}, {}
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        cols = line.strip().split('\t')
        BC = cols[-1]
        Histone = cols[0]
        assert BC not in BC_Histone
        assert Histone not in Histone_BC
        BC_Histone[BC] = Histone
        Histone_BC[Histone] = BC
    return BC_Histone, Histone_BC

#BC_Histone, Histone_BC = read_table("oncohistoneTable.csv")
BC_Histone, Histone_BC = read_PTMtable("PTMhistoneTable.csv")

used_BCs = set(BC_Histone.keys())
all_BCs = set(all_path(6))
unused_BCs = all_BCs | used_BCs

# pick ones distant from used BCs
candidates = []
max_dist = 0 
for BC in list(unused_BCs):
    dist = min([hamming_dist(BC, bc) for bc in list(used_BCs)])
    if dist > max_dist:
        max_dist = copy.deepcopy(dist)
        candidates = [BC]
    elif dist == max_dist:
        candidates.append(BC)

# pick ones 50% GC content
new_candidates = []
for bc in candidates:
    if GC_content(bc) == 50:
        new_candidates.append(bc)

# sanity check
for BC in new_candidates:
    assert GC_content(BC) == 50
    for bc in list(used_BCs):
        assert hamming_dist(BC, bc) >= max_dist

final_candidates = []
min_polyscore = sys.maxint
for bc in new_candidates:
    polyscore = poly_score(bc)
    if polyscore < min_polyscore:
        min_polyscore = copy.deepcopy(polyscore)
        final_candidates = [bc]
    elif polyscore == min_polyscore:
        final_candidates.append(bc)



seq1 = "CTCTTTCCCTACACGACG"
seq2 = "CTGGAGAATCCCGGTG"
new_seq1 = "GATCGAACTCGCTACTGG"
new_seq2 = "TGCCTCACCAGCAACG"

        
#seq1 = "TACGGCGACCACCGAG"
#seq2 = "CATCAGAGACCTAGG"
#new_seq1 = "TCTGCTGCAGGGTGCG"
#new_seq2 = "ATTTGGTCGGGAACG"

#print random_DNA (len(seq1), GC_content(seq1))
#print random_DNA (len(seq2), GC_content(seq2))

print len(seq1), len(new_seq1)
print GC_content(seq1), GC_content(new_seq1)
print
print len(seq2), len(new_seq2)
print GC_content(seq2), GC_content(new_seq2)
