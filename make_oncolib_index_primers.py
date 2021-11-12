def rev_comp (seq):
    rev_seq=seq[::-1]; new_seq=''
    for nt in rev_seq.upper():
        if nt == 'A':
            new_seq += 'T'
        elif nt == 'T':
            new_seq += 'A'
        elif nt == 'C':
            new_seq += 'G'
        elif nt == 'G':
            new_seq += 'C'
        else:
            new_seq += nt
    return new_seq

def read_index (fname):
    index_seq = {}
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            First = False
            continue
        index, seq = cols
        assert len(seq) == 6
        index = int(index)
        assert index not in index_seq
        index_seq[index] = seq
    return index_seq
index_seq = read_index("NEB_6N_single_index.csv")

assert len(set(index_seq.values())) == len(index_seq)

## Oncohistone library
#left = 'CAAGCAGAAGACGGCATACGAGAT'
#right = 'CCTAGGTCTCTGATGCTG'
#prefix = "Oncolib_index"


# PTM library
left = 'CAAGCAGAAGACGGCATACGAGAT'
right = 'CTGGAGAATCCCGGTG'
prefix = "PTMlib_index"


f = open(prefix+"_primers.csv", 'w')
print >> f, "%s\t%s\t%s" % ("Well Position", "Name", "Sequence")

pt = 1
Done = False
for letter in list('ABCDEFGH'):
    for number in range(1, 13):
        well_position = letter + str(number)
        name = prefix + ' ' + str(pt)
        try:
            seq = left + rev_comp(index_seq[pt]) + right
            assert len(seq) == len(left) + len(right) + 6
            print >> f, "%s\t%s\t%s" % (well_position, name, seq)
        except:
            Done = True
            break
        pt +=1
    if Done:
        break

f.close()


        
