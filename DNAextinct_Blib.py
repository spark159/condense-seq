def count_cmp(a, b):
    if a[1] < b[1]:
        return -1
    else:
        return 1

def rev_cmp (seq):
    seq = seq.strip().replace(' ', '').upper()
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        assert nt in dic
        output+=dic[nt]
    return output[::-1]

def get_EC (seq):
    seq = seq.strip().replace(' ', '').upper()
    di_table = {'AA':27400, 'AC':21200, 'AG':25000, 'AT':22800, 'CA':21200, 'CC':14600, 'CG':18000, 'CT':15200, 'GA':25200, 'GC':17600, 'GG':21600, 'GT':20000, 'TA':23400, 'TC':16200, 'TG':19000, 'TT':16800}
    mono_table = {'A':15400, 'C':7400, 'G':11500, 'T':8700}
    total = 0.0
    for i in range(len(seq)-1):
        distep = seq[i:i+2]
        assert distep in di_table
        total += di_table[distep]
        if i > 0 :
            monostep = seq[i]
            assert monostep in mono_table
            total -= mono_table[monostep]
    return total

def get_fullEC (seq):
    seq = seq.strip().replace(' ', '').upper()    
    return get_EC(seq) + get_EC(rev_cmp(seq))

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

def get_seq_ncount (sort_fname, id_seq):
    seq_ncount = {}
    total = 0.0
    for line in open(sort_fname):
        cols = line.strip().split()
        if len(cols) != 3:
            continue
        _, id, _ = cols
        if not id.startswith('id_'):
            continue
        id = int(id[3:])
        seq = id_seq[id]
        if seq not in seq_ncount:
            seq_ncount[seq] = 0.0
        seq_ncount[seq] +=1
        total += 1.0
    for seq in seq_ncount:
        seq_ncount[seq] = seq_ncount[seq] / total
    return seq_ncount

def get_seq_EC (seq_ncount):
    seq_EC = []
    for seq in seq_ncount:
        seq_EC.append([seq, get_EC(seq)])
    return seq_EC

def get_libEC (seq_ncount):
    result = 0.0
    for seq, ncount in seq_ncount.items():
        result += get_EC(seq) * ncount
    return result

def get_fulllibEC (seq_ncount):
    result = 0.0
    for seq, ncount in seq_ncount.items():
        result += get_fullEC(seq) * ncount
    return result



id_seq = read_ref("Blib.ref")
seq_ncount = get_seq_ncount('sp1_S1_L001_R.sort', id_seq)
seqtoncount = [[seq,ncount] for seq, ncount in seq_ncount.items()]
seqtoncount = sorted(seqtoncount, cmp=count_cmp, reverse=True)
#for i in range(50):
#    print seqtoncount[i][0], seqtoncount[i][1]
#seq_EC = get_seq_EC(seq_ncount)
#seq_EC = sorted(seq_EC, cmp=count_cmp)
#for i in range(50):
#    print seq_EC[i][0], seq_EC[i][1]

print get_fulllibEC(seq_ncount)
