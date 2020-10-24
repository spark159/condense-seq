import random

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

def GC_content(seq, percent=False):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    if percent:
        return (num/float(len(seq)))*100
    return (num/float(len(seq)))

def all_path (N, states='ATCG'):
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output

def hamming_dist (seq1, seq2):
    assert len(seq1) == len(seq2)
    dist = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist +=1
    return dist

# make DNA barcodes with minimized crosstalk 
def make_barcodes (alphabet, length, num, GCcontent=False):
    pool = []
    # check GC content
    if GCcontent:
        for barcode in all_path(length, states=alphabet):
            if GC_content(barcode) == GCcontent:
                pool.append(barcode)
    else:
        pool = all_path(length, states=alphabet)
    random.seed(0)
    random.shuffle(pool)

    threshold = length
    while threshold > 0:
        output = []
        for barcode in pool:
            i = 0
            while i < len(output): 
                if hamming_dist(barcode, output[i]) < threshold:
                    break
                i +=1
            if i == len(output):
                output.append(barcode)
        if len(output) >= num:
            print "Hamming dist: " + str(threshold)
            print "total barcode number: " + str(len(output))
            break
        threshold -=1
    return output, threshold

barcodes, threshold = make_barcodes('ATCG', 6, 100, GCcontent=0.5)
random.seed(0)
random.shuffle(barcodes)

barcode_list = []
for barcode in barcodes[:10]:
    print barcode
    barcode_list.append(barcode)

f = open("barcoded_primers.txt", 'w')
for i in range(len(barcode_list)):
    id = 'BC' + str(i+1) + '_0N0_p2'
    seq = barcode_list[i] + "TGGAGAATCCCGGTGCCG"
    print >> f, '>' + id
    print >> f, seq
f.close()

f = open("BC1_10.ref", 'w')
for i in range(len(barcode_list)):
    id = 'BC' + str(i+1)
    seq = "CAGGATGTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGATTCTCCA" + rev_comp(barcode_list[i])
    print >> f, '>' + id
    print >> f, seq
f.close()


    
