import os, sys, subprocess, re

def rev_comp (seq):
    rev_seq=seq[::-1]; new_seq=''
    for base in rev_seq:
        if base == 'A':
            new_seq += 'T'
        if base == 'T':
            new_seq += 'A'
        if base == 'C':
            new_seq += 'G'
        if base == 'G':
            new_seq += 'C'
    return new_seq

def map_to_genome (lib_fname, ref_fname):
    aligner_cmd=["bowtie2", '-x', ref_fname, '-U', lib_fname, '-f' ]
    align_proc = subprocess.Popen(aligner_cmd, stdout=subprocess.PIPE,stderr=open("/dev/null", 'w'))
    check = {'in':[], 'out':[]}
    pos_dic = {}
    for line in align_proc.stdout:
        if line.startswith('@'):
            continue
        #print line
    
        cols = line.strip().split()
        read_id, flag, ref_id, pos, mapQ, cigar_str = cols[:6]
        #read_id=":".join(read_id.split(':')[3:7])
        read_id = read_id.strip()
        flag, pos = int(flag), int(pos)
        pos -= 1

        # invalid: mapping failure
        if pos < 0:
            check['out'].append(read_id)
            pos_dic[read_id] = [None, None]
            continue
        if flag & 0x4 != 0:
            check['out'].append(read_id)
            pos_dic[read_id] = [None, None]
            continue
        
        read_seq, qual =cols[9:11]
        ref_id = ref_id.strip()
        if ref_id == '*':
            pos_dic[read_id] = [None, None]
            continue
        
        AS,NM,MD = None, None, None
        for i in range(11, len(cols)):
            col = cols[i]
            if col.startswith('AS'):
                AS = int(col[5:])
            elif col.startswith('NM'):
                NM = int(col[5:])
            elif col.startswith('MD'):
                MD = col[5:]

        # invalid: any mismatch or indel (edit distance > 0)
        if NM > 0:
            check['out'].append(read_id)
            pos_dic[read_id] = [None, None]
            continue

        # possible indel check
        M_count = 0
        cigar_str=re.split('(\d+)',cigar_str)[1:]
        for i in range(len(cigar_str)/2):
            s = cigar_str[2*i+1]
            num = int(cigar_str[2*i])
            if s == 'M':
                M_count += num

        if M_count != len(read_seq):
            check['out'].append(read_id)
            pos_dic[read_id] = [None, None]
            continue

        check['in'].append(read_id)
        pos_dic[read_id] = [ref_id, pos]
                
    return check, pos_dic

def read_Blib(fname):
    Blib = {}
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            id = line[1:]
            continue
        assert id not in Blib
        Blib[id] = line[23:len(line)-23]
    return Blib
Blib = read_Blib("data/Blib.ref")
#print Blib

# build small Blib sequence fasta file
f = open("data/Blibseq.fa", 'w')
for id in Blib:
    seq = Blib[id]
    print >> f, ">%s" % (id)
    print >> f, "%s" % (seq)
f.close()

Blib_check, Bpos_dic = map_to_genome("data/Blibseq.fa", "data/sacCer2")
#print len(Blib_check['in']), len(Blib_check['out'])

def read_Llib(fname):
    Llib = {}
    i = 0
    for line in open(fname):
        line = line.strip()
        try:
            seq, flex = line.split()
            Llib[i] = seq[2:len(seq)-2]
            i +=1
        except ValueError:
            continue    
    return Llib
Llib = read_Llib("data/DataReady.txt")
#print Llib

# build small Llib sequence fasta file
f = open("data/Llibseq.fa", 'w')
for id in Llib:
    seq = Llib[id]
    print >> f, ">%s" % (id)
    print >> f, "%s" % (seq)
f.close()

Llib_check, Lpos_dic = map_to_genome("data/Llibseq.fa", "data/sacCer2")
#print len(Llib_check['in']), len(Llib_check['out'])

def read_genome(fname):
    chr_dic = {}
    chr_name, sequence = "", ""
    for line in open(fname):
        if line.startswith(">"):
            if chr_name and sequence:
                chr_dic[chr_name] = sequence
            chr_name = line.strip().split()[0][1:]
            sequence = ""
        else:
            sequence += line.strip()
    if chr_name and sequence:
        chr_dic[chr_name] = sequence
    return chr_dic
ygenome = read_genome("data/sacCer2.fa")
#ygenome = read_genome("data/sacCer3.fa")
#print ygenome.keys()

def read_NPS(fname):
    NPS_info = {}
    for line in open(fname):
        line = line.strip()
        try:
            chr, pos, score, snr = line.split()
            pos, score, snr = int(pos) - 1, float(score), float(snr)
            if snr <= 0.0:
                continue
            if chr not in NPS_info:
                NPS_info[chr] = []
            NPS_info[chr].append([pos, score, snr])
            #NPS_info.append([chr, pos, score, snr])
        except ValueError:
            continue
    return NPS_info
NPS_info = read_NPS("data/nature11142-s2.txt")
#NPS_info = read_NPS("data/nature11142-s3.txt")
#print NPS_info.keys()

# check all pos of Blib are in NPS pos
check = {"in":0, "out":0}
for id, loc in Bpos_dic.items():
    chr, pos = loc
    Find = False
    for i in range(len(NPS_info[chr])):
        Npos, _, _ = NPS_info[chr][i]
        if pos + 50 == Npos:
            check['in'] +=1
            Find = True
            break
    if not Find:
        check['out'] +=1
print check

# check all pos of Llib are in NPS pos
check = {"in":0, "out":0}
for id, loc in Lpos_dic.items():
    chr, pos = loc
    Find = False
    for i in range(len(NPS_info[chr])):
        Npos, _, _ = NPS_info[chr][i]
        if pos + 50 == Npos or pos - 1 == Npos:
            check['in'] +=1
            Find = True
            break
    if not Find:
        check['out'] +=1
print check

# build small NPS sequence fasta file
f = open("data/NPSseq.fa", 'w')
num = 0
for chr in NPS_info:
    for i in range(len(NPS_info[chr])):
        Npos, score, snr = NPS_info[chr][i]
        Nst, Nen = Npos - 147/2, Npos + 147/2
        Nseq = ygenome[chr][Nst:Nen+1]
        Nseq = Nseq[23:len(Nseq)-23]
        print >> f, ">%d" % (num)
        print >> f, "%s" % (Nseq)
        num += 1
f.close()

def read_TSS(fname):
    TSS_info = {}
    for line in open(fname):
        line = line.strip()
        try:
            _, name, chr, strand, st, en, _, _, _, _, _, _ = line.split()
            st, en = int(st) - 1, int(en) - 1
            if chr not in TSS_info:
                TSS_info[chr] = {}
            assert name not in TSS_info[chr]
            TSS_info[chr][name] = [strand, st, en]
        except ValueError:
            continue
    return TSS_info
TSS_info = read_TSS("data/yeastgenes.txt")
#print TSS_info.keys()

def read_plusone (fname):
    num_dic = {1:"I", 2:"II", 3:"III", 4:"IV", 5:"V", 6:"VI", 7:"VII", 8:"VIII", 9:"IX", 10:"X", 11:"XI", 12:"XII", 13:"XIII", 14:"XIV", 15:"XV", 16:"XVI", 17:"XVII", 18:"XVIII"} 
    plusone = {}
    for line in open(fname):
        line = line.strip()
        try:
            pos, _, chr, _ = line.split()
            pos, chr = int(pos) - 1, int(chr)
            chr = "chr" + num_dic[chr]
            if chr not in plusone:
                plusone[chr] = []
            plusone[chr].append(pos)
        except ValueError:
            continue
    return plusone
plusone = read_plusone ("data/plusone.txt")
#print plusone.keys()

# check all plusone position in NPS
check = {'in':0, 'out':0}
for chr in plusone:
    for i in range(len(plusone[chr])):
        pos = plusone[chr][i]
        Find = False
        for j in range(len(NPS_info[chr])):
            Npos, _, _ = NPS_info[chr][j]
            if pos == Npos :
                check['in'] +=1
                Find = True
                break
        if not Find:
            check['out'] +=1
#print check

# check all plusone sequences are in the Blib
check = {'in':0, 'out':0}
for chr in plusone:
    for i in range(len(plusone[chr])):
        pos = plusone[chr][i]
        seq = ygenome[chr][pos-147/2:pos+147/2+1][23:len(seq)-23]
        if seq in Blib.values():
            check['in'] +=1
        else:
            check['out'] +=1
            #print seq
#print check

# check all plusone sequences are in the Llib
check = {'in':0, 'out':0}
for chr in plusone:
    for i in range(len(plusone[chr])):
        pos = plusone[chr][i]
        seq1 = ygenome[chr][pos-50:pos]
        seq2 = ygenome[chr][pos+1:pos+50+1]
        if seq1 in Llib.values() or seq2 in Llib.values():
            check['in'] +=1
        else:
            check['out'] +=1
            #print seq
#print check

