import os, sys, subprocess, re
import re
import random
import matplotlib.pyplot as plt
import numpy as np
from pyliftover import LiftOver

lo = LiftOver("data/liftover.chn")
num_dic = {1:"I", 2:"II", 3:"III", 4:"IV", 5:"V", 6:"VI", 7:"VII", 8:"VIII", 9:"IX", 10:"X", 11:"XI", 12:"XII", 13:"XIII", 14:"XIV", 15:"XV", 16:"XVI", 17:"XVII", 18:"XVIII"} 

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

def key_cmp(a, b):
    if a[0] > b[0]:
        return -1
    elif a[0] == b[0]:
        if a[1] >= b[1]:
            return -1
        else:
            return 1
    else:
        return 1

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

def get_new_pos (old_pos_dic, genome1, genome2):    
    new_pos_dic = {}
    for id, loc in old_pos_dic.items():
        chr, pos = loc
        seq = genome1[chr]
        st, en = max(0,pos-100), min(pos+100,len(seq))
        #st, en = max(0,pos-10000), min(pos+10000,len(seq))
        #st, en = pos, min(pos+10000,len(seq))
        seq = seq[st:en]
        text = genome2[chr]
        st, en = max(0,pos-10000), min(pos+10000,len(text))
        text = text[st:en]
        #st = 0
        temp = []
        for m in re.finditer(seq, text):
            new_pos = m.start() + st
            temp.append(new_pos)
        if len(temp) <= 0:
            continue
        elif len(temp) == 1:
            new_pos_dic[id] = new_pos
        else:
            new_pos_dic[id] = "$"
            print chr, pos
    return new_pos_dic        

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

#Blib_check, Bpos_dic = map_to_genome("data/Blibseq.fa", "data/sacCer2")
#print len(check['in']), len(check['out'])

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
#print len(Llib.keys())
#print Llib

# build small Llib sequence fasta file
f = open("data/Llibseq.fa", 'w')
for id in Llib:
    seq = Llib[id]
    print >> f, ">%s" % (id)
    print >> f, "%s" % (seq)
f.close()

#Llib_check, Lpos_dic = map_to_genome("data/Llibseq.fa", "data/sacCer2")
#print len(check['in']), len(check['out'])

def read_genome(fname, tricky=False):
    chr_dic = {}
    chr_name, sequence = "", ""
    for line in open(fname):
        if line.startswith(">"):
            if chr_name and sequence:
                chr_dic[chr_name] = sequence
            if not tricky:
                chr_name = line.strip().split()[0][1:]
            else:
                chr_name = line.strip().split()[5]
                #print chr_name
                chr_name = chr_name[1:len(chr_name)-1].split('=')[1]
                chr_name = 'chr' + chr_name
                #print chr_name
            sequence = ""
        else:
            sequence += line.strip()
    if chr_name and sequence:
        chr_dic[chr_name] = sequence
    return chr_dic
ygenome = read_genome("data/resacCer2.fa")
oldgenome = read_genome("data/rescAll.fa")
for i in range(1,17):
    name = "chr" + num_dic[i]
    if len(ygenome[name]) != len(oldgenome[name]):
        print name
#ygenome = read_genome("data/sacCer3.fa")
#print ygenome.keys(), len(ygenome.keys())

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
#print NPS_info.keys()

"""
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
"""
def read_TSS(fname, old=False):
    TSS_info = {}
    for line in open(fname):
        line = line.strip()
        try:
            if not old:
                _, name, chr, strand, st, en, _, _, _, _, _, _ = line.split()
                st, en = int(st) - 1, int(en) - 1
                if chr not in TSS_info:
                    TSS_info[chr] = {}
                assert name not in TSS_info[chr]
                TSS_info[chr][name] = [strand, st, en]
            else:
                _, chr, strand, st, en, _, name,  _, _, _ = line.split()
                st, en = int(st) - 1, int(en) - 1
                chr = 'chr' + num_dic[int(chr)]
                st = lo.convert_coordinate(chr, st)[0][1]
                en = lo.convert_coordinate(chr, en)[0][1]
                if chr not in TSS_info:
                    TSS_info[chr] = {}
                while name in TSS_info[chr]:
                    name += "p"
                assert name not in TSS_info[chr]
                TSS_info[chr][name] = [strand, st, en]
        except ValueError:
            continue
    return TSS_info
#TSS_info = read_TSS("data/yeastgenes.txt")
TSS_info = read_TSS("data/TSS.csv", old=True)
#print TSS_info.keys()
        
def read_plusone (fname):
    plusone = {}
    for line in open(fname):
        line = line.strip()
        try:
            pos, _, chr, _ = line.split()
            pos, chr = int(pos) - 1, int(chr)
            chr = "chr" + num_dic[chr]
            if chr not in plusone:
                plusone[chr] = set([])
            plusone[chr].add(pos)
        except ValueError:
            continue
    for chr in plusone:
        plusone[chr] = list(plusone[chr])
    return plusone
#plusone = read_plusone("data/plusone.txt")
plusone = read_plusone("data/plusone_ver2.txt")
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

def combine_all (genome, NPS, TSS):
    genome_NPS = {}
    for chr in TSS:
        if chr not in NPS:
            continue
        for gene in TSS[chr]:
            strand, st, en = TSS[chr][gene]
            for i in range(len(NPS[chr])):
                Npos, score, snr = NPS[chr][i]
                #Nst, Nen = Npos + 147/2, Npos - 147/2
                Nst, Nen = Npos, Npos
                assert st < en
                if Nst > st and Nen < en  :
                    Nseq = genome[chr][Npos-147/2:Npos+147/2+1]
                    assert len(Nseq) == 147
                    if chr not in genome_NPS:
                        genome_NPS[chr] = {}
                    if gene not in genome_NPS[chr]:
                        genome_NPS[chr][gene] = []
                    genome_NPS[chr][gene].append([strand, Npos, score, snr, Nseq])
                if Nst > en:
                    break
            if chr in genome_NPS and gene in genome_NPS[chr] and strand == '-':
                genome_NPS[chr][gene] = genome_NPS[chr][gene][::-1]
    return genome_NPS
11
genome_NPS = combine_all (ygenome, NPS_info, TSS_info)
#print genome_NPS

# check plusone pos are in genome_NPS +1
check = {"in":0,"out":0}
for chr in plusone:
    for i in range(len(plusone[chr])):
        pos = plusone[chr][i]
        Find = False
        for gene in genome_NPS[chr]:
            temp = []
            for j in range(min(len(genome_NPS[chr][gene]), len(genome_NPS[chr][gene]))):
                _, Npos, _, _, _ = genome_NPS[chr][gene][j]
                temp.append(Npos)
            if pos in temp:
                check['in'] +=1
                Find = True
                break
        if not Find:
            check['out'] +=1
            print chr, pos
print check

# gather all plus-one NPS in genome and remove redundancy
allplus = []
Nseq_list= []
for chr in genome_NPS:
    for gene in genome_NPS[chr]:
        strand, Npos, score, snr, Nseq = genome_NPS[chr][gene][0]
        data = [score, snr, chr, strand, Npos, Nseq]
        if Nseq not in Nseq_list:
            Nseq_list.append(Nseq)
            allplus.append(data)
        #if data not in allplus:
        #    allplus.append(data)

# check allplus is in plusone
check = {"in":0,"out":0}
for _,  _, chr, _, Npos, _ in allplus:
    Find = False
    for i in range(len(plusone[chr])):
        pos = plusone[chr][i]
        if pos == Npos:
            check['in'] +=1
            Find = True
            break
    if not Find:
        check['out'] +=1
print check

# check top 600 plus +1 are in Llib and pick top 500 included in Llib
f = open("plusonelib_table.txt", 'w')
print >> f, "%s\t%s\t%s\t%s\t%s\t%s" % ("chr", "strand", "pos", "score", "snr", "seq")
count = 0
check = {"in":0,"out":0}
allplus = sorted(allplus, cmp=key_cmp)[:600]
final_list = []
num_plus, num_minus = 0, 0
for score, snr, chr, strand, Npos, Nseq in allplus:
    Find1, Find2 = 0, 0
    Nseq1 = Nseq[len(Nseq)/2 - 50 :len(Nseq)/2]
    Nseq2 = Nseq[len(Nseq)/2 + 1 :len(Nseq)/2 + 51]
    assert len(Nseq1) == len(Nseq2) == 50
    for id in Llib:
        seq = Llib[id]
        if Find1 == 0 and seq == Nseq1:
            Find1 = 1
        if Find2 == 0 and seq == Nseq2:
            Find2 = 1
        if Find1 == 1 and Find2 == 1:
            check['in'] += 1
            if count < 500:
                print >> f, "%s\t%s\t%d\t%f\t%f\t%s" % (chr, strand, Npos, score, snr, Nseq)
                final_list.append([chr, strand, Npos, score, snr, Nseq])
                count += 1
                if strand == '+':
                    num_plus += 1
                else:
                    assert strand == '-'
                    num_minus += 1
            break
    if Find1 == 0 or Find2 == 0:
        check['out'] += 1
f.close()
print check
print num_plus, num_minus

# check final list are in plus-one
check = {"in":0,"out":0}
for chr, strand, Npos, score, snr, Nseq in final_list:
    Find = False
    for i in range(len(plusone[chr])):
        pos = plusone[chr][i]
        if pos == Npos:
            check['in'] +=1
            Find = True
            break
    if not Find:
        check['out'] +=1
        #print chr, strand, Npos, score, snr, Nseq
print check

# check final list are in Llib
check = {"in":0,"out":0}
for chr, strand, Npos, score, snr, Nseq in final_list:
    Nseq1 = Nseq[len(Nseq)/2 - 50 :len(Nseq)/2]
    Nseq2 = Nseq[len(Nseq)/2 + 1 :len(Nseq)/2 + 51]
    for sequence in [Nseq1, Nseq2]:
        Find = False
        for i in range(len(Llib)):
            seq = Llib[i]
            if Nseq1 == seq:
                Find = True
        if not Find:
            check['out'] += 1
            break
    if Find:
        check['in'] += 1
print check

# make plusonelib FASTA file
f = open("plusonelib.fa", 'w')
#left = "ATCCGACTGGCACCGGCAAGGTCGCTGTTCAATACATGCA"
#right = "GGGCGGCCGCGTATAGGGTCCATCACATAAGGGATGAACT"
#left = "ATCCGACTGGCACCGGCAAGGTCGCTGT TCGCCACATGCG"
#right = "GGGCGTCCTCGT ATAGGGTCCATCACATAAGGGATGAACT"

left = "CTGTTCGCCACATGCG"
right = "GGGCGTCCTCGTATAGG"

count = 0
First = True
for line in open("plusonelib_table.txt"):
    if First:
        First = False
        continue
    chr, strand, Npos, score, snr, Nseq = line.strip().split()
    seq = left + Nseq[1:-1] + right
    assert len(seq) == 178
    print >> f, ">%d" % (count)
    print >> f, seq
    count += 1
f.close()

def genomic_NCP_profile(NPS, TSS, win_size = 1000):
    profile = [0] * len(range(-win_size, win_size+1))
    for chr in TSS:
        if chr not in NPS:
            continue
        for gene in TSS[chr]:
            strand, st, en = TSS[chr][gene]
            for i in range(len(NPS[chr])):
                Npos, score, snr = NPS[chr][i]
                if Npos < st - win_size - 147/2:
                    continue
                if Npos > st + win_size + 147/2:
                    break
                for Ncc in range(Npos - 147/2, Npos + 147/2 + 1):
                    if Ncc >= st-win_size and Ncc <= st+win_size:
                        idx = Ncc - st + win_size
                        if strand == '-':
                            idx += 1
                            idx *= -1
                        profile[idx] += 1
    profile = [float(e)/sum(profile) for e in profile]
    return profile
                
#profile = genomic_NCP_profile(NPS_info, TSS_info)
#plt.plot(profile)
#plt.axvline(x=len(profile)/2, color='k', linestyle='--')
#plt.show()
                

def lib_check (lib_fname, ref_fname):
    aligner_cmd = ["bowtie2", '-x', ref_fname, '-U', lib_fname, '-f' ]
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
            continue
        if flag & 0x4 != 0:
            check['out'].append(read_id)
            continue
        
        read_seq, qual =cols[9:11]
        ref_id = ref_id.strip()
        if ref_id == '*':
            check['out'].append(read_id)
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

        # invalid: mismatch or indel (edit distance > 0)
        if NM > 0:
            check['out'].append(read_id)
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
            continue

        # record all aligned position
        if pos not in pos_dic:
            pos_dic[pos] = 0
        pos_dic[pos] +=1

        # aligned position check
        #if pos != 23:
        #    check['out'].append(read_id)
        #    continue
            
        check['in'].append(read_id)
                
    return check, pos_dic
    
# check all Blib sequences in NPS seq
#Blib_check, pos_dic = lib_check("data/Blibseq.fa", "data/NPSseq")
#print pos_dic
#print len(Blib_check['in']), len(Blib_check['out'])

# check all Llib sequences in NPS seq
#Llib_check, pos_dic = lib_check("data/Llibseq.fa", "data/NPSseq")
#print pos_dic
#print len(Llib_check['in']), len(Llib_check['out'])
