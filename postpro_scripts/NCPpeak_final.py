import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import copy
import time

def Peak_calling (filenames,
                  min_len,
                  max_len,
                  mm_cutoff,
                  NCP_len,
                  chr_list,
                  out_fname):

    # collect filenames
    labels = [filename.rsplit('.',1)[0] for filename in filenames]
    # make a dictionary for NCP position and score
    chr_NPSs = [{} for i in range(len(labels))]

    # open the SAM/BAM files
    for filename, chr_NPS in zip(filenames, chr_NPSs):
        print >> sys.stderr, "reading %s" % (filename)
        samtools_cmd = ["samtools",  "view", "-F 0x10", filename]
        samtools_proc = subprocess.Popen(samtools_cmd,
                                         stdout=subprocess.PIPE,
                                         stderr=open("/dev/null", 'w'))

        for line in samtools_proc.stdout:
            # skip the header
            if line.startswith('@'):
                continue

            cols = line.strip().split()
            read_id, flag, ref_id, pos, mapQ, cigar_str = cols[:6]
            read_id=":".join(read_id.split(':')[3:7])
            
            tlen = int(cols[8])
            flag, pos = int(flag), int(pos)
            ref_id = ref_id.strip()
            pos-=1
            
            # non target chromosome
            if chr_list and ref_id not in chr_list:
                continue

            # invalid: mapping failure
            if pos < 0:
                #type = 'invalid:mutant'
                continue
            if flag & 0x4 != 0:
                #type = 'invalid:mutant'
                continue
            if ref_id == '*':
                continue

            # invalid: mate pairing failture
            if flag & 0x8 != 0:
                #type = 'invalid:mutant'
                continue
            if flag & 0x2 == 0:
                continue

            # invalid: ambiguous mapping
            #mapQ = float(mapQ)
            #if mapQ < 10:
            #    type = 'invalid:multimap'

            AS,NM,MD = None, None, None
            for i in range(11, len(cols)):
                col = cols[i]
                if col.startswith('AS'):
                    AS = int(col[5:])
                elif col.startswith('NM'):
                    NM = int(col[5:])
                elif col.startswith('MD'):
                    MD = col[5:]

            # invalid: too large edit distance 
            if NM > mm_cutoff:
                #type = 'invalid:mutant'
                continue

            # invalid: not right size of read (mono-NCP selection)
            if abs(tlen) < min_len or abs(tlen) > max_len:
                #type = 'invalid:mutant'
                continue

            # get NCP position and score
            if tlen > 0: # left-end read
                end_pos = pos + tlen
                if tlen % 2 != 0:
                    NCP_info = [[pos+tlen/2, 1]]
                else:
                    NCP_info = [[pos+tlen/2-1, 0.5],
                                [pos+tlen/2, 0.5]]
                    
            else: # skip right-end read to count the read-pair only once
                continue
            
            # record NCP position and score
            for Npos, Nscore in NCP_info:
                if ref_id not in chr_NPS:
                    chr_NPS[ref_id] = {}
                if Npos not in chr_NPS[ref_id]:
                    chr_NPS[ref_id][Npos] = 0
                chr_NPS[ref_id][Npos] += Nscore


    # need to implement
    """
    # Kernel smoothing of NCP positioning signal
    print >> sys.stderr, "peak calling:Kernel smoothing of NCP positions"
    def gauss (x):
        return np.exp(-(x**2)*0.5)/np.sqrt(2*np.pi)

    chr_KDE = {}
    for chr in sorted(chr_NCP.keys()):
        for NCPpos in sorted(chr_NCP[chr].keys()):
            try:
                score = chr_NCP[chr][NCPpos][label[-1]]
            except:
                continue
            if chr not in chr_KDE:
                chr_KDE[chr] = {}
            for k in range(NCPpos-3*bandwidth, NCPpos+3*bandwidth+1):
                value = score * gauss(float(k - NCPpos)/bandwidth)/bandwidth
                if k not in chr_KDE[chr]:
                    chr_KDE[chr][k] = 0.0
                chr_KDE[chr][k] += value
    del chr_NCP

    # find the local maximum of signals
    print >> sys.stderr, "find peaks of NCP positions"
    """
    
    # select NCPs with the highest score but minimizing overlapping
    print >> sys.stderr, "peak calling: filtering NCP positions"

    # sorting rule for tuple list
    def tuple_cmp (a,b):
        if a[0] <= b[0]:
            return -1
        else:
            return 1

    # find the index where the target would be inserted in right order
    def binary_insert (sortlist, target):
        st, ed = 0, len(sortlist)-1
        while st <= ed:
            mid = (st+ed) / 2
            if sortlist[mid] == target:
                return mid
            elif sortlist[mid] > target:
                ed = mid - 1 
            elif sortlist[mid] < target:
                st = mid + 1
        return st

    chr_peaks = [{} for i in range(len(labels))]
    for chr_NPS, label in zip(chr_NPSs, labels):
        for chr in sorted(chr_NPS.keys()):
            NPSlist = sorted([(Nscore, Npos) for Npos, Nscore in chr_NPS[chr].items()], cmp=tuple_cmp, reverse=True)
            selected = []
            
            
            

    chr_peak = {}
    for chr in sorted(chr_NCP.keys()):
        label_temp = {}
        for NCPpos in sorted(chr_NCP[chr].keys()):
            for label in labels:
                try:
                    score = chr_NCP[chr][NCPpos][label]
                    if label not in label_temp:
                        label_temp[label] = []
                    label_temp[label].append([score, NCPpos])
                except:
                    continue
        for label in label_temp:
            selected = []
            temp = sorted(label_temp[label], cmp=tuple_cmp, reverse=True)
            for k in range(len(temp)):
                score, NCPpos = temp[k]
                idx = binary_insert (selected, NCPpos)
                neighbor = selected[max(0,idx-1):min(idx+1,len(selected))]
                check = True
                for pos in neighbor:
                    if abs(pos - NCPpos) < NCP_len - overlap:
                        check = False
                        break
                if check:
                    selected.insert(idx, NCPpos)
                    if chr not in chr_peak:
                        chr_peak[chr] = {}
                    if NCPpos not in chr_peak[chr]:
                        chr_peak[chr][NCPpos] = {}
                    assert label not in chr_peak[chr][NCPpos]
                    chr_peak[chr][NCPpos][label] = score


    # summarize the output    
    print >> sys.stderr, "writing peak file"
    #print "writing peak file"
    
    f = open(out_fname + '_peak.cn', 'w')
    s = 'SNP\tChromosome\tPhysicalPosition'
    for i in range(len(label)):
        s += '\t' + label[i]
    print >> f, s

    ID = 0
    for chr in sorted(chr_peak.keys()):
        for NCPpos in sorted(chr_peak[chr].keys()):
            s = str(ID) + "\t" + chr + "\t" + str(NCPpos)
            for name in label:
                try:
                    score = chr_peak[chr][NCPpos][name]
                except:
                    score = 0
                s += "\t" + str(score)
            print >> f, s
            ID += 1

    f.close()

    print >> sys.stderr, "Done"
    

if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = ArgumentParser(description='Nucleosome peak calling')
    parser.add_argument(metavar='-f',
                        dest="filenames",
                        type=str,
                        nargs='+',
                        help='SAM/Bam filenames')
    parser.add_argument('-m',
                        dest="mm_cutoff",
                        type=int,
                        default=10,
                        help='mismatch cut-off in bp')
    parser.add_argument('--min',
                        dest="min_len",
                        type=int,
                        default=126,
                        help='minimum length for selection in bp (default:126bp)')
    parser.add_argument('--max',
                        dest="max_len",
                        type=int,
                        default=184,
                        help='maximum length for selection in bp (default:184bp)')
    parser.add_argument('--Nlen',
                        dest="NCP_len",
                        type=int,
                        default=147,
                        help='Mono-nucleosomal length in bp')
    parser.add_argument('--ovlap',
                        dest="overlap",
                        type=int,
                        default=40,
                        help='maximum allowed overlap between NCPS in bp')
    parser.add_argument('--chr',
                        dest="chr_list",
                        type=str,
                        nargs='+',
                        help='tagert chromosome list')
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()
