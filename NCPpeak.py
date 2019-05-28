import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import copy
import time

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

def NCP_count (fnames,
               genome_size,
               min_len,
               max_len,
               mm_cutoff,
               NCP_len,
               overlap,
               cov_choice,
               skip_zero,
               chr_list,
               out_fname):

    # gather whole file names
    data, label = [], []
    for fname in fnames:
        data.append(fname)
        label.append(fname.split('.')[0])

    # make genome counts dictionary
    chr_NCP = {}
    
    if cov_choice:
        chr_profile = {}
    
    # open the sam file
    for i in range(len(data)):
        name = label[i]
        #chr_cov = out_cov[i]
        filename = data[i]
        print >> sys.stderr, "reading %s" % (filename)
        #print "reading %s" % (filename)
        samtools_cmd = ["samtools",  "view", "-F 0x10", filename]
        samtools_proc = subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE, stderr=open("/dev/null", 'w'))

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

            # get NCP position and coverage range
            if tlen > 0: # left read
                end_pos = pos + tlen
                if tlen % 2 != 0:
                    NCPscore = [[pos+tlen/2, 1]]
                else:
                    NCPscore = [[pos+tlen/2-1, 0.5], [pos+tlen/2, 0.5]]
            else: # right read
                end_pos = pos
                cigar_str=re.split('(\d+)',cigar_str)[1:]
                for i in range(len(cigar_str)/2):
                    s = cigar_str[2*i+1]
                    num = int(cigar_str[2*i])
                    if s == 'M' or s == 'D':
                        end_pos += num
                pos = end_pos + tlen
                if tlen % 2 != 0:
                    NCPscore = [[end_pos+tlen/2-1, 1]]
                else:
                    NCPscore = [[end_pos+tlen/2-1, 0.5], [end_pos+tlen/2, 0.5]]
            
            # record coverage range
            if cov_choice:           
                if ref_id not in chr_profile:
                    chr_profile[ref_id] = {}
                if pos not in chr_profile[ref_id]:
                    chr_profile[ref_id][pos] = {}
                if name not in chr_profile[ref_id][pos]:
                    chr_profile[ref_id][pos][name] = 0
                if end_pos not in chr_profile[ref_id]:
                    chr_profile[ref_id][end_pos] = {}
                if name not in chr_profile[ref_id][end_pos]:
                    chr_profile[ref_id][end_pos][name] = 0
                chr_profile[ref_id][pos][name] += 1
                chr_profile[ref_id][end_pos][name] -= 1

            # record NCP position
            for NCPpos, score in NCPscore:
                if ref_id not in chr_NCP:
                    chr_NCP[ref_id] = {}
                if NCPpos not in chr_NCP[ref_id]:
                    chr_NCP[ref_id][NCPpos] = {}
                if name not in chr_NCP[ref_id][NCPpos]:
                    chr_NCP[ref_id][NCPpos][name] = 0
                chr_NCP[ref_id][NCPpos][name] += score

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
    
    # select NCP positions with highest score and minimized overlapping    
    print >> sys.stderr, "peak calling: filtering NCP positions"
    #print "filtering NCP positions"
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

    chr_peak = {}
    for chr in sorted(chr_NCP.keys()):
        name_temp = {}
        for NCPpos in sorted(chr_NCP[chr].keys()):
            for name in label:
                try:
                    score = chr_NCP[chr][NCPpos][name]
                    if name not in name_temp:
                        name_temp[name] = []
                    name_temp[name].append([score, NCPpos])
                except:
                    continue
        for name in name_temp:
            selected = []
            temp = sorted(name_temp[name], cmp=tuple_cmp, reverse=True)
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
                    assert name not in chr_peak[chr][NCPpos]
                    chr_peak[chr][NCPpos][name] = score


    # summarize the output
    print >> sys.stderr, "writing NPS file"
    #print "writing NPS file"
    
    f = open(out_fname + '_NPS.cn', 'w')
    s = 'SNP\tChromosome\tPhysicalPosition'
    for i in range(len(label)):
        s += '\t' + label[i]
    print >> f, s

    ID = 0
    for chr in sorted(chr_NCP.keys()):
        for NCPpos in sorted(chr_NCP[chr].keys()):
            s = str(ID) + "\t" + chr + "\t" + str(NCPpos)
            for name in label:
                try:
                    score = chr_NCP[chr][NCPpos][name]
                except:
                    score = 0
                s += "\t" + str(score)
            print >> f, s
            ID += 1
        
    f.close()

    
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

    if cov_choice:
        print >> sys.stderr, "writing coverage file"
        #print "writing coverage file"
    
        f = open(out_fname + '_cov.cn', 'w')
        s = 'SNP\tChromosome\tPhysicalPosition'
        for i in range(len(label)):
            s += '\t' + label[i]
        print >> f, s
            
        ID = 0
        for chr in sorted(chr_profile.keys()):
            previous = [0]*len(label)
            if skip_zero:
                start = min(chr_profile[chr])
                end = max(chr_profile[chr]) + 1
            else:
                start = 0
                end = genome_size[chr]
            for i in xrange(start, end):
                s = str(ID) + "\t" + chr + "\t" + str(i)
                for j in range(len(label)):
                    name = label[j]
                    past = previous[j]
                    try:
                        current = past + chr_profile[chr][i][name]
                    except:
                        current = past + 0
                    s += "\t" + str(current)
                    previous[j] = current
                if skip_zero and sum(previous) == 0:
                    continue
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

    parser = ArgumentParser(description='Calling NCP peaks and scores')
    parser.add_argument(metavar='-f1',
                        dest="fnames",
                        type=str,
                        nargs='+',
                        help='SAM/Bam filenames (last file used as control)')
    parser.add_argument(metavar='-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence filename')
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
    parser.add_argument('--cov',
                        dest="cov_choice",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='calculate coverage option')
    parser.add_argument('--skip',
                        dest="skip_zero",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='skip the zero coverage positions')
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

    # get length for each chromosome
    genome_size = {}

    if args.cov_choice:
        if not args.ref_fname:
            print >> sys.stderr, "Error: there is no reference file input."
            sys.exit(1)
            
        for line in open(args.ref_fname):
            line = line.strip()
            if line.startswith('>'):
                key = line[1:]
                assert key not in genome_size
                genome_size[key] = 0
                continue
            genome_size[key] += len(line)

    
    chr_list = []
    if not args.chr_list:
        chr_list = genome_size.keys()
    else:
        chr_list = sorted(args.chr_list)

    NCP_count (args.fnames,
               genome_size,
               args.min_len,
               args.max_len,
               args.mm_cutoff,
               args.NCP_len,
               args.overlap,
               args.cov_choice,
               args.skip_zero,
               chr_list,
               args.out_fname
               )
