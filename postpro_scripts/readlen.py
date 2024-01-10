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

def get_rlen (fnames,
              mm_cutoff,
              chr_list,
              out_fname):

    # gather whole file names
    data, label = [], []
    for fname in fnames:
        data.append(fname)
        label.append(fname.split('.')[0].split('/')[-1])

    # make genome read length dictionary
    chr_rlen = {}
        
    # open the sam file
    for i in range(len(data)):
        name = label[i]
        filename = data[i]
        print >> sys.stderr, "reading %s" % (filename)
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

            # record read length
            rlen = abs(tlen)
            
            # record mid-point of reads
            if tlen > 0: # left read
                end_pos = pos + tlen
            else: # right read
                end_pos = pos
                cigar_str=re.split('(\d+)',cigar_str)[1:]
                for i in range(len(cigar_str)/2):
                    s = cigar_str[2*i+1]
                    num = int(cigar_str[2*i])
                    if s == 'M' or s == 'D':
                        end_pos += num
                pos = end_pos + tlen
            mid = (pos + end_pos)/2
            
            if ref_id not in chr_rlen:
                chr_rlen[ref_id] = {}
            if mid not in chr_rlen[ref_id]:
                chr_rlen[ref_id][mid] = {}
            if name not in chr_rlen[ref_id][mid]:
                chr_rlen[ref_id][mid][name] = {"total":0, "count":0}
            chr_rlen[ref_id][mid][name]["total"] += rlen
            chr_rlen[ref_id][mid][name]["count"] += 1

    # summarize the output
    print >> sys.stderr, "writing read length file"
    
    f = open(out_fname + '_rlen.cn', 'w')
    s = 'SNP\tChromosome\tPhysicalPosition'
    for i in range(len(label)):
        s += '\t' + label[i]
    print >> f, s

    ID = 0
    for chr in sorted(chr_rlen.keys()):
        #start = min(chr_rlen[chr])
        #end = max(chr_rlen[chr]) + 1
        #for pos in xrange(start, end):
        for pos in sorted(chr_rlen[chr].keys()):
            s = str(ID) + "\t" + chr + "\t" + str(pos)
            write = False
            for name in label:
                try:
                    count = chr_rlen[chr][pos][name]['count']
                    if count <= 0:
                        rlen = '-'
                    else:
                        total = chr_rlen[chr][pos][name]['total']
                        rlen = int(round(total/float(count)))
                        write = True
                except:
                    rlen = '-'
                s += "\t" + str(rlen)
            if write:
                print >> f, s
                ID += 1

    f.close()

    

if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = ArgumentParser(description='Get read length along the genome')
    
    parser.add_argument(metavar='-f1',
                        dest="fnames",
                        type=str,
                        nargs='+',
                        help='SAM/Bam filenames (last file used as control)')
    parser.add_argument('-m',
                        dest="mm_cutoff",
                        type=int,
                        default=10,
                        help='mismatch cut-off in bp')
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

    chr_list = []
    if not args.chr_list:
        chr_list = None
    else:
        chr_list = sorted(args.chr_list)

    get_rlen (args.fnames,
              args.mm_cutoff,
              chr_list,
              args.out_fname)
