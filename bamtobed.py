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
               mm_cutoff,
               chr_list,
               out_fname):

    # gather whole file names
    data, label = [], []
    for fname in fnames:
        data.append(fname)
        label.append(fname.split('.')[0])

    # make genome counts dictionary
    chr_tag = {}
        
    # open the sam file
    for i in range(len(data)):
        name = label[i]
        #chr_cov = out_cov[i]
        filename = data[i]
        print >> sys.stderr, "reading %s" % (filename)
        #print "reading %s" % (filename)
        samtools_cmd = ["samtools",  "view", "-F 0x10", filename] # see only Forward strands
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

            # get NCP position and coverage range
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

            if ref_id not in chr_tag:
                chr_tag[ref_id] = []
            chr_tag[ref_id].append((pos, end_pos-1))
            
    # summarize the output
    def tuple_cmp (t1, t2):
        if t1[0] < t2[0]:
            return -1
        elif t1[0] > t2[0]:
            return 1
        else:
            if t1[1] < t2[1]:
                return -1
            elif t1[1] < t2[1]:
                return 1
            else:
                return 0
    
    print >> sys.stderr, "writing BED file"
    #print "writing NPS file"
    
    f = open(out_fname + '.bed', 'w')
    s = 'chromosome\tstart\tend'
    print >> f, s

    for chr in sorted(chr_tag.keys()):
        for tag in sorted(chr_tag[chr]):
            print >> f, "%s\t%d\t%d" % (chr, tag[0], tag[1])

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

    parser = ArgumentParser(description='Convert pair-end BAM file to BED like file')
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

    NCP_count (args.fnames,
               args.mm_cutoff,
               chr_list,
               args.out_fname
               )
