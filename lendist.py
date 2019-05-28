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
        label.append(fname.split('.')[0].split('/')[-1])

    # make genome read length dictionary
    output_list = []
        
    # open the sam file
    for i in range(len(data)):
        name = label[i]
        filename = data[i]
        print >> sys.stderr, "reading %s" % (filename)
        samtools_cmd = ["samtools",  "view", "-F 0x10", filename]
        samtools_proc = subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE, stderr=open("/dev/null", 'w'))

        chr_rlen = {}
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
            if ref_id not in chr_rlen:
                chr_rlen[ref_id] = {}
            if rlen not in chr_rlen[ref_id]:
                chr_rlen[ref_id][rlen] = 0
            chr_rlen[ref_id][rlen] += 1
            
        output_list.append(chr_rlen)

    # summarize the output
    print >> sys.stderr, "writing rlen file"

    for i in range(len(output_list)):
        chr_rlen = output_list[i]
        name = label[i]
        
        f = open(out_fname + '_' + name + '_rlen.txt', 'w')
        s = 'Chromosome\tReadLength\tCounts'
        print >> f, s

        for chr in sorted(chr_rlen.keys()):
            for rlen in sorted(chr_rlen[chr].keys()):
                count = chr_rlen[chr][rlen]
                s = chr + "\t" + str(rlen) + "\t" + str(count)
                print >> f, s

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

    parser = ArgumentParser(description='Get read length distribution')
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
