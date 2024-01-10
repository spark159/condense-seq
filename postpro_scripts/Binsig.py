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

def Bin_sig (cov_fname,
             genome_size,
             bin_size,
             bin_step,
             skip_zero,
             chr_list,
             sample_choice,
             out_fname):

    # mark start and end of each bin along genome
    print >> sys.stderr, "binning the genome"
    chr_st, chr_ed = {}, {}
    for chr in chr_list:
        if chr not in chr_st:
            chr_st[chr] = []
            chr_ed[chr] = []
        st = 0
        while st < genome_size[chr]:
            ed = min(st+bin_size-1, genome_size[chr]-1)
            chr_st[chr].append(st)
            chr_ed[chr].append(ed)
            st += bin_step
    assert len(chr_st) == len(chr_ed)

    # reading signal file and get binned siganl
    print >> sys.stderr, "reading signal file"
    chr_Bsig = {}  
    First = True
    order = 2
    for line in open(cov_fname):
        cols = line.strip().split()
        if First:
            label = cols[3:]
            if sample_choice == 'all':
                target = label
            else:
                target = [cols[-1]] # control only
            First = False
            continue
        _, chr, pos = cols[0], cols[1], int(cols[2])
        if chr not in chr_st:
            continue
        if math.log10(pos+1) > order:
            print >> sys.stderr, chr + " pos " + str(pos) +" is reading"
            order += 1
        idx = pos / bin_step 
        #offset = pos % bin_step
        st, ed = chr_st[chr][idx], chr_ed[chr][idx]
        while pos >= st and pos <= ed:
            if chr not in chr_Bsig:
                chr_Bsig[chr] = {}
            if (st, ed) not in chr_Bsig[chr]:
                chr_Bsig[chr][(st, ed)] = {}
            counts = cols[3:]
            for i in range(len(counts)):
                name = label[i]
                if name not in target:
                    continue
                count = round(float(counts[i]), 3)
                if name not in chr_Bsig[chr][(st, ed)]:
                    chr_Bsig[chr][(st, ed)][name] = 0
                chr_Bsig[chr][(st, ed)][name] += count
            idx -= 1
            if idx < 0:
                break
            st, ed = chr_st[chr][idx], chr_ed[chr][idx]
    
    # write bin signal file
    print >> sys.stderr, "writing bin signal file"
    
    f = open(out_fname + '_Bsig.cn', 'w')
    s = 'SNP\tChromosome\tPhysicalPosition'
    for i in range(len(target)):
        s += '\t' + target[i]
    print >> f, s

    ID = 0
    for chr in sorted(chr_st.keys()):
        for i in xrange(len(chr_st[chr])):
            st, ed = chr_st[chr][i], chr_ed[chr][i]
            if skip_zero:
                try: 
                    total = sum(chr_Bsig[chr][(st, ed)].values())
                except:
                    continue
                if total <= 0:
                    continue
            BINpos = (st + ed)/2
            s = str(ID) + "\t" + chr + "\t" + str(BINpos)
            for name in target:
                try:
                    score = chr_Bsig[chr][(st, ed)][name]
                except:
                    score = 0
                s += "\t" + str(score)
            print >> f, s
            ID += 1
    f.close()

    print >> sys.stderr, "Done"

    #end = time.time()
    #print end - start

if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = ArgumentParser(description='Calculate bin signal')
    parser.add_argument(metavar = '--cov',
                        dest='cov_fname',
                        type=str,
                        help='signal cn file')
    parser.add_argument(metavar='-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence filename')
    parser.add_argument('--Bsize',
                        dest="bin_size",
                        type=int,
                        default=167,
                        help='Bin window in bp')
    parser.add_argument('--Bstep',
                        dest="bin_step",
                        type=int,
                        default=25,
                        help='Bin moving step in bp')
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
    parser.add_argument('--sample-choice',
                        dest="sample_choice",
                        type=str,
                        default='all',
                        help='sample data choice (default:all)')
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()

    # get length for each chromosome
    genome_size = {}

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
        chr_list = sorted(genome_size.keys())
    else:
        chr_list = sorted(args.chr_list)

    Bin_sig (args.cov_fname,
             genome_size,
             args.bin_size,
             args.bin_step,
             args.skip_zero,
             chr_list,
             args.sample_choice,
             args.out_fname)
