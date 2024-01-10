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

def NCP_count (peak_fname,
               cov_fname,
               peak_choice,
               genome_size,
               NCP_len,
               sigma,
               chr_list,
               skip_zero,
               out_fname):

    # read peak file
    print >> sys.stderr, "reading peak file"
    chr_peak = {}
    chr_peakpos = {}
    chr_Pcov = {}
    First = True
    for line in open(peak_fname):
        cols = line.strip().split()
        if First:
            label = cols[3:]
            First = False
            # make a list of target NCP peaks
            if peak_choice == 'input':
                target = [label[-1]] # target only control sample
            elif peak_choice == 'all':
                target = label # target all sample
            continue
        _, chr, pos = cols[0], cols[1], int(cols[2])
        if chr_list and chr not in chr_list:
            continue
        counts = cols[3:]
        for i in range(len(label)):
            name = label[i]
            if name not in target:
                continue
            count = int(round(float(counts[i])))
            if count <= 0:
                continue
            if chr not in chr_peak:
                chr_peak[chr] = {}
                chr_peakpos[chr] = []
                chr_Pcov[chr] = {}
            if pos not in chr_peak[chr]:
                chr_peak[chr][pos] = {}
                chr_peakpos[chr].append(pos)
                chr_Pcov[chr][pos] = {}
            assert name not in chr_peak[chr][pos]
            chr_peak[chr][pos][name] = count
            chr_Pcov[chr][pos][name] = 0
    for chr in chr_peakpos:
        chr_peakpos[chr] = sorted(chr_peakpos[chr])

    # reading coverage file and get NCP peak coverage if necessary
    if cov_fname:
        print >> sys.stderr, "reading coverage file"
        pt = 0
        prev_chr = None
        First = True
        order = 2
        for line in open(cov_fname):
            cols = line.strip().split()
            if First:
                label = cols[3:]
                First = False
                continue
            _, chr, pos = cols[0], cols[1], int(cols[2])
            if chr not in chr_peakpos:
                continue
            if prev_chr != chr:
                pt = 0
                prev_chr = chr
            if math.log10(pos+1) > order:
                print >> sys.stderr, chr + " pos " + str(pos) +" is reading"
                order += 1
            if pt >= len(chr_peakpos[chr]):
                continue
            if pos < chr_peakpos[chr][pt]:
                continue
            while pos > chr_peakpos[chr][pt] and pt < len(chr_peakpos[chr]):
                pt += 1
            if pos != chr_peakpos[chr][pt]:
                continue
            counts = cols[3:]
            for i in range(len(counts)):
                name = label[i]
                if name not in target:
                    continue
                count = int(round(float(counts[i])))
                chr_Pcov[chr][pos][name] += count
            pt += 1
        chr_peak = copy.deepcopy(chr_Pcov)
        

    # check start-end position of nucleosome with gaussian smoothing to get profile
    if sigma:
        chr_profile = {}
        for chr in sorted(chr_peakpos.keys()):
            for pos in chr_peakpos[chr]:
                for name in chr_peak[chr][pos]:
                    score = chr_peak[chr][pos][name]
                    for i in range(1, score+1):
                        offset = int(round(sigma* math.sqrt(2*math.log(float(score)/i))))
                        st = max(0, pos - offset)
                        ed = min(genome_size[chr], pos + offset + 1)
                        if chr not in chr_profile:
                            chr_profile[chr] = {}
                        if st not in chr_profile[chr]:
                            chr_profile[chr][st] = {}
                        if ed not in chr_profile[chr]:
                            chr_profile[chr][ed] = {}
                        if name not in chr_profile[chr][st]:
                            chr_profile[chr][st][name] = 0
                        if name not in chr_profile[chr][ed]:
                            chr_profile[chr][ed][name] = 0
                        chr_profile[chr][st][name] += 1
                        chr_profile[chr][ed][name] -= 1
    
    if NCP_len:
        chr_profile = {}
        for chr in sorted(chr_peakpos.keys()):
            for pos in chr_peakpos[chr]:
                st = pos - NCP_len/2
                ed = pos + NCP_len/2 + 1
                for name in chr_peak[chr][pos]:
                    score = chr_peak[chr][pos][name]
                    if chr not in chr_profile:
                        chr_profile[chr] = {}
                    if st not in chr_profile[chr]:
                        chr_profile[chr][st] = {}
                    if ed not in chr_profile[chr]:
                        chr_profile[chr][ed] = {}
                    if name not in chr_profile[chr][st]:
                        chr_profile[chr][st][name] = 0
                    if name not in chr_profile[chr][ed]:
                        chr_profile[chr][ed][name] = 0
                    chr_profile[chr][st][name] += score
                    chr_profile[chr][ed][name] -= score
            
    # summarize the output
    print >> sys.stderr, "writing occupancy file"
    
    f = open(out_fname + '_occ.cn', 'w')
    s = 'SNP\tChromosome\tPhysicalPosition'
    for i in range(len(target)):
        s += '\t' + target[i]
    print >> f, s

    ID = 0
    for chr in chr_profile:
        previous = [0]*len(target)
        if skip_zero:
            start = min(chr_profile[chr])
            end = max(chr_profile[chr]) + 1
        else:
            start = 0
            end = genome_size[chr]
        for i in xrange(start, end):
            s = str(ID) + "\t" + chr + "\t" + str(i)
            for j in range(len(target)):
                name = target[j]
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

    parser = ArgumentParser(description='Calculate nucleosome occupancy')
    parser.add_argument(metavar='--peak',
                        dest="peak_fname",
                        type=str,
                        help='peak cn file')
    parser.add_argument(metavar='-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence filename')
    parser.add_argument('--cov',
                        dest='cov_fname',
                        type=str,
                        help='coverage cn file')
    parser.add_argument('--peak-choice',
                        dest="peak_choice",
                        type=str,
                        default='all',
                        help='NCP peak data choice (default:all data, input:control only)')
    parser.add_argument('--Nlen',
                        dest="NCP_len",
                        type=int,
                        default=147,
                        help='Mono-nucleosomal length in bp')
    parser.add_argument('--sigma',
                        dest="sigma",
                        type=int,
                        nargs='?',
                        const=20,
                        help='Gaussain smoothing binwidth, default:20bp')
    parser.add_argument('--chr',
                        dest="chr_list",
                        type=str,
                        nargs='+',
                        help='tagert chromosome list')
    parser.add_argument('--skip',
                        dest="skip_zero",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='skip the zero occupancy positions')
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()

    # get length for each chromosome
    genome_size = {}

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

    if not args.cov_fname:
        cov_fname = None
    else:
        cov_fname = args.cov_fname

    if args.sigma:
        NCP_len = None
    else:
        NCP_len = args.NCP_len
        
    NCP_count (args.peak_fname,
               cov_fname,
               args.peak_choice,
               genome_size,
               NCP_len,
               args.sigma,
               chr_list,
               args.skip_zero,
               args.out_fname
               )
