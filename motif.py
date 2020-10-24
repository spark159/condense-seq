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
               ref_fname,
               cov_fname,
               peak_choice,
               motif_len,
               chr_list,
               out_fname):

    # read peak file
    print >> sys.stderr, "reading peak file"
    chr_peak = {}
    chr_peakpos = {}
    chr_Pcov = {}
    chr_info = {}
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
        ID, chr, pos = cols[0], cols[1], int(cols[2])
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
                chr_info[chr] = []
            if pos not in chr_peak[chr]:
                chr_peak[chr][pos] = {}
                chr_peakpos[chr].append(pos)
                chr_Pcov[chr][pos] = {}
                chr_info[chr].append([pos, (ID, chr, pos)])
            assert name not in chr_peak[chr][pos]
            chr_peak[chr][pos][name] = count
            chr_Pcov[chr][pos][name] = 0
    
    chr_st = {}
    for chr in chr_peakpos:
        chr_peakpos[chr] = sorted(chr_peakpos[chr])
        chr_info[chr] = sorted(chr_info[chr])
        chr_info[chr] = [info for pos, info in chr_info[chr]]
        chr_st[chr] = [pos - motif_len/2 for pos in chr_peakpos[chr]]

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

    # start writing motif file
    print >> sys.stderr, "writing motif file"
    f = open(out_fname + '_motif.txt', 'w')
    s = 'Sample\tID\tChromosome\tPhysicalPosition\tStrand\tWeight\tSequence'
    print >> f, s
OA
    # read genome to get sequence
    print >> sys.stderr, "reading reference genome file"
    chrhit_count = 0
    Find = False
    for line in open(ref_fname):
        line = line.strip()
        if line.startswith(">"):
            while Find and len(left) > 0:
                idx, seq = left.pop(0)
                seq += 'N'*(motif_len-len(seq))
                ID, chr, peakpos = chr_info[chr][idx]
                for j in range(len(target)):
                    name = target[j]
                    try:
                        weight = chr_peak[chr][peakpos][name]
                    except:
                        continue
                    s = name + "\t" + ID + "\t" + chr + "\t" + str(peakpos) + "\t" + '+' + "\t" + str(weight)
                    s += '\t' + seq
                    print >> f, s
            Find = False
            if chrhit_count >= len(chr_st):
                break
            if line[1:] in chr_st:
                Find = True
                chrhit_count += 1
                k = 0
                chr = line[1:]
                pos = chr_st[chr][k]
                pt = -1
                seq = ""
                left = []
            continue
        if Find:
            if len(left) == 0 and pt + len(line) < pos:
                pt += len(line)
                continue
            for i in range(len(left)):
                idx, seq = left.pop(0)
                ed = min(len(line), motif_len-len(seq))
                seq += line[:ed]
                if len(seq) == motif_len:
                    ID, chr, peakpos = chr_info[chr][idx]
                    for j in range(len(target)):
                        name = target[j]
                        try:
                            weight = chr_peak[chr][peakpos][name]
                        except:
                            continue
                        s = name + "\t" + ID + "\t" + chr + "\t" + str(peakpos) + "\t" + '+' + "\t" + str(weight)
                        s += '\t' + seq
                        print >> f, s
                else:
                    left.append([idx, seq])
            while pt + len(line) >= pos and k < len(chr_st[chr]):
                if pos < 0:
                    seq = 'N' * abs(pos)
                    seq += line[0:min(motif_len-abs(pos),len(line))]
                else:
                    loc = pos - pt - 1
                    seq = line[loc:min(loc+motif_len,len(line))]
                if len(seq) == motif_len:
                    ID, chr, peakpos = chr_info[chr][k]
                    for j in range(len(target)):
                        name = target[j]
                        try:
                            weight = chr_peak[chr][peakpos][name]
                        except:
                            continue
                        s = name + "\t" + ID + "\t" + chr + "\t" + str(peakpos) + "\t" + '+' + "\t" + str(weight)
                        s += '\t' + seq
                        print >> f, s
                else:
                    left.append([k, seq])
                k += 1
                try:
                    pos = chr_st[chr][k]
                except:
                    None
            if chrhit_count >= len(chr_st) and k >= len(chr_st[chr]) and len(left) == 0:
                break
            pt += len(line)
    while Find and len(left) > 0:
        idx, seq = left.pop(0)
        seq += 'N'*(motif_len-len(seq))
        ID, chr, peakpos = chr_info[chr][idx]
        for j in range(len(target)):
            name = target[j]
            try:
                weight = chr_peak[chr][peakpos][name]
            except:
                continue
            s = name + "\t" + ID + "\t" + chr + "\t" + str(peakpos) + "\t" + '+' + "\t" + str(weight)
            s += '\t' + seq
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

    parser = ArgumentParser(description='Get motif sequences')
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
    parser.add_argument('--Mlen',
                        dest="motif_len",
                        type=int,
                        default=147,
                        help='motif length in bp (should be odd number)')
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

    if args.motif_len % 2 == 0:
        print >> sys.stderr, "Error: motif length should be odd number."
        sys.exit(1)
    
    chr_list = []
    if not args.chr_list:
        chr_list = None
    else:
        chr_list = sorted(args.chr_list)

    if not args.cov_fname:
        cov_fname = None
    else:
        cov_fname = args.cov_fname

    NCP_count (args.peak_fname,
               args.ref_fname,
               cov_fname,
               args.peak_choice,
               args.motif_len,
               chr_list,
               args.out_fname
               )
