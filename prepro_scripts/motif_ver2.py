import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import copy
import gzip

def chr_cmp (chr_name1, chr_name2):
    assert chr_name1.startswith('chr')
    assert chr_name2.startswith('chr')
    chr_num1 = chr_name1[3:]
    try:
        chr_num1 = int(chr_num1)
    except:
        pass
    chr_num2 = chr_name2[3:]
    try:
        chr_num2 = int(chr_num2)
    except:
        pass
    if chr_num1 < chr_num2:
        return -1
    elif chr_num1 > chr_num2:
        return 1
    return 0

def gzopen (fname):
    if fname.endswith('.gz'):
        reading_file = gzip.open(fname, 'rb')
    else:
        reading_file = open(fname, 'r')
    return reading_file

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

def NCP_motif (peak_fname,
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
    for line in gzopen(peak_fname):
        cols = line.strip().split()
        if First:
            label = cols[2:]
            First = False
            # make a list of target NCP peaks
            if peak_choice == 'input':
                target = [label[-1]] # target only control sample
            elif peak_choice == 'all':
                target = label # target all sample
            continue
        chr, pos = cols[0], int(cols[1])
        if chr_list and chr not in chr_list:
            continue
        counts = cols[2:]
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
        First = True
        prev_chr = None
        order = None
        for line in gzopen(cov_fname):
            cols = line.strip().split()
            if First:
                label = cols[2:]
                First = False
                continue
            chr, pos = cols[0], int(cols[1])
            if chr not in chr_peakpos:
                continue
            if prev_chr != chr:
                pt = 0
                prev_chr = chr
            if order == None:
                print >> sys.stderr, chr + " pos " + str(pos) +" is reading"
                order = int(math.log10(max(pos, 1)))
            elif int(math.log10(max(pos, 1))) > order:
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
            counts = cols[2:]
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
    f = gzip.open(out_fname + '_motif.txt.gz', 'wb')
    s = 'Sample\tChromosome\tPosition\tStrand\tWeight\tSequence'
    print >> f, s

    # read genome to get sequence
    print >> sys.stderr, "reading reference genome file"
    chrhit_count = 0
    Find = False
    for line in gzopen(ref_fname):
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
                    s = name + "\t" + chr + "\t" + str(peakpos) + "\t" + '+' + "\t" + str(weight)
                    s += '\t' + seq
                    print >> f, s
            Find = False
            if chrhit_count >= len(chr_st):
                break
            if line.split()[0][1:] in chr_st:
                Find = True
                chrhit_count += 1
                k = 0
                chr = line.split()[0][1:]
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
                        s = name + "\t" + chr + "\t" + str(peakpos) + "\t" + '+' + "\t" + str(weight)
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
                        s = name + "\t" + chr + "\t" + str(peakpos) + "\t" + '+' + "\t" + str(weight)
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
            s = name + "\t" + chr + "\t" + str(peakpos) + "\t" + '+' + "\t" + str(weight)
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
                        help='peak gtab file')
    parser.add_argument(metavar='-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence filename')
    parser.add_argument('--cov',
                        dest='cov_fname',
                        type=str,
                        help='coverage gtab file')
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
        chr_list = sorted(args.chr_list, cmp=chr_cmp)

    if not args.cov_fname:
        cov_fname = None
    else:
        cov_fname = args.cov_fname

    NCP_motif (args.peak_fname,
               args.ref_fname,
               cov_fname,
               args.peak_choice,
               args.motif_len,
               chr_list,
               args.out_fname
               )
