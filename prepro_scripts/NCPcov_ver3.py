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

def NCP_count (peak_fname,
               cov_fname,
               peak_choice,
               NCP_len,
               chr_list,
               out_fname):

    # read peak file
    print >> sys.stderr, "reading peak file"
    chr_peak = {}
    First = True
    
    for line in gzopen(peak_fname):
        cols = line.strip().split()
        if First:
            label = cols[2:]
            First = False
            continue
        chr, pos = cols[0], int(cols[1])
        if chr_list and chr not in chr_list:
            continue
        if chr not in chr_peak:
            chr_peak[chr] = {}
        if pos not in chr_peak[chr]:
            chr_peak[chr][pos] = {}
        counts = cols[2:]
        for i in range(len(label)):
            name = label[i]
            count = float(counts[i])
            if count <=0:
                continue
            if name not in chr_peak[chr][pos]:
                chr_peak[chr][pos][name] = 0
            chr_peak[chr][pos][name] += count

    # make a list of target NCP peaks
    if peak_choice == 'input':
        target = [label[-1]] # target only control peaks
    elif peak_choice == 'all':
        target = label # target all peaks
        
    chr_st, chr_ed = {}, {}
    for chr in sorted(chr_peak.keys(), cmp=chr_cmp):
        for NCPpos in sorted(chr_peak[chr].keys()):
            for name in target:
                try:
                    score = chr_peak[chr][NCPpos][name]
                    if chr not in chr_st:
                        chr_st[chr] = []
                        chr_ed[chr] = []
                    st = NCPpos - NCP_len/2
                    ed = NCPpos + NCP_len/2
                    if st not in chr_st:
                        chr_st[chr].append(st)
                        chr_ed[chr].append(ed)
                except:
                    continue           

    # reading coverage file and get NCP coverage
    print >> sys.stderr, "reading coverage file"
    chr_Ncov = {}
    NCPcovs = []
    
    First = True
    order = None
    prev_chr = None

    for line in gzopen(cov_fname):
        cols = line.strip().split()
        if First:
            label = cols[2:]
            First = False
            continue
        chr, pos = cols[0], int(cols[1])
        if chr not in chr_st:
            continue
        if prev_chr != None and chr != prev_chr:
            while chr_st[prev_chr]:
                chr_st[prev_chr].pop(0)
                NCPcovs.insert(0, [0]*len(label))
            while chr_ed[prev_chr]:
                ed = chr_ed[prev_chr].pop(0)
                dyad = ed - NCP_len/2
                covs = NCPcovs.pop()
                if prev_chr not in chr_Ncov:
                    chr_Ncov[prev_chr] = {}
                chr_Ncov[prev_chr][dyad] = {}
                for k in range(len(label)):
                    name = label[k]
                    chr_Ncov[prev_chr][dyad][name] = covs[k]
            prev_chr = chr
        if prev_chr == None:
            prev_chr = chr
        if order == None:
            print >> sys.stderr, chr + " pos " + str(pos) +" is reading"
            order = int(math.log10(max(pos, 1)))
        elif int(math.log10(max(pos, 1))) > order:
            print >> sys.stderr, chr + " pos " + str(pos) +" is reading"
            order += 1
        while chr_st[chr] and pos >= chr_st[chr][0]:
            chr_st[chr].pop(0)
            NCPcovs.insert(0, [0]*len(label))
        while chr_ed[chr] and pos > chr_ed[chr][0]:
            ed = chr_ed[chr].pop(0)
            dyad = ed - NCP_len/2
            covs = NCPcovs.pop()
            if chr not in chr_Ncov:
                chr_Ncov[chr] = {}
            chr_Ncov[chr][dyad] = {}
            for k in range(len(label)):
                name = label[k]
                chr_Ncov[chr][dyad][name] = covs[k]
        if len(NCPcovs) <= 0:
            continue
        counts = cols[2:]
        for k in range(len(NCPcovs)):
            for u in range(len(label)):
                count = int(counts[u])
                NCPcovs[k][u] += count
    while chr_st[chr]:
        chr_st[chr].pop(0)
        NCPcovs.insert(0, [0]*len(label))
    while chr_ed[chr]:
        ed = chr_ed[chr].pop(0)
        dyad = ed - NCP_len/2
        covs = NCPcovs.pop()
        if chr not in chr_Ncov:
            chr_Ncov[chr] = {}
        chr_Ncov[chr][dyad] = {}
        for k in range(len(label)):
            name = label[k]
            chr_Ncov[chr][dyad][name] = covs[k]
    
    # write NCP coverage file
    print >> sys.stderr, "writing NCP coverage file"
    #print "writing NCP file"
    
    f = gzip.open(out_fname + '_Ncov.gtab.gz', 'wb')
    s = 'Chromosome\tPosition'
    for i in range(len(label)):
        s += '\t' + label[i]
    print >> f, s

    #ID = 0
    for chr in sorted(chr_Ncov.keys(), cmp=chr_cmp):
        for NCPpos in sorted(chr_Ncov[chr].keys()):
            #s = str(ID) + "\t" + chr + "\t" + str(NCPpos)
            s = chr + "\t" + str(NCPpos)
            for name in label:
                try:
                    score = chr_Ncov[chr][NCPpos][name]
                except:
                    score = 0
                s += "\t" + str(score)
            print >> f, s
            #ID += 1
        
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

    parser = ArgumentParser(description='Calculate NCP coverage')
    parser.add_argument(metavar='--peak',
                        dest="peak_fname",
                        type=str,
                        help='peak gtab file')
    parser.add_argument(metavar = '--cov',
                        dest='cov_fname',
                        type=str,
                        help='coverage gtab file')
    parser.add_argument('--peak-choice',
                        dest="peak_choice",
                        type=str,
                        default='input',
                        help='NCP peak data choice (default:input control only, "all":all data)')
    parser.add_argument('--Nlen',
                        dest="NCP_len",
                        type=int,
                        default=171,
                        help='Mono-nucleosomal window in bp')
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
        chr_list = sorted(args.chr_list, cmp=chr_cmp)

    NCP_count (args.peak_fname,
               args.cov_fname,
               args.peak_choice,
               args.NCP_len,
               chr_list,
               args.out_fname
               )
