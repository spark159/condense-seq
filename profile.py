import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import copy
import time
#import random

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

def NCP_count (sig_fname,
               table_fname,
               feature_choice,
               genome_size,
               chr_list,
               up_win,
               down_win,
               data_choice,
               skip_zero,
               out_fname):
    
    # read UCSC gene annotation table file
    print >> sys.stderr, "reading gene annotation table file"
    chr_start, chr_end = {}, {}
    chr_info = {}
    chr_win = {}
    feature_count = 0
    First = True
    for line in open(table_fname):
        cols = line.strip().split()
        if First:
            First = False
        _, ID, chr, strand = cols[:4]
        if chr not in chr_list:
            continue
        assert set(feature_choice) <= set(['TSS', 'TTS', 'CSS', 'CTS', 'ESS', 'ETS'])
        pos_list = []
        feat_list = []
        if 'TSS' in feature_choice:
            if strand == '+':
                pos_list.append(int(cols[4]))
            else:
                pos_list.append(int(cols[5]))
            feat_list.append('TSS')
        if 'TTS' in feature_choice:
            if strand == '+':
                pos_list.append(int(cols[5]))
            else:
                pos_list.append(int(cols[4]))
            feat_list.append('TTS')
        if 'CSS' in feature_choice:
            if strand == '+':
                pos_list.append(int(cols[6]))
            else:
                pos_list.append(int(cols[7]))
            feat_list.append('CSS')
        if 'CTS' in feature_choice:
            if strand == '+':
                pos_list.append(int(cols[7]))
            else:
                pos_list.append(int(cols[6]))
            feat_list.append('CTS')
        if 'ESS' in feature_choice:
            if strand == '+':
                for pos in cols[9].strip(',').split(','):
                    pos_list.append(int(pos))
                    feat_list.append('ESS')
            else:
                for pos in cols[10].strip(',').split(','):
                    pos_list.append(int(pos))
                    feat_list.append('ESS')            
        if 'ETS' in feature_choice:
            if strand == '+':
                for pos in cols[10].strip(',').split(','):
                    pos_list.append(int(pos))
                    feat_list.append('ETS')
            else:
                for pos in cols[9].strip(',').split(','):
                    pos_list.append(int(pos))
                    feat_list.append('ETS')
        assert len(pos_list) == len(feat_list)
        if chr not in chr_start:
            chr_start[chr] = []
            chr_end[chr] = []
            chr_info[chr] = []
            chr_win[chr] = []
        for k in range(len(pos_list)):
            #pos = pos_list[k] + random.randint(0,1000000)
            pos = pos_list[k]
            feature = feat_list[k]
            if strand == '+':
                start = max(0, pos - up_win)
                end = min(genome_size[chr], pos + down_win)
            else:
                start = max(0, pos - down_win)
                end = min(genome_size[chr], pos + up_win)
            if (start, end) not in chr_win[chr]:
                chr_win[chr].append((start, end))
                chr_start[chr].append(start)
                chr_end[chr].append(end)
                chr_info[chr].append([start,(ID, feature, chr, str(pos), strand)])
                feature_count += 1
            
    for chr in chr_start:
        chr_start[chr] = sorted(chr_start[chr])
        chr_end[chr] = sorted(chr_end[chr])
        chr_info[chr] = sorted(chr_info[chr])

    for chr in chr_info:
        temp = []
        for info in chr_info[chr]:
            temp.append(info[1])
        chr_info[chr] = temp

    # start writing profile file
    print >> sys.stderr, "writing profile file"
    f = open(out_fname + '_profile.txt', 'w')
    s = 'Sample\tID\tFeature\tChromosome\tPhysicalPosition\tStrand'
    for i in range(-up_win, down_win+1):
        s += '\t' + str(i)
    print >> f, s
    
    # read genomic signal cn file
    print >> sys.stderr, "reading signal cn file"
    stpt, edpt = 0, 0
    prev_chr = None
    domain = []
    domain_st = []
    infos = []
    First = True
    order = 2
    for line in open(sig_fname):
        cols = line.strip().split()
        if First:
            label = cols[3:]
            if data_choice == 'all':
                target = label
            elif data_choice == 'input':
                target = [label[-1]]
            First = False
            continue
        _, chr, pos = cols[0], cols[1], int(cols[2])
        if chr not in chr_start:
            continue
        if prev_chr != chr:
            while len(domain) > 0:
                temp = domain.pop()
                domain_st.pop()
                ID, feature_choice, chr, position, strand = infos.pop()
                if skip_zero and sum([sum(temp[u]) for u in range(len(temp))]) <= 0:
                    continue
                for i in range(len(target)):
                    name = target[i]
                    profile = temp[i]
                    if strand == '-':
                        profile = profile[::-1]
                    s = name + "\t" + ID + "\t" + feature_choice + "\t" + chr + "\t" + position + "\t" + strand + "\t"
                    for u in range(len(profile)-1):
                        value = profile[u]
                        s += str(value) + "\t"
                    s += str(profile[-1])
                    print >> f, s
            stpt = 0 
            edpt = 0
            prev_chr = chr
        if math.log10(pos+1) > order:
            print >> sys.stderr, chr + " pos " + str(pos) +" is reading"
            order += 1
        while stpt < len(chr_start[chr]) and pos >= chr_start[chr][stpt]:
            temp = []
            for i in range(len(target)):
                temp.append([0.0]*(up_win + down_win + 1))
            domain.insert(0, temp)
            domain_st.insert(0, chr_start[chr][stpt])
            info = chr_info[chr][stpt]
            infos.insert(0, info)
            stpt += 1
        while edpt < len(chr_end[chr]) and pos > chr_end[chr][edpt]:
            temp = domain.pop()
            domain_st.pop()
            ID, feature_choice, chr, position, strand = infos.pop()
            edpt += 1
            if skip_zero and sum([sum(temp[u]) for u in range(len(temp))]) <= 0:
                continue
            for i in range(len(target)):
                name = target[i]
                profile = temp[i]
                if strand == '-':
                    profile = profile[::-1]
                s = name + "\t" + ID + "\t" + feature_choice + "\t" + chr + "\t" + position + "\t" + strand + "\t"
                for u in range(len(profile)-1):
                    value = profile[u]
                    s += str(value) + "\t"
                s += str(profile[-1])
                print >> f, s
        if len(domain) <= 0:
            continue
        counts = cols[3:]
        for i in range(len(counts)):
            name = label[i]
            if name not in target:
                continue
            count = float(counts[i])
            for j in range(len(domain)):
                idx = pos - domain_st[j]
                domain[j][i][idx] = count
    while len(domain) > 0:
        temp = domain.pop()
        domain_st.pop()
        ID, feature_choice, chr, position, strand = infos.pop()
        if skip_zero and sum([sum(temp[u]) for u in range(len(temp))]) <= 0:
            continue
        for i in range(len(target)):
            name = target[i]
            profile = temp[i]
            if strand == '-':
                profile = profile[::-1]
            s = name + "\t" + ID + "\t" + feature_choice + "\t" + chr + "\t" + position + "\t" + strand + "\t"
            for u in range(len(profile)-1):
                value = profile[u]
                s += str(value) + "\t"
            s += str(profile[-1])
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

    parser = ArgumentParser(description='Make a profile matrix of input signals with respect to annotation')
    parser.add_argument(metavar='--sig',
                        dest="sig_fname",
                        type=str,
                        help='signal cn file')
    parser.add_argument(metavar='--table',
                        dest='table_fname',
                        type=str,
                        help='UCSC annotation table/Custom region file')
    parser.add_argument(metavar='-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence filename')
    parser.add_argument('--feature',
                        dest="feature_choice",
                        type=str,
                        nargs='+',
                        help='TSS: tx start site \nTTS: tx terminal site\nCSS: coding start site\nCTS: conding terminal site\nESS: exon start site\nETS: exon terminal site')
    parser.add_argument('--chr',
                        dest="chr_list",
                        type=str,
                        nargs='+',
                        help='tagert chromosome list')
    parser.add_argument('--up',
                        dest="up_win",
                        type=int,
                        default=1000,
                        help='down stream window size')
    parser.add_argument('--down',
                        dest="down_win",
                        type=int,
                        default=2000,
                        help='up stream window size')
    parser.add_argument('--data',
                        dest="data_choice",
                        type=str,
                        default='all',
                        help='sample data choice, input: control only, all:all data (default:all)')
    parser.add_argument('--skip',
                        dest="skip_zero",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='skip the zero profile')
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

    if not args.feature_choice:
        feature_choice = ['TSS']
    else:
        feature_choice = args.feature_choice

    NCP_count (args.sig_fname,
               args.table_fname,
               feature_choice,
               genome_size,
               chr_list,
               args.up_win,
               args.down_win,
               args.data_choice,
               args.skip_zero,
               args.out_fname
               )
