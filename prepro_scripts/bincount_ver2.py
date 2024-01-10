import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import copy

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

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

def bin_count (fnames,
               genome_size,
               win_size,
               skip_zero,
               bin_GC,
               tlen_option,
               mm_cutoff,
               min_len,
               max_len,
               chr_list,
               out_fname):

    # gather whole file names
    data, label = [], []
    for fname in fnames:
        data.append(fname)
        label.append(fname.rsplit('/', 1)[-1].split('.')[0])

    # partition the genome
    chr_bins = {}
    for chr in chr_list:
        chr_bins[chr] = [0] * int(math.ceil(float(genome_size[chr]) / win_size))
        
    # output genome count dictionary
    out = [copy.deepcopy(chr_bins) for i in range(len(data))]

    # mean tlen per each bin for control reads
    if tlen_option:
        g_tlen = copy.deepcopy(chr_bins)
    
    # open the sam file
    for k in range(len(data)):
        g_count = out[k]      
        filename = data[k]
        print >> sys.stderr, "reading %s" % (filename)
        samtools_cmd = ["samtools",  "view", "-F 0x10", filename]
        samtools_proc = subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE, stderr=open("/dev/null", 'w'))
        for line in samtools_proc.stdout:
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
            if ref_id not in chr_list:
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

            # invalid: mate pairing failure
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

            # invalid: not right size of read 
            if abs(tlen) < min_len or abs(tlen) > max_len:
                #type = 'invalid:mutant'
                continue

            # get NCP position as the center of full read
            if tlen > 0: # left read
                # split the center if tlen is even
                if tlen % 2 != 0:
                    NCPscore = [[pos+tlen/2, 1]]
                else:
                    #NCPscore = [[pos+tlen/2-1, 0.5], [pos+tlen/2, 0.5]]
                    NCPscore = [[pos+tlen/2-1, 1]] # don't worry about even tlen case
            else: # right read
                end_pos = pos
                cigar_str=re.split('(\d+)',cigar_str)[1:]
                for i in range(len(cigar_str)/2):
                    s = cigar_str[2*i+1]
                    num = int(cigar_str[2*i])
                    if s == 'M' or s == 'D':
                        end_pos += num

                # split the center if tlen is even
                if tlen % 2 != 0:
                    NCPscore = [[end_pos+tlen/2-1, 1]]
                else:
                    #NCPscore = [[end_pos+tlen/2-1, 0.5], [pos+tlen/2, 0.5]]
                    NCPscore = [[end_pos+tlen/2-1, 1]] # don't worry about even tlen case

            # collect valid data
            for NCPpos, score in NCPscore:
                index = int(NCPpos) / int(win_size)
                g_count[ref_id][index] += score

                # record tlen for control
                if tlen_option and k == len(data)-1:
                    g_tlen[ref_id][index] += abs(tlen)

    # summarize the output
    print >> sys.stderr, "writing bin file"

    f = open(out_fname + '_bin.cn', 'w')
    #s = 'SNP\tChromosome\tPhysicalPosition'
    s = 'BinID\tChromosome\tStart\tEnd'
    for i in range(len(label)):
        s += '\t' + label[i]
    if bin_GC != None:
        s += '\t' + 'GCcontent'
    if tlen_option:
        s += '\t' + 'Meantlen'
    print >> f, s
        

    ID = 0
    for chr in sorted(g_count.keys(), cmp=chr_cmp):
        for i in range(len(g_count[chr])):
            #Binpos = i + win_size/2
            #s = str(ID) + "\t" + chr + "\t" + str(Binpos)
            st, ed = i*win_size, (i+1)*win_size 
            s = str(ID) + "\t" + chr + "\t" + str(st) + '\t' + str(ed)
            total = 0
            for j in range(len(out)):
                total += out[j][chr][i]
                s += '\t%s' % (out[j][chr][i])
            if skip_zero and total <= 0:
                continue
            if bin_GC != None:
                BinID = chr + '_' + str(i)
                s += '\t%f' % (bin_GC[BinID])
            if tlen_option:
                if out[-1][chr][i] > 0:
                    mean_tlen = float(g_tlen[chr][i])/out[-1][chr][i]
                else:
                    mean_tlen = 0
                s += '\t%f' % (mean_tlen)
            ID += 1
            print >> f, s
    f.close()
    

if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = ArgumentParser(description='Divide the genome into bins and counts the read number')
    parser.add_argument(metavar='-f',
                        dest="fnames",
                        type=str,
                        nargs='+',
                        help='SAM/Bam filenames')
    parser.add_argument(metavar='-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence filename')
    parser.add_argument('-w',
                        dest="win_size",
                        type=int,
                        default=1001,
                        help='bin window size in bp')
    parser.add_argument('-m',
                        dest="mm_cutoff",
                        type=int,
                        default=10,
                        help='mismatch cut-off in bp')
    parser.add_argument('--min',
                        dest="min_len",
                        type=int,
                        default=0,
                        help='minimum length for selection in bp')
    parser.add_argument('--max',
                        dest="max_len",
                        type=int,
                        default=sys.maxint,
                        help='maximum length for selection in bp')
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
                        help='skip the zero count bins')
    parser.add_argument('--gc',
                        dest="GC_option",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='GC content option')
    parser.add_argument('--tlen',
                        dest="tlen_option",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='tlen option')    
    parser.add_argument('-o',
                        dest='out_fname',
                        default='bc_out',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()

    # get length for each chromosome
    genome_size = {}

    for line in open(args.ref_fname):
        line = line.strip()
        if line.startswith('>'):
            key = line.split()[0][1:]
            assert key not in genome_size
            genome_size[key] = 0
            continue
        genome_size[key] += len(line)

    # get GC content of the binned genome
    def GC_content(seq):
        seq = seq.upper()
        output = 0.0
        for nt in seq:
            if nt in "GC":
                output += 1.0
        return output/len(seq)
    def add (chr, num, seq, dic):
        BinID = chr + "_" + str(num)
        GC = GC_content(seq)
        assert BinID not in dic
        dic[BinID] = GC

    win_size = args.win_size
    bin_GC = None
    if args.GC_option:
        bin_GC = {}
        for line in open(args.ref_fname):
            line = line.strip()
            if line.startswith(">"):
                if len(bin_GC) > 0 and len(seq) > 0:
                    add(chr, num, seq, bin_GC)
                chr, num = line.split()[0][1:], 0
                seq = ""
                remain = ""
                continue
            if len(seq) < win_size:
                seq += line
            if len(seq) > win_size:
                seq, remain = seq[:win_size], seq[win_size:]
            if len(seq) == win_size:
                add(chr, num, seq, bin_GC)
                seq = remain[:]
                remain = ""
                num += 1
        if len(seq) > 0:
            add(chr, num, seq, bin_GC)    
    
    chr_list = []
    if not args.chr_list:
        chr_list = sorted(genome_size.keys(), cmp=chr_cmp)
    else:
        chr_list = sorted(args.chr_list, cmp=chr_cmp)

    #print chr_list

    bin_count (args.fnames,
               genome_size,
               win_size,
               args.skip_zero,
               bin_GC,
               args.tlen_option,
               args.mm_cutoff,
               args.min_len,
               args.max_len,
               chr_list,
               args.out_fname
               )
