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

def NCP_occ (fnames,
             genome_size,
             min_len,
             max_len,
             mm_cutoff,
             NCP_len,
             sigma,
             skip_zero,
             chr_list,
             out_fname):

    # gather whole file names
    data, label = [], []
    for fname in fnames:
        data.append(fname)
        label.append(fname.rsplit('/', 1)[-1].split('.')[0])

    # make genome start-end profile dictionary
    chr_mid = {}
    
    # open the sam file
    for i in range(len(data)):
        name = label[i]
        #chr_cov = out_cov[i]
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

            # invalid: not right size of read 
            if abs(tlen) < min_len or abs(tlen) > max_len:
                #type = 'invalid:mutant'
                continue

            # record mid-point of reads
            if tlen > 0: # left read
                end_pos = pos + tlen
                if tlen % 2 != 0:
                    midscore_list = [[pos+tlen/2, 1]]
                else:
                    midscore_list = [[pos+tlen/2-1, 0.5], [pos+tlen/2, 0.5]]
            else: # right read
                end_pos = pos
                cigar_str=re.split('(\d+)',cigar_str)[1:]
                for i in range(len(cigar_str)/2):
                    s = cigar_str[2*i+1]
                    num = int(cigar_str[2*i])
                    if s == 'M' or s == 'D':
                        end_pos += num
                pos = end_pos + tlen
                if tlen % 2 != 0:
                    midscore_list = [[end_pos+tlen/2-1, 1]]
                else:
                    midscore_list = [[end_pos+tlen/2-1, 0.5], [end_pos+tlen/2, 0.5]]

            if ref_id not in chr_mid:
                chr_mid[ref_id] = {}
            for mid, score in midscore_list:
                if mid not in chr_mid[ref_id]:
                    chr_mid[ref_id][mid] = {}
                if name not in chr_mid[ref_id][mid]:
                    chr_mid[ref_id][mid][name] = 0
                chr_mid[ref_id][mid][name] += score

    chr_profile = {}
    # mark start-end position of nucleosome
    # no smoothing
    if NCP_len:
        for chr in sorted(chr_mid.keys()):
            for pos in chr_mid[chr]:
                st = pos - NCP_len/2
                ed = pos + NCP_len/2 + 1
                for name in chr_mid[chr][pos]:
                    score = chr_mid[chr][pos][name]
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
                    
    # gaussain smoothing
    if sigma:
        for chr in sorted(chr_mid.keys()):
            for pos in chr_mid[chr]:
                for name in chr_mid[chr][pos]:
                    score = chr_mid[chr][pos][name]
                    for i in range(int(math.floor(score))):
                        offset = int(round(sigma* math.sqrt(2*math.log(float(score)/(i+1)))))
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
                
    # summarize the output
    print >> sys.stderr, "writing occupancy file"

    f = gzip.open(out_fname + '_occ.gtab.gz', 'wb')
    s = 'Chromosome\tPosition'
    for i in range(len(label)):
        s += '\t' + label[i]
    print >> f, s

    #ID = 0
    for chr in sorted(chr_profile.keys(), cmp=chr_cmp):
        previous = [0]*len(label)
        if skip_zero:
            start = min(chr_profile[chr])
            end = max(chr_profile[chr]) + 1
        else:
            start = 0
            end = genome_size[chr]
        for i in xrange(start, end):
            #s = str(ID) + "\t" + chr + "\t" + str(i)
            s = chr + "\t" + str(i)
            for j in range(len(label)):
                name = label[j]
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

    parser = ArgumentParser(description='Calculate occupancy along the genome')
    parser.add_argument(metavar='-f1',
                        dest="fnames",
                        type=str,
                        nargs='+',
                        help='SAM/Bam filenames')
    parser.add_argument(metavar='-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence filename')
    parser.add_argument('--min',
                        dest="min_len",
                        type=int,
                        default=126,
                        help='minimum length for selection in bp, default:126bp')
    parser.add_argument('--max',
                        dest="max_len",
                        type=int,
                        default=184,
                        help='maximum length for selection in bp, default:184bp')
    parser.add_argument('-m',
                        dest="mm_cutoff",
                        type=int,
                        default=10,
                        help='mismatch cut-off in bp')
    parser.add_argument('--Nlen',
                        dest="NCP_len",
                        type=int,
                        default=73,
                        help='Mono-nucleosomal occupancy length in bp, default:73bp')
    parser.add_argument('--sigma',
                        dest="sigma",
                        type=int,
                        nargs='?',
                        const=20,
                        default=None,
                        help='Gaussain smoothing binwidth, default:20bp')
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
            key = line.split()[0][1:]
            assert key not in genome_size
            genome_size[key] = 0
            xcontinue
        genome_size[key] += len(line)

    chr_list = []
    if not args.chr_list:
        chr_list = sorted(genome_size.keys(), cmp=chr_cmp)
    else:
        chr_list = sorted(args.chr_list, cmp=chr_cmp)

    if args.sigma:
        NCP_len = None
    else:
        NCP_len = args.NCP_len

    NCP_occ (args.fnames,
             genome_size,
             args.min_len,
             args.max_len,
             args.mm_cutoff,
             NCP_len,
             args.sigma, 
             args.skip_zero,
             chr_list,
             args.out_fname
             )
