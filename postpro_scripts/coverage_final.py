import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import copy

def get_coverage (filenames,
                  chr_length,
                  min_len,
                  max_len,
                  mm_cutoff,
                  skip_zero,
                  chr_list,
                  out_fname):

    # collect filenames
    labels = [filename.rsplit('.',1)[0] for filename in filenames]

    # make genome start-end profile dictionary
    chr_profile = {}
    
    # open the SAM/BAM files
    for filename, label in zip(filenames, labels):
        print >> sys.stderr, "reading %s" % (filename)
        samtools_cmd = ["samtools",  "view", "-F 0x10", filename]
        samtools_proc = subprocess.Popen(samtools_cmd,
                                         stdout=subprocess.PIPE,
                                         stderr=open("/dev/null", 'w'))

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

            # get NCP position and coverage range
            if tlen > 0: # left-end read
                end_pos = pos + tlen
            else: # skip right-end read to count the read-pair only once
                continue
                #end_pos = copy.deepcopy(pos)
                #cigar_str=re.split('(\d+)',cigar_str)[1:]
                #for i in range(len(cigar_str)/2):
                #    s = cigar_str[2*i+1]
                #    num = int(cigar_str[2*i])
                #    if s == 'M' or s == 'D':
                #        end_pos += num
                #assert abs(tlen) == end_pos - pos

            # record start-end positions of reads
            if ref_id not in chr_profile:
                chr_profile[ref_id] = {}
            if pos not in chr_profile[ref_id]:
                chr_profile[ref_id][pos] = {}
            if label not in chr_profile[ref_id][pos]:
                chr_profile[ref_id][pos][label] = 0
            if end_pos not in chr_profile[ref_id]:
                chr_profile[ref_id][end_pos] = {}
            if label not in chr_profile[ref_id][end_pos]:
                chr_profile[ref_id][end_pos][label] = 0
            chr_profile[ref_id][pos][label] += 1
            chr_profile[ref_id][end_pos][label] -= 1
            
    # summarize the output
    print >> sys.stderr, "writing coverage file"
    
    f = open(out_fname + '_cov.pdata', 'w')
    s = 'Chromosome\tPosition'
    for i in range(len(labels)):
        s += '\t' + labels[i]
    print >> f, s

    for chr in chr_list:
        previous = [0]*len(labels)
        if skip_zero:
            start = min(chr_profile[chr])
            end = max(chr_profile[chr]) + 1
        else:
            start = 0
            end = chr_length[chr]
        for i in xrange(start, end):
            s = chr + "\t" + str(i)
            for j in range(len(labels)):
                label = labels[j]
                past = previous[j]
                try:
                    current = past + chr_profile[chr][i][label]
                except:
                    current = past + 0
                s += "\t" + str(current)
                previous[j] = current
            if skip_zero and sum(previous) == 0:
                continue
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

    parser = ArgumentParser(description='Calculate coverage along the genome')
    parser.add_argument(metavar='-f',
                        dest="filenames",
                        type=str,
                        nargs='+',
                        help='SAM/Bam filenames')
    parser.add_argument(metavar='-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence filename')
    parser.add_argument('-m',
                        dest="mm_cutoff",
                        type=int,
                        default=10,
                        help='mismatch cut-off in bp')
    parser.add_argument('--min',
                        dest="min_len",
                        type=int,
                        nargs='?',
                        default=0,
                        const=120,
                        help='minimum length for selection in bp')
    parser.add_argument('--max',
                        dest="max_len",
                        type=int,
                        nargs='?',
                        default=sys.maxint,
                        const=170,
                        help='maximum length for selection in bp')
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
    chr_length = {}

    if not args.ref_fname:
        print >> sys.stderr, "Error: there is no reference file input."
        sys.exit(1)

    for line in open(args.ref_fname):
        line = line.strip()
        if line.startswith('>'):
            chr = line[1:]
            assert chr not in chr_length
            chr_length[chr] = 0
            continue
        chr_length[chr] += len(line)

    chr_list = []
    if not args.chr_list:
        chr_list = chr_length.keys()
    else:
        chr_list = sorted(args.chr_list)

    #print chr_list

    get_coverage (args.filenames,
                  chr_length,
                  args.min_len,
                  args.max_len,
                  args.mm_cutoff,
                  args.skip_zero,
                  chr_list,
                  args.out_fname
                  )
