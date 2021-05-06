import os, sys, subprocess, re
from argparse import ArgumentParser, FileType

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

def sort_seqs (read_fname,
               ref_fname,
               out_fname,
               mm_cutoff,
               bab_seq,
               win_sted):

    # get parameters
    all_wins = sorted(win_sted.keys()) # all sorted window names

    sted_wins = {} # window positions
    for win, pos_list in win_sted.items():
        for pos in pos_list:
            if pos not in sted_wins:
                sted_wins[pos] = []
            sted_wins[pos].append(win)

    N_len = 0 # total ambiguous nt counts
    for nt in bab_seq:
        if nt not in 'ATCG':
            N_len +=1

    
    # sequence alignment by Bowtie2
    aligner_cmd = ["bowtie2", '-x', 'temporal', '-U', read_fname ]
    aligner_cmd += ['--n-ceil', 'L,'+str(N_len)+','+str(0.15)] #turn off ambiguous nt limit
    aligner_cmd += ['--score-min', 'L,' + str(-100) + ',' +str(-100)] #turn off min score limit
    aligner_cmd += ['-N', str(1), '-L', str(1), '-i',  'S,1,0'] #turn off heuristic seeding
    aligner_cmd += ['--np', '0'] #turn off penalty for ambiguous nt
    align_proc = subprocess.Popen(aligner_cmd, stdout=subprocess.PIPE, stderr=open("/dev/null", 'w'))
    
    # start sort the reads
    seq_sort_fname = open(out_fname+".sort",'w')
    
    read_count = 0
    for line in align_proc.stdout:
        if line.startswith('@'):
            continue
        
        #print line
        defects = []

        cols = line.strip().split()
        read_id, flag, ref_id, pos, mapQ, cigar_str = cols[:6]
        read_id=":".join(read_id.split(':')[3:7])
        read_count +=1
        flag, pos = int(flag), int(pos)
        pos-=1

        # invalid: mapping failure
        if pos < 0 or flag & 0x4 != 0:
            defects.append('Mappingfail')
            #continue
        
        read_seq, qual =cols[9:11]
        ref_id = ref_id.strip()    

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
        if NM > mm_cutoff + N_len:
            defects.append('largeEdit')
            #continue

        # collect invalid data 
        if len(defects) > 0:
            type = 'invalid:' + '/'.join(defects)
            print >> seq_sort_fname, "@%s::%s" % (read_id, type)
            print >> seq_sort_fname, read_seq
            continue        
                    
        # get aligned read position at the start and end of the target windows
        # include a hit on the first nt of reads (start, end)
        # include a hit on the last nt of reads (only start)
        ref_pt, read_pt = pos-1, -1
        win_readpos = {}
        cigar_str=re.split('(\d+)',cigar_str)[1:]

        for i in range(len(cigar_str)/2):
            s = cigar_str[2*i+1]
            num = int(cigar_str[2*i])
            for j in range(num):
                if s=='M':
                    ref_pt += 1
                    read_pt += 1
                elif s=='I':
                    read_pt += 1
                elif s=='D':
                    ref_pt += 1
                    
                # record hit positions
                if ref_pt in sted_wins:
                    for win in sted_wins[ref_pt]:
                        if win not in win_readpos:
                            win_readpos[win] = []
                        if s == 'D':
                            win_readpos[win].append(read_pt+1)
                        else:
                            win_readpos[win].append(read_pt)

        # get bracodes or indexes of reads
        win_seq = {}
        for win in all_wins:
            try:
                st, ed = sorted(win_readpos[win])
                win_seq[win] = read_seq[st:ed]
            except:
                defects.append("Missing-WIN%d" % (win))

        # collect invalid data (missing windows)
        if len(defects) > 0:
            type = 'invalid:' + '/'.join(defects)
            print >> seq_sort_fname, "@%s::%s" % (read_id, type)
            print >> seq_sort_fname, read_seq
            continue        

        # otherwise, valid data
        type = 'valid:'
        for win in all_wins:
            type += '[' + win_seq[win] + ']'
        print >> seq_sort_fname, "@%s::%s" % (read_id, type)
        print >> seq_sort_fname, read_seq
             
    seq_sort_fname.close()
    
    # remove temporal files
    subprocess.call(['rm', 'temporal.ref'])
    subprocess.call('rm ' + 'temporal' + '*bt2', shell=True)


if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')
                
    parser = ArgumentParser(description='barcoding/indexing the Oncohistone library data')
    parser.add_argument(dest="filename",
                        type=str,
                        help='read file name (single read fastq file)')
    parser.add_argument(dest='ref_fname',
                        type=str,
                        help='reference prefix filename')
    parser.add_argument('-o',
                        dest='out_fname',
                        type=str,
                        help='output prefix filename')
    parser.add_argument('-m',
                        dest="mm_cutoff",
                        type=int,
                        default=10,
                        help='minimun cutoff mismatch or indel length')

    args = parser.parse_args()

    read_fname = args.filename
    
    # read reference file
    if not args.ref_fname:
        print >> sys.stderr, "Error: there is no reference file input."
        sys.exit(1)

    ref_seq = {}
    bab_seq = None
    for line in open(args.ref_fname + '.ref'):
        line = line.strip()
        if line.startswith('>'):
            ref_id = line[1:].strip()
            continue
        if ref_id == 'BACKBONE':
            bab_seq = line.strip().upper()
            continue
        assert ref_id not in ref_seq
        ref_seq[ref_id] = line.strip().upper()

    if not bab_seq:
        print >> sys.stderr, "Error: there is no backbone sequence in reference file."
        sys.exit(1)

        
    # find the locations of target windows
    win_sted = {}
    i = 0
    w_count = 0
    while i < len(bab_seq):
        if bab_seq[i] not in 'ATCG':
            nt = bab_seq[i]
            j = i + 1
            while j < len(bab_seq):
                if bab_seq[j] != nt:
                    break
                j +=1
            win_sted[w_count] = [i, j]
            w_count +=1
            i = j
        else:
            i +=1

    if not win_sted:
        print >> sys.stderr, "Error: there is no target windows in backbone sequence."
        sys.exit(1)

        
    # make a temporal build files for bowtie2 alignment
    f = open("temporal.ref", "w")
    print >> f, ">BACKBONE"
    print >> f, bab_seq
    f.close()
    subprocess.call(["bowtie2-build", "temporal.ref", "temporal"], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

    
    # set output file name
    if not args.out_fname:
        out_fname = read_fname.rsplit('.', 1)[0]
    else:
        out_fname = args.out_fname

        
    sort_seqs (read_fname,
               args.ref_fname,
               out_fname,
               args.mm_cutoff,
               bab_seq,
               win_sted)
