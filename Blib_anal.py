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
               ref_length):

    #print read_fname
    # seq alignment
    aligner_cmd=["bowtie2", '-x', ref_fname, '-U', read_fname ]
    #aligner_cmd=["bowtie2", '-x', ref_fname, '-1', read_fname+'1_001.fastq', '-2', read_fname+'2_001.fastq' ]

    #aligner_cmd += ['--n-ceil', 'L,'+str(insert_len)+','+str(0.15)]
    align_proc = subprocess.Popen(aligner_cmd, stdout=subprocess. PIPE,stderr=open("/dev/null", 'w'))
    
    # start sort the reads
    seq_sort_fname=open(out_fname+".sort",'w')
    seq_sort={}
    seq_sort_count = {}
    
    read_count=0

    for line in align_proc.stdout:
        if line.startswith('@'):
            continue
        #print line
        type, insert, cut_loc, read_seq = 'NA', 'NA', 'NA', 'NA'

        cols = line.strip().split()
        read_id, flag, ref_id, pos, mapQ, cigar_str = cols[:6]
        read_id=":".join(read_id.split(':')[3:7])
        read_count +=1
        flag, pos = int(flag), int(pos)
        pos-=1

        if pos <0:
            type = 'mutant'
            #continue
        if flag & 0x4 != 0:
            type = 'mutant'
            #continue

        read_seq, qual =cols[9:11]
        ref_id = ref_id.strip()

        if ref_id != '*':
            if len(read_seq) > ref_length[ref_id] + 10 or len(read_seq) < ref_length[ref_id] - 10:
                type = 'mutant'
                #continue

        AS,NM,MD = None, None, None
        for i in range(11, len(cols)):
            col = cols[i]
            if col.startswith('AS'):
                AS=int(col[5:])
            elif col.startswith('NM'):
                NM=int(col[5:])
            elif col.startswith('MD'):
                MD=col[5:]

        if NM > mm_cutoff:
            type = 'mutant'
            #continue

        # mutant:alignment failture or too many mismatch
        if type != 'NA':
            assert type == 'mutant'
            if type not in seq_sort_count:
                seq_sort_count[type] = 0.0
            seq_sort_count[type] +=1
            assert read_id not in seq_sort
            seq_sort[read_id]=[type, read_seq]
            print >> seq_sort_fname, "%s\t%s\t%s" % (read_id, type, read_seq)
            continue

        type = ref_id
        if type not in seq_sort_count:
            seq_sort_count[type] = 0.0
        seq_sort_count[type] +=1
        assert read_id not in seq_sort
        seq_sort[read_id]=[type, read_seq]

        print >> seq_sort_fname, "%s\t%s\t%s" % (read_id, type, read_seq)

    print seq_sort_count
    seq_sort_fname.close()

if __name__ == '__main__':
    parser = ArgumentParser(description='Sort the slide_seq data and extract valid data')
    parser.add_argument('-f1',
                        dest='read_fname1',
                        type=str,
                        help='one of pair-ends read fname1')
    parser.add_argument('-f2',
                        dest='read_fname2',
                        type=str,
                        help='one of pair-ends read fname2')
    parser.add_argument('-f',
                        dest='read_fname',
                        type=str,
                        help='read filename')
    parser.add_argument('-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence prefix filename')
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

    read_type=None
    
    if not args.read_fname:
        if not args.read_fname1 and not args.read_fname2:
            print >> sys.stderr, "Error: there is no input read file."
            sys.exit(1)
        elif not (args.read_fname1 and args.read_fname2):
            print >> sys.stderr, "Error: one of read-pair is missing."
            sys.exit(1)
        else:
            read_type='pair'
            
    else:
        if args.read_fname1 or args.read_fname2:
            print >> sys.stderr, "Error: too much inputs."
            sys.exit(1)
        else:
            read_type='single'

    assert read_type


    
    if read_type == 'pair':
        # set combined fastq filename
        def common (str1, str2):
            N=min(len(str1), len(str2))
            name=''
            for i in range(N):
                if str1[i] == str2[i]:
                    name += str1[i]
                else:
                    break
            name=name.strip()
            if not name:
                return 'out'
            return name

        common_name = common(args.read_fname1.rsplit('.')[0], args.read_fname2.rsplit('.')[0])
        #read_fname = common_name
        read_fname = common_name  + '.combined.fastq'

        
        # combine pair-end reads
        flash_cmd = ["flash", '-c', args.read_fname1, args.read_fname2, '-m2']
        flash_proc = subprocess.Popen(flash_cmd, stdout=open(read_fname,'w'), stderr=open("/dev/null", 'w'))
        flash_proc.communicate()
        
    else:
        assert read_type == 'single'
        read_fname = args.read_fname
    
    if not args.ref_fname:
        print >> sys.stderr, "Error: there is no reference file input."
        sys.exit(1)

    ref_seq, ref_length = {}, {}
    for line in open(args.ref_fname + '.ref'):
        line = line.strip()
        if line.startswith('>'):
            ref_id = line[1:].split()[0]
            continue
        if line:
            assert ref_id not in ref_seq
            assert ref_id not in ref_length
            ref_seq[ref_id] = line
            ref_length[ref_id] = len(line)
    
    
    if not args.out_fname:
        out_fname = read_fname.rsplit('.')[0]
    else:
        out_fname = args.out_fname
    
    sort_seqs (read_fname,
               args.ref_fname,
               out_fname,
               args.mm_cutoff,
               ref_length)
