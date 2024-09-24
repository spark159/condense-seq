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

# read titration file
def read_titration (fname, bg=False):
    all_fracs = {}
    tnum_tfrac = {}
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        cols = line.strip().split()
        try:
            tnum = int(cols[-1])
        except:
            continue
        total_frac = float(cols[-3])
        assert tnum not in tnum_tfrac
        tnum_tfrac[tnum] = total_frac
    return tnum_tfrac

# metric for condensabiltiy
def get_score (test, control, check_control=False):
    if check_control and control <= 0:
        return "NA"
    rcount = float(test)/float(control)
    return -math.log(rcount)


def NCP_score (fnames_list,
               control_fnames,
               numc_choice,
               tfname,
               tnums,
               chr_list,
               out_fname):

    # combining all file list
    all_fnames_list = copy.deepcopy(fnames_list)
    all_fnames_list.append(control_fnames)

    # get total coverage/counts for normalization
    print >> sys.stderr, "summing data for normalization"
    total_covs_list = []
    for fnames in all_fnames_list:
        total_covs = []
        for fname in fnames:            
            print >> sys.stderr, "\t reading %s" % (fname.rsplit('/', 1)[-1])
                            
            First = True
            for line in open(fname):
                line = line.strip()
                if not line:
                    continue
                cols = line.split()
                if First:
                    # set data range
                    if cols[2] == "PhysicalPosition":
                        data_type = 'point'
                        col_st = 3
                        col_ed = len(cols)
                    else:
                        assert cols[2] == 'Start'
                        assert cols[3] == 'End'
                        data_type = 'binned'
                        col_st = 4
                        try:
                            col_ed = cols.index('GCcontent')
                        except:
                            col_ed = len(cols)
                    
                    if len(total_covs) == 0:
                        total_covs = [0 for i in range(col_st, col_ed)]
                    else:
                        # all files has same data range
                        assert len(total_covs) == col_ed - col_st
                    First = False
                    continue

                # non target chromosome
                chr_name = cols[1]
                if chr_list and chr_name not in chr_list:
                    continue

                for j in range(col_st, col_ed):
                    cov = float(cols[j]) + 1 # add pesudo count
                    total_covs[j-col_st] += cov

        total_covs_list.append(total_covs)
        print >> sys.stderr, ""

    control_covs = total_covs_list.pop()
    assert len(total_covs_list) == len(fnames_list)
    
    # compute condensabilty scores    
    if numc_choice:
        # read titration file
        tnum_tfrac = read_titration (tfname)
    
    print >> sys.stderr, "computing scores"
    for k in range(len(fnames_list)):
        fnames = fnames_list[k]
        total_covs = total_covs_list[k]

        if numc_choice:
            tfrac = tnum_tfrac[tnums[k]]
            offset = get_score(tfrac, 1)
        else:
            offset = 0
        
        for fname, control_fname in zip(fnames, control_fnames):
            print >> sys.stderr, "processing %s" % (fname.rsplit('/', 1)[-1])

            test_f = open(fname)
            control_f = open(control_fname)

            out_fname = fname.rsplit('_', 1)[0] + '_score.cn'

            f = open(out_fname, 'w')

            First = True
            test_EOF, control_EOF = False, False
            while True:
                # read test file line by line but skipping empty line
                while True:
                    test_line = test_f.readline()
                    if not test_line:
                        test_EOF = True
                        break
                    test_line = test_line.strip()
                    if test_line:
                        break

                # read control file line by line but skipping empty line
                while True:
                    control_line = control_f.readline()
                    if not control_line:
                        control_EOF = True
                        break
                    control_line = control_line.strip()
                    if control_line:
                        break

                # break out the loop if one of files reach EOF
                if test_EOF or control_EOF:
                    break                

                test_cols = test_line.split()
                control_cols = control_line.split()

                if First:
                    if test_cols[2] == "PhysicalPosition":
                        data_type = 'point'
                        col_st = 3
                        col_ed = len(test_cols)
                    else:
                        assert test_cols[2] == 'Start'
                        assert test_cols[3] == 'End'
                        data_type = 'binned'
                        col_st = 4
                        try:
                            col_ed = test_cols.index('GCcontent')
                        except:
                            col_ed = len(test_cols)
                    assert test_cols[:col_st] == control_cols[:col_st] # same field
                    print >> f, test_line
                    First = False
                    continue

                assert test_cols[:col_st] == control_cols[:col_st] # same range

                # non target chromosome
                chr_name = test_cols[1]
                if chr_list and chr_name not in chr_list:
                    continue

                s = []
                s += test_cols[:col_st]
                for i in range(col_st, col_ed):
                    count = float(test_cols[i]) + 1 # add pseudo count
                    control = float(control_cols[i]) + 1 # add pseudo count
                    ncount = float(count)/total_covs[i-col_st]
                    ncontrol = float(control)/control_covs[i-col_st]
                    score = get_score(ncount, ncontrol)
                    score += offset # offset correction
                    s.append(str(round(score, 5)))

                print >> f, '\t'.join(s)

            assert test_EOF == control_EOF # same line number

            f.close()

        print >> sys.stderr, ""    
                    
    print >> sys.stderr, "Done"


if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = ArgumentParser(description='Get condensability scores')
    parser.add_argument(metavar='-f',
                        dest="fnames_list",
                        type=str,
                        nargs='+',
                        action='append',
                        help='Ncov/bin/Bdata.cn file list')
    parser.add_argument('-i',
                        dest="control_fnames",
                        type=str,
                        nargs='+',
                        help='input control files (in same order of data files)')
    parser.add_argument('-t',
                        dest='tfname',
                        type=str,
                        help='titration filename')
    parser.add_argument('--tnum',
                        dest="tnums",
                        type=int,
                        nargs='+',
                        help='titration number for each data (in same order of data files)')
    parser.add_argument('--numc',
                        dest="numc_choice",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='correct scores using titration file')    
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

    
    # check the control filenames
    for fnames in args.fnames_list:
        if len(fnames) != len(args.control_fnames):
            print >> sys.stderr, "Error: mismatch of data and control files"
            sys.exit(1)

            
    # check number correction option
    if args.numc_choice:
        numc_choice = True
        if not args.tnums or not args.tfname:
            print >> sys.stderr, "No titration information"
            print >> sys.stderr, "starting score computation without correction"
        # check tnums for filenames
        elif len(args.fnames_list) != len(args.tnums):
            print >> sys.stderr, "Error: mismatch of file and titration number."
            sys.exit(1)
        else:
            print >> sys.stderr, "starting score computation with correction"
            
    else:
        print >> sys.stderr, "starting score computation without correction"
        numc_choice = False
            
            
    # list target chromosomes
    chr_list = []
    if not args.chr_list:
        chr_list = None
    else:
        chr_list = sorted(args.chr_list, cmp=chr_cmp)

        
    NCP_score (args.fnames_list,
               args.control_fnames,
               numc_choice,
               args.tfname,
               args.tnums,
               chr_list,
               args.out_fname
               )
