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

    if numc_choice:
        # read titration file
        tnum_tfrac = read_titration (tfname)

    # first-loop: get total coverage/counts for normalization
    # second-loop: compute condensabilty scores
    for u in range(2):

        if u == 0:
            print >> sys.stderr, "summing data for normalization"
            total_covs_list = []
        else:
            print >> sys.stderr, "computing scores"

        for k in range(len(fnames_list)):
            fnames = fnames_list[k]

            if u == 0:
                total_covs = []
                if k == 0:
                    control_covs = []
            else:
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

                if u == 1:
                    out_fname = fname.rsplit('_', 1)[0] + '_score.cn'
                    f = open(out_fname, 'w')

                test_read, control_read = True, True
                test_EOF, control_EOF = False, False
                First = True
                ID = 0
                while True:
                    # read test file line by line but skipping empty line
                    while test_read and not test_EOF:
                        test_line = test_f.readline()

                        if not test_line:
                            test_EOF = True
                            break

                        test_line = test_line.strip()
                        if test_line:
                            test_cols = test_line.split()

                            if not First:
                                test_chr = test_cols[1]
                                test_pos = test_cols[2:col_st]

                                # skip non-target chromosome
                                if chr_list and test_chr not in chr_list:
                                    continue

                                try:
                                    test_chr = int(test_chr[3:])
                                except:
                                    test_chr = test_chr[3:]

                                test_pos = [int(value) for value in test_pos]
                                test_pt = [test_chr] + test_pos
                                test_pt = tuple(test_pt)    

                            break

                    # read control file line by line but skipping empty line
                    while control_read and not control_EOF:
                        control_line = control_f.readline()

                        if not control_line:
                            control_EOF = True
                            break

                        control_line = control_line.strip()
                        if control_line:
                            control_cols = control_line.split()

                            if not First:
                                control_chr = control_cols[1]
                                control_pos = control_cols[2:col_st]

                                # skip non-target chromosome
                                if chr_list and control_chr not in chr_list:
                                    continue

                                try:
                                    control_chr = int(control_chr[3:])
                                except:
                                    control_chr = control_chr[3:]

                                control_pos = [int(value) for value in control_pos]
                                control_pt = [control_chr] + control_pos
                                control_pt = tuple(control_pt)

                            break

                    ## break out the loop if control files reach EOF
                    #if control_EOF:
                    #    break

                    # break out the loop if both files reach EOF
                    if test_EOF and control_EOF:
                        break

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

                        if u == 0:
                            if len(total_covs) == 0:
                                total_covs = [0] * (col_ed - col_st)
                                
                                if k == 0 and len(control_covs) == 0:
                                    control_covs = [0] * (col_ed - col_st)

                            else:
                                # all files has same data range
                                assert len(total_covs) == col_ed - col_st
                        else:
                            print >> f, test_line

                        First = False
                        continue

                    if test_EOF: # test data is empty
                        assert not control_EOF
                        data_pt = control_pt
                        test_data = [0] * (col_ed - col_st)
                        control_data = control_cols[col_st:col_ed]
                        test_read = False
                        control_read = True

                    elif control_EOF: # skip when control data is empty
                        assert not test_EOF
                        #data_pt = test_pt
                        #txest_data = test_cols[col_st:col_ed]
                        #control_data = [0] * (col_ed - col_st)
                        test_read = True
                        control_read = False
                        continue

                    else:
                        assert not test_EOF
                        assert not control_EOF

                        # skip when control data is empty
                        if test_pt < control_pt:
                            #data_pt = test_pt
                            #test_data = test_cols[col_st:col_ed]
                            #control_data = [0] * (col_ed - col_st)
                            test_read = True
                            control_read = False
                            continue

                        # test data is empty
                        if test_pt > control_pt:
                            data_pt = control_pt
                            test_data = [0] * (col_ed - col_st)
                            control_data = control_cols[col_st:col_ed]
                            test_read = False
                            control_read = True

                        # both test and control data are available
                        if test_pt == control_pt:
                            data_pt = test_pt
                            test_data = test_cols[col_st:col_ed]
                            control_data = control_cols[col_st:col_ed]
                            test_read = True
                            control_read = True

                    # skp when whole row of control is empty
                    if sum([float(count) for count in control_data]) <=0:
                        continue

                    # compare the test and control data to compute the score
                    s = []
                    for i in range(col_ed - col_st):
                        count = float(test_data[i])
                        control = float(control_data[i])

                        if control <= 0: # score is not defined when control data is empty
                            score = 'NA'
                            s.append(score)
                            continue

                        # add pseudo count
                        count+=1
                        control+=1

                        if u == 0:
                            total_covs[i] += count

                            if k == 0:
                                control_covs[i] += control

                        else:
                            # normalize by total
                            ncount = float(count)/total_covs[i]
                            ncontrol = float(control)/control_covs[i]

                            # get score
                            score = get_score(ncount, ncontrol)
                            score += offset # offset correction

                            s.append(str(round(score, 5)))

                    if u == 1:
                        row = [str(ID)]
                        row += ['chr' + str(data_pt[0])]
                        row += [str(value) for value in data_pt[1:]]
                        row += s
                        print >> f, '\t'.join(row)
                        ID +=1

                if u == 1:
                    f.close()

        if u == 0:
            total_covs_list.append(total_covs)

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
