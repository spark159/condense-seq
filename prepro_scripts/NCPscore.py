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

# metric for condensabiltiy
def get_score (test, control, check_control=False):
    if check_control and control <= 0:
        return "NA"
    test +=1
    control +=1
    rcount = float(test)/float(control)
    return -math.log(rcount)


def NCP_score (fnames_list,
               input_fnames,
               concat,
               labels,
               out_fnames):

    for fnames in fnames_list:
        for fname, input_fname in zip(fnames, input_fnames):
            print >> sys.stderr, "processing %s" % (fname.rsplit('/', 1)[-1])
            
            test_f = open(fname)
            input_f = open(input_fname)

            output_fname = fname.rsplit('_num.cn', 1)[0] + '_score.cn'
            f = open(output_fname, 'w')

            First = True
            test_EOF, input_EOF = False, False
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

                # read input file line by line but skipping empty line
                while True:
                    input_line = input_f.readline()
                    if not input_line:
                        input_EOF = True
                        break
                    input_line = input_line.strip()
                    if input_line:
                        break

                # break out the loop if one of files reach EOF
                if test_EOF or input_EOF:
                    break                

                test_cols = test_line.split()
                input_cols = input_line.split()

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
                    assert test_cols[:col_st] == input_cols[:col_st] # same field
                    print >> f, test_line
                    First = False
                    continue

                assert test_cols[:col_st] == input_cols[:col_st] # same range

                s = []
                s += test_cols[:col_st]
                for i in range(col_st, col_ed):
                    count = int(test_cols[i])
                    control = int(input_cols[i])
                    score = get_score(count, control)
                    score = round(score, 5)
                    s.append(str(score))

                print >> f, '\t'.join(s)

            assert test_EOF == input_EOF # same line number

            f.close()

    # concatenating files option
    if concat:        
        assert True

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
                        help='num.cn file list')
    parser.add_argument('--input',
                        dest="input_fnames",
                        type=str,
                        nargs='+',
                        help='input control num.cn files (in same order of data files)')
    parser.add_argument('--concat',
                        dest="concat",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='concatenate output files')
    parser.add_argument('--label',
                        dest='labels',
                        type=str,
                        nargs='+',
                        help='field label names')
    parser.add_argument('-o',
                        dest='out_fnames',
                        type=str,
                        nargs='+',
                        help='output prefix filenames')
    
    args = parser.parse_args()


    # check the input is correct
    for fnames in args.fnames_list:
        if len(fnames) != len(args.input_fnames):
            print >> sys.stderr, "Error: mismatch of data and input files"
            sys.exit(1)

    if args.concat:
        if len(args.out_fnames) != len(args.input_fnames):
            print >> sys.stderr, "Error: incorrect output filenames"
            sys.exit(1)
    
    NCP_score (args.fnames_list,
               args.input_fnames,
               args.concat,
               args.labels,
               args.out_fnames
               )
