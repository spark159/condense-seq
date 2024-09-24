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

def NCP_number (tnum_fnames,
                tfname,
                input_mnum,
                chr_list,
                concat,
                out_fname):

    # read titration file
    tnum_tfrac = read_titration (tfname)

    # sort by titration number
    tnums = sorted(tnum_fnames.keys())

    for i in range(len(tnums)):
        tnum = tnums[i]
        print >> sys.stderr, "processing tnum %d" % (tnum) 

        mnum = input_mnum * tnum_tfrac[tnum] # molecule number of titraion

        # get total of all total coverage/counts
        total_covs = []
        fnames = tnum_fnames[tnum]
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

                for i in range(col_st, col_ed):
                    cov = float(cols[i])
                    total_covs[i-col_st] += cov

        # convert to molecule number
        for fname in fnames:
            print >> sys.stderr, "\t converting %s" % (fname.rsplit('/', 1)[-1])

            out_fname = fname.rsplit('_', 1)[0] + '_num.cn'
            f = open(out_fname, 'w')

            First = True
            for line in open(fname):
                line = line.strip()
                if not line:
                    continue
                cols = line.split()
                if First:
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

                    if data_type == 'point':
                        fields = ['SNP', 'Chromosome', 'PhysicalPosition']
                    elif data_type == 'binned':
                        fields = ['BinID', 'Chromosome', 'Start', 'End']                            
                    fields += cols[col_st:col_ed]
                    print >> f, '\t'.join(fields)
                    
                    First = False
                    continue

                # non target chromosome
                chr_name = cols[1]
                if chr_list and chr_name not in chr_list:
                    continue

                nums = []
                for i in range(col_st, col_ed):
                    cov = float(cols[i])
                    frac = float(cov)/total_covs[i-col_st]
                    num = int(round(mnum * frac))
                    nums.append(str(num))

                print >> f, '\t'.join(cols[:col_st] + nums)

            f.close()

    # To dos
    if concat:
        print >> sys.stderr, "Concatenating output files over titrations"

        fnum = len(tnum_fnames[tnums[0]])
        groups = []
        i = 0
        while i < fnum:
            group = []
            for tnum in tnums:
                fname = tnum_fnames[tum][i]
                if fname.endswith('_bin.cn'):
                    out_fname = fname.rsplit('_bin.cn', 1)[0] + '_num.cn'
                elif fname.endswith('_Ncov.cn'):
                    out_fname = fname.rsplit('_Ncov.cn', 1)[0] + '_num.cn'
                group.append(out_fname)
            groups.append(group)
            i +=1
                    
    print >> sys.stderr, "Done"


if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = ArgumentParser(description='Estimate physical number of molecules')
    parser.add_argument(metavar='-f',
                        dest="fnames_list",
                        type=str,
                        nargs='+',
                        action='append',
                        help='Ncov/bin/Bdata.cn file list')
    parser.add_argument('-t',
                        dest='tfname',
                        type=str,
                        help='titration filename')
    parser.add_argument('--tnum',
                        dest="tnums",
                        type=int,
                        nargs='+',
                        help='titration number of each data')
    parser.add_argument('--mscale',
                        dest="mscale",
                        default = 1,
                        type=int,
                        help='total molecule number scale of input')
    parser.add_argument('--chr',
                        dest="chr_list",
                        type=str,
                        nargs='+',
                        help='tagert chromosome list')
    parser.add_argument('--concat',
                        dest="concat",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='concatenate output files')
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()


    # map tnum to filenames
    if not args.tnums:
        tnums = range(len(args.fnames_list))
    else:
        if len(args.fnames_list) != len(args.tnums):
            print >> sys.stderr, "Error: mismatch of input file and titration number."
            sys.exit(1)
        tnums = args.tnums
        
    tnum_fnames = {}
    for tnum, fnames in zip(tnums, args.fnames_list):
        tnum_fnames[tnum] = fnames

    # map file to file over titrations
    if args.concat:
        for tnum, fnames in tnum_fnames.items():
            if len(fnames) != len(args.fnames_list[0]):
                print >> sys.stderr, "Error: file number is not same over titration."
                sys.exit(1)

    # list target chromosomes
    chr_list = []
    if not args.chr_list:
        chr_list = None
    else:
        chr_list = sorted(args.chr_list, cmp=chr_cmp)

    input_mnum = args.mscale * 1.6*(10**12)
    
    NCP_number (tnum_fnames,
                args.tfname,
                input_mnum,
                chr_list,
                args.concat,
                args.out_fname
                )
