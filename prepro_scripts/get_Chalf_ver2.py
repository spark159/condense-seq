import glob
import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import copy
import random
import numpy as np
import scipy
from scipy.optimize import curve_fit

# logistic function
def logistic_func(x, L ,x0, k):
    y = L / (1 + np.exp(k*(x-x0)))
    return (y)

# read titration file
def read_titration (fname):
    tnum_conc = {}
    tnum_frac = {}
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        conc, frac, tnum = cols[0], cols[7], cols[-1]
        try:
            tnum = int(tnum)
        except:
            continue
        conc = float(conc)
        frac = float(frac)
        tnum_conc[tnum] = conc
        tnum_frac[tnum] = frac
    return tnum_conc, tnum_frac

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

def get_Chalf (fnames,
               tfname,
               tnums,
               min_rsq,
               chr_list,
               out_fname):
    
    # read titration file
    tnum_conc, tnum_tfrac = read_titration (tfname)

    # covert titration number to physical concentration
    concs = [tnum_conc[tnum] for tnum in tnums]

    # check input titration point
    input_index = tnums.index(0) # input data column index
    assert concs[input_index] == 0 # conc is zero at input

    # read files and compute C-half
    rsq_list = []
    total_count, fail_count = 0, 0
    
    f = open(out_fname + '_Chalf.cn', 'w')
    
    for fname in fnames:
        print("Processing %s" % (fname.rsplit('/', 1)[-1]),
              file=sys.stderr)
        First = True
        for line in open(fname):
            line = line.strip()
            if not line:
                continue
            cols = line.split()
            
            if First:
                # find data type and range
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

                assert col_ed - col_st == len(tnums)  # data col len == titration num

                s = cols[:col_st]
                s += ['C-half',
                      'height',
                      'rate',
                      'R-squared']
                
                print('\t'.join(s), end='\n', file=f) # write header
                First = False
                continue

            # non target chromosome
            chr_name = cols[1]
            if chr_list and chr_name not in chr_list:
                continue

            # get molecule numbers
            nums = []
            for i in range(col_st, col_ed):
                num = int(cols[i])
                nums.append(num)

            input_num = nums[input_index] # input molecule number

            if input_num <=0: # skip if input has no molecules
                continue

            fracs = [float(num/input_num) for num in nums] # convert num to fraction
            
            # fitting the data with a logistic function
            X, Y = copy.deepcopy(concs), fracs

            try:
                p0 = [max(Y), np.median(X), 1]
                bounds = ([0.0, 0.0, 0.0], [max(Y)+max(Y)*0.1, np.inf, np.inf])
                
                popt, pcov = curve_fit(logistic_func,
                                       X,
                                       Y,
                                       p0,
                                       bounds = bounds,
                                       method='dogbox')

                residuals = np.asarray(Y)- logistic_func(X, *popt)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((np.asarray(Y)-np.mean(Y))**2)

                r_squared = 1 - (ss_res / ss_tot) # fitting quality
                height, Chalf, rate = popt[:3] # fitting parameters
            except:
                # fitting failure
                r_squared = 'NA'
                height, Chalf, rate = 'NA', 'NA', 'NA'
                fail_count +=1

            if r_squared != 'NA' and r_squared < min_rsq: # poor fitting
                r_squared = 'NA'
                height, Chalf, rate = 'NA', 'NA', 'NA'
                fail_count +=1

            if r_squared != 'NA':
                rsq_list.append(r_squared)

            s = cols[:col_st]
            s += [str(Chalf),
                  str(height),
                  str(rate),
                  str(r_squared)]
        
            print('\t'.join(s), end='\n', file=f)
            total_count +=1

    f.close()

    # summarize output
    print("fitting failure %d/%d" % (fail_count, total_count), file=sys.stderr)
    print("fitting failure %.2f %%" % (100*float(fail_count)/total_count), file=sys.stderr)
    print("Median R-squared of fitting %.2f" % (np.meidian(rsq_list)), file=sys.stderr)
    print("Done", file=sys.stderr)
    

if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = ArgumentParser(description='Compute the condensation point (C 1/2) by fitting data with a logistic function')
    parser.add_argument(metavar='-f',
                        dest="fnames",
                        type=str,
                        nargs='+',
                        help='concatenated num.cn file list')
    parser.add_argument('-t',
                        dest='tfname',
                        type=str,
                        help='titration filename')
    parser.add_argument('--tnum',
                        dest="tnums",
                        type=int,
                        nargs='+',
                        help='titration number of each columns in num data')
    parser.add_argument('--min_tnum',
                        dest="min_tnum",
                        type=int,
                        default=3,
                        help='minimum titration data number for fitting')
    parser.add_argument('--min_rsq',
                        dest="min_rsq",
                        type=float,
                        default=0.5,
                        help='minimum R-squared value for fitting quality')
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

    # check titration point of input control is included
    if 0 not in args.tnums:
        print("input control data should be included (titration number is zero)",
              file=sys.stderr)
        sys.exit(1)
        
    # check at least 3 titration points for fitting
    if len(args.tnums) < args.min_tnum:
        print("Need more data points for fitting",
              file=sys.stderr)
        sys.exit(1)
    
    # list target chromosomes
    chr_list = []
    if not args.chr_list:
        chr_list = None
    else:
        chr_list = sorted(args.chr_list, cmp=chr_cmp)
    
    get_Chalf (args.fnames,
               args.tfname,
               args.tnums,
               args.min_rsq,
               chr_list,
               args.out_fname
               )
