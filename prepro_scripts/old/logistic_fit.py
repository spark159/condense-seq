import glob
import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import copy
import random
import functools
import warnings
import numpy as np
from scipy.optimize import differential_evolution

# 4 parameter logistic function
def fourPL (x, a, b, c, d):
    y = d + (a-d)/(1.0 + (x/float(c))**b)
    return (y)

# compute CP value from 4PL data
def get_CP (top, hill, chalf, bottom, percent):
    surv_frac = 1 - percent/100.0
    CP = chalf*(((float(top-bottom)/(surv_frac-bottom)) - 1)**(1.0/float(hill)))
    return CP

# objective function to optimize parameters
def obj_func (parms, func, X, Y):
    warnings.filterwarnings("ignore")
    Y_pred = func(X, *parms)
    return np.sum((Y_pred - Y)**2)

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

def logistic_fit (fnames,
                  tfname,
                  tnums,
                  min_rsq,
                  min_top,
                  max_top,
                  min_bottom,
                  max_bottom,
                  min_hill,
                  max_hill,
                  min_chalf,
                  max_chalf,
                  chr_list,
                  graph_option,
                  out_fname):
    
    # read titration file
    tnum_conc, tnum_tfrac = read_titration (tfname)

    # covert titration number to physical concentration
    concs = [tnum_conc[tnum] for tnum in tnums]

    # check input titration point
    input_index = tnums.index(0) # input data column index
    assert concs[input_index] == 0 # conc is zero at input

    # set C-half bound values
    if not min_chalf:
        min_chalf = min(concs)
    if not max_chalf:
        max_chalf = max(concs)

    # read files and compute C-half
    rsq_list = []
    total_count, fail_count = 0, 0
    
    f = open(out_fname + '_4PL.cn', 'w')
    
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
                s += ['Top',
                      'Hill',
                      'C-half',
                      'Bottom',
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

            X = np.asarray(X)
            Y = np.asarray(Y)

            # set boundary values for four parameters
            bounds = [(min_top, max_top),
                      (min_hill, max_hill),
                      (min_chalf, max_chalf),
                      (min_bottom, max_bottom)]

            result = differential_evolution(obj_func,
                                            args=(fourPL, X, Y),
                                            bounds=bounds,
                                            seed=3)

            p_opt = result.x
            success = result.success

            # check fitting quality
            if success:
                residuals = np.asarray(Y)- fourPL(X, *p_opt)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((np.asarray(Y)-np.mean(Y))**2)
                r_squared = 1 - (ss_res / ss_tot)

                if r_squared < min_rsq:
                    success = False
                    pass

            if success:
                top, hill, chalf, bottom = p_opt
                rsq_list.append(r_squared)

                # temporal
                CP80 = get_CP (top, hill, chalf, bottom, 80)

            else:
                top, hill, chalf, bottom = 'NA', 'NA', 'NA', 'NA'
                r_squared = 'NA'
                fail_count +=1

            if graph_option and success and abs(np.log10(CP80)) > 3 :

                print (top, hill, chalf, bottom)
                print (r_squared)
                print (CP80)
                
                fig = plt.figure()
                plt.plot(X, Y, '.', markersize=10, alpha=0.2)

                if success:
                    X_pred = np.linspace(min(X), max(X), 1000)
                    Y_pred = fourPL(X_pred, *p_opt)
                    plt.plot(X_pred, Y_pred, 'k-', alpha=0.2)
                    plt.axhline(y=bottom + (top-bottom)*0.5, linestyle='--', color='b')
                    plt.axvline(x=chalf, linestyle='--', color='r')

                plt.xlabel("Concentration")
                plt.ylabel("Soluble fractin")
                #if agent in ['HP1a']:
                #    plt.xscale('log', basex=2)
                #elif agent in ['sp', 'spd', 'CoH']:
                #    plt.xscale('log', basex=10)
                #plt.xscale('log', base=10)
                #plt.title("%s" % (chalf))
                plt.title("%s" % (r_squared))
                #plt.title("%s binID:%d" % (chr, binID))
                #plt.xlim([min(X)-1, max(X)+1])
                #plt.ylim([-0.2, 1.2])
                plt.show()
                plt.close()

            s = cols[:col_st]
            s += [str(top),
                  str(hill),
                  str(chalf),
                  str(bottom),
                  str(r_squared)]
        
            print('\t'.join(s), end='\n', file=f)
            total_count +=1

    f.close()

    # summarize output
    print("fitting failure %d/%d" % (fail_count, total_count), file=sys.stderr)
    print("fitting failure %.2f %%" % (100*float(fail_count)/total_count), file=sys.stderr)
    print("Median R-squared of fitting %.2f" % (np.median(rsq_list)), file=sys.stderr)
    print("Done", file=sys.stderr)
    

if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = ArgumentParser(description='fitting condense-seq data with four parameter logistic function')
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
    parser.add_argument('--min_top',
                        dest="min_top",
                        type=float,
                        default=1.0,
                        help='lower bound of Top parameter in 4PL model')
    parser.add_argument('--max_top',
                        dest="max_top",
                        type=float,
                        default=1.0,
                        help='upper bound of Top parameter in 4PL model')
    parser.add_argument('--min_bottom',
                        dest="min_bottom",
                        type=float,
                        default=0.0,
                        help='lower bound of Bottom parameter in 4PL model')
    parser.add_argument('--max_bottom',
                        dest="max_bottom",
                        type=float,
                        default=0.0,
                        help='upper bound of Bottom parameter in 4PL model')
    parser.add_argument('--min_hill',
                        dest="min_hill",
                        type=float,
                        default=0.0,
                        help='lower bound of Hill parameter in 4PL model')
    parser.add_argument('--max_hill',
                        dest="max_hill",
                        type=float,
                        default=100.0,
                        help='upper bound of Hill parameter in 4PL model')
    parser.add_argument('--min_chalf',
                        dest="min_chalf",
                        type=float,
                        help='lower bound of C-half parameter in 4PL model')
    parser.add_argument('--max_chalf',
                        dest="max_chalf",
                        type=float,
                        help='upper bound of C-half parameter in 4PL model')
    parser.add_argument('--chr',
                        dest="chr_list",
                        type=str,
                        nargs='+',
                        help='tagert chromosome list')
    parser.add_argument('--graph',
                        dest="graph_option",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='plot the grpahs fitting the data')
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
        chr_list = sorted(args.chr_list, key=functools.cmp_to_key(chr_cmp))

    if args.graph_option:
        import matplotlib.pyplot as plt

    logistic_fit (args.fnames,
                  args.tfname,
                  args.tnums,
                  args.min_rsq,
                  args.min_top,
                  args.max_top,
                  args.min_bottom,
                  args.max_bottom,
                  args.min_hill,
                  args.max_hill,
                  args.min_chalf,
                  args.max_chalf,
                  chr_list,
                  args.graph_option,
                  args.out_fname
                  )
