import glob
import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import copy
import random
import functools
import warnings
import gzip
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution

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
        reading_file = gzip.open(fname, 'rt')
    else:
        reading_file = open(fname, 'r')
    return reading_file

# 4-parameter logistic function (sigmoid type)
def sigmoid_func (x, top, rate, chalf, bottom):
    y = bottom + float(top-bottom)/(1+np.exp(rate*(x-chalf)))
    return y

# 4-parameter logistic function (Hill type)
def hill_func (x, top, rate, chalf, bottom):
    y = bottom + float(top-bottom)/(1.0 + (x/float(chalf))**rate)
    return y

# objective function to optimize parameters
def obj_func (parms, func, X, Y):
    warnings.filterwarnings("ignore")
    Y_pred = func(X, *parms)
    return np.sum((Y_pred - Y)**2)

# compute CP value from 4PL data
def get_CP (top, hill, chalf, bottom, percent):
    surv_frac = 1 - percent/100.0
    CP = chalf*(((float(top-bottom)/(surv_frac-bottom)) - 1)**(1.0/float(hill)))
    return CP

# read titration file
def read_titration (fname):
    tnum_conc = {}
    tnum_frac = {}
    for line in gzopen(fname):
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

def logistic_fit (fnames,
                  tfname,
                  tnums,
                  model,
                  method,
                  min_rsq,
                  min_top,
                  max_top,
                  min_bottom,
                  max_bottom,
                  min_rate,
                  max_rate,
                  min_chalf,
                  max_chalf,
                  chr_list,
                  graph_option,
                  out_fname):
    
    # read titration file
    tnum_conc, tnum_tfrac = read_titration (tfname)

    # covert titration number to physical concentration
    concs = [tnum_conc[tnum] for tnum in tnums]

    # set indices of data by the order of concentration
    conc_idx = sorted([(concs[i], i) for i in range(len(concs))])
    idx_list = [idx for conc, idx in conc_idx]
    
    # check input titration point    
    input_index = idx_list[0] # input data column index
    assert concs[input_index] == 0 # conc is zero at input

    # read files and start logistic regression
    rsq_list = []
    total_count, fail_count = 0, 0
    
    f = gzip.open(out_fname + '_4PL.gtab.gz', 'wt')
    
    for fname in fnames:
        print("Processing %s" % (fname.rsplit('/', 1)[-1]),
              file=sys.stderr)
        First = True
        for line in gzopen(fname):
            line = line.strip()
            if not line:
                continue
            cols = line.split()
            
            if First:
                # find data type and range
                if cols[1] == "Position":
                    data_type = 'point'
                    col_st = 2
                    col_ed = len(cols)
                else:
                    assert cols[1] == 'Start'
                    assert cols[2] == 'End'
                    data_type = 'binned'
                    col_st = 3
                    try:
                        col_ed = cols.index('GCcontent')
                    except:
                        col_ed = len(cols)

                assert col_ed - col_st == len(tnums)  # data col len == titration num

                s = cols[:col_st]
                s += ['Model',
                      'Method',
                      'Top',
                      'Rate',
                      'C-half',
                      'Bottom',
                      'R-squared']
                
                print('\t'.join(s), end='\n', file=f) # write header
                First = False
                continue

            # non target chromosome
            chr_name = cols[0]
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
            
            # organize titration data for fitting
            X = np.asarray([concs[idx] for idx in idx_list])
            Y = np.asarray([fracs[idx] for idx in idx_list])

            # set guess and boundary of parameters
            top_guess = max(Y)
            bottom_guess = min(Y)
            chalf_guess = np.mean(X)
            rate_guess = 1.0

            if min_top == None:
                top_min = top_guess * (1-0.01)
            else:
                top_min = min_top

            if max_top == None:
                top_max = top_guess * (1+0.01)
            else:
                top_max = max_top

            if min_bottom == None:
                bottom_min = bottom_guess * (1-0.01)
            else:
                bottom_min = min_bottom

            if max_bottom == None:
                bottom_max = bottom_guess * (1+0.01)
            else:
                bottom_max = max_bottom

            if min_chalf == None:
                chalf_min = min(X)
            else:
                chalf_min = min_chalf

            if max_chalf == None:
                chalf_max = max(X)
            else:
                chalf_max = max_chalf

            if min_rate == None:
                rate_min = 0.0
            else:
                rate_min = min_rate

            if max_rate == None:
                rate_max = 100.0
            else:
                rate_max = max_rate
                
            # fitting the data with a logistic function
            if method == 'curve_fit':
                # set initial guess of parameters
                p0 = [top_guess, rate_guess, chalf_guess, bottom_guess]
                # set boundary values of parameters
                bounds = ([top_min, rate_min, chalf_min, bottom_min],
                          [top_max, rate_max, chalf_max, bottom_max])

                try:
                    if model == 'sigmoid':                        
                        p_opt, p_cov = curve_fit(sigmoid_func,
                                                 X,
                                                 Y,
                                                 p0,
                                                 bounds = bounds,
                                                 method='dogbox')
                    elif model == 'hill':
                        p_opt, p_cov = curve_fit(hill_func,
                                                 X,
                                                 Y,
                                                 p0,
                                                 bounds = bounds,
                                                 method='dogbox')
                    success = True
                except:
                    success = False

            elif method == 'evolution':
                # set boundary values of parameters
                bounds = [(top_min, top_max),
                          (rate_min, rate_max),
                          (chalf_min, chalf_max),
                          (bottom_min, bottom_max)]

                if model == 'sigmoid':
                    result = differential_evolution(obj_func,
                                                    args=(sigmoid_func, X, Y),
                                                    bounds=bounds,
                                                    seed=3)
                elif model == 'hill':
                    result = differential_evolution(obj_func,
                                                    args=(hill_func, X, Y),
                                                    bounds=bounds,
                                                    seed=3)

                p_opt = result.x
                success = result.success

            # check fitting quality
            if success:
                if model == 'sigmoid':
                    residuals = np.asarray(Y)- sigmoid_func(X, *p_opt)
                elif model == 'hill':
                    residuals = np.asarray(Y)- hill_func(X, *p_opt)

                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((np.asarray(Y)-np.mean(Y))**2)
                r_squared = 1 - (ss_res / ss_tot)

                # too poor fitting
                if r_squared < min_rsq:
                    success = False
                    pass

            if success:
                top, rate, chalf, bottom = p_opt
                rsq_list.append(r_squared)

            else:
                top, rate, chalf, bottom = 'NA', 'NA', 'NA', 'NA'
                r_squared = 'NA'
                fail_count +=1

            if graph_option:
            #if graph_option and not success:
            #if graph_option and success and rate>100:
            #if graph_option and success and r_squared < 0.7:

                print (top, rate, chalf, bottom)
                print (r_squared)
                #print (CP80)
                
                fig = plt.figure()
                plt.plot(X, Y, '.', markersize=10, alpha=0.2)

                if success:
                    X_pred = np.linspace(min(X), max(X), 1000)
                    if model == 'sigmoid':
                        Y_pred = sigmoid_func(X_pred, *p_opt)
                    elif model == 'hill':
                        Y_pred = hill_func(X_pred, *p_opt)
                    plt.plot(X_pred, Y_pred, 'k-', alpha=0.2)
                    plt.axhline(y=bottom + (top-bottom)*0.5, linestyle='--', color='b')
                    plt.axvline(x=chalf, linestyle='--', color='r')

                plt.xlabel("Concentration")
                plt.ylabel("Soluble fractin")
                plt.title("%s" % (r_squared))
                plt.show()
                plt.close()

            s = cols[:col_st]
            s += [model,
                  method,
                  str(top),
                  str(rate),
                  str(chalf),
                  str(bottom),
                  str(r_squared)]
        
            print('\t'.join(s), end='\n', file=f)
            total_count +=1

    f.close()

    # summarize output
    print("fitting failure %d/%d (%.2f %%)"
          % (fail_count,
             total_count,
             100*float(fail_count)/total_count),
          file=sys.stderr)
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
                        help='concatenated num.gtab file list')
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
    parser.add_argument('--model',
                        dest="model",
                        type=str,
                        default='sigmoid',
                        help='logistic model for fitting data (sigmoid or hill)')
    parser.add_argument('--method',
                        dest="method",
                        type=str,
                        default='evolution',
                        help='logistic regression method (curve_fit or evolution)')    
    parser.add_argument('--min_rsq',
                        dest="min_rsq",
                        type=float,
                        default=0.5,
                        help='minimum R-squared value for fitting quality')
    parser.add_argument('--min_top',
                        dest="min_top",
                        type=float,
                        help='lower bound of Top parameter in 4PL model')
    parser.add_argument('--max_top',
                        dest="max_top",
                        type=float,
                        help='upper bound of Top parameter in 4PL model')
    parser.add_argument('--min_bottom',
                        dest="min_bottom",
                        type=float,
                        help='lower bound of Bottom parameter in 4PL model')
    parser.add_argument('--max_bottom',
                        dest="max_bottom",
                        type=float,
                        help='upper bound of Bottom parameter in 4PL model')
    parser.add_argument('--min_rate',
                        dest="min_rate",
                        type=float,
                        help='lower bound of Rate parameter in 4PL model')
    parser.add_argument('--max_rate',
                        dest="max_rate",
                        type=float,
                        help='upper bound of Rate parameter in 4PL model')
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

    # set logistic model
    model = args.model.lower()

    # set regression method
    method = args.method.lower()

    logistic_fit (args.fnames,
                  args.tfname,
                  args.tnums,
                  model,
                  method,
                  args.min_rsq,
                  args.min_top,
                  args.max_top,
                  args.min_bottom,
                  args.max_bottom,
                  args.min_rate,
                  args.max_rate,
                  args.min_chalf,
                  args.max_chalf,
                  chr_list,
                  args.graph_option,
                  args.out_fname
                  )
