import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import copy
import gzip

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
        reading_file = gzip.open(fname, 'rb')
    else:
        reading_file = open(fname, 'r')
    return reading_file

def standardize (fnames,
                 chr_list,
                 colnums,
                 out_fname):

    count_list = None
    mean_list = None
    std_list = None
    fieldnames = None

    ## read all files 2 times to get statistics
    ## 1st: get mean, 2nd: get std
    for i in range(2):
        if i == 0:
            print >> sys.stderr, "summing up data to compute means"
        elif i == 1:
            print >> sys.stderr, "summing up data to compute stds"
            
        for fname in fnames:            
            print >> sys.stderr, "\t reading %s" % (fname.rsplit('/', 1)[-1])

            First = True
            for line in gzopen(fname):
                line = line.strip()
                if not line:
                    continue
                cols = line.split()

                if First:
                    # check field names
                    if fieldnames == None:
                        fieldnames = cols
                    else:
                        # all files has same field names
                        assert cols == fieldnames

                    # set data range
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

                    if colnums == None:
                        colnums = range(col_st, col_ed)

                    if i == 0 and mean_list == None:
                        count_list = [0] * len(colnums)
                        mean_list = [0] * len(colnums)

                    if i == 1 and std_list == None:
                        std_list = [0] * len(colnums)

                    First = False
                    continue

                # non target chromosome
                chr_name = cols[0]
                if chr_list and chr_name not in chr_list:
                    continue

                # sum up the data
                for k in range(len(colnums)):
                    colnum = colnums[k]
                    try:
                        value = float(cols[colnum])
                    except:
                        continue

                    if i == 0:
                        count_list[k] +=1
                        mean_list[k] += value
                    elif i == 1:
                        mean = mean_list[k]
                        std = (value - mean)**2
                        std_list[k] += std

        # take average of sums at last
        for k in range(len(colnums)):
            count = count_list[k]
            if count <= 0:
                continue
            if i == 0:            
                mean_list[k] = float(mean_list[k]) / count
            elif i == 1:
                std_list[k] = math.sqrt(float(std_list[k]) / count)

    ##compute z-score and writing down output files
    print >> sys.stderr, "data standardization and writing output"
    for fname in fnames:
        out_fname = fname.rsplit('_', 1)[0] + '_zscore.gtab.gz'
        print >> sys.stderr, "\t writing %s" % (out_fname.rsplit('/', 1)[-1])
        
        f = gzip.open(out_fname, 'wb')
        
        First = True
        for line in gzopen(fname):
            line = line.strip()
            if not line:
                continue
            cols = line.split()

            if First:
                # set data range
                if cols[1] == "Position":
                    data_type = 'point'
                    col_st = 2
                else:
                    assert cols[1] == 'Start'
                    assert cols[2] == 'End'
                    data_type = 'binned'
                    col_st = 3
                    
                row = cols[:col_st]
                for colnum in colnums:
                    row.append(cols[colnum])
                print >> f, '\t'.join(row)
                First = False
                continue

            # non target chromosome
            chr_name = cols[0]
            if chr_list and chr_name not in chr_list:
                continue

            # compute z-score and writing to output
            row = cols[:col_st]
            for k in range(len(colnums)):
                colnum = colnums[k]
                try:
                    value = float(cols[colnum])
                    mean = mean_list[k]
                    std = std_list[k]
                    zscore = round(float(value-mean)/std, 5)
                except: # can't compute
                    zscore = 'NA'
                row.append(str(zscore))

            print >> f, '\t'.join(row)

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

    parser = ArgumentParser(description='standardize values in gtab file')
    parser.add_argument(metavar='-f',
                        dest="fnames",
                        type=str,
                        nargs='+',
                        help='gtab file list')
    parser.add_argument('--chr',
                        dest="chr_list",
                        type=str,
                        nargs='+',
                        help='tagert chromosome list')
    parser.add_argument('-c',
                        dest="colnums",
                        type=str,
                        nargs='+',
                        help='picked column numbers of files for standardization')
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()

    # list target chromosomes
    chr_list = []
    if not args.chr_list:
        chr_list = None
    else:
        chr_list = sorted(args.chr_list, cmp=chr_cmp)

    # list target column numbers
    if not args.colnums:
        colnums = None
    else:
        colnums = [int(colnum) for colnum in colnums.split(',')]

    standardize (args.fnames,
                 chr_list,
                 colnums,
                 args.out_fname
                 )
