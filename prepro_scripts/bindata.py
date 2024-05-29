import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import copy
import time

# chromosome comparison sort
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

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

def Bin_data (data_fname,
              chr_size,
              bin_size,
              bin_step,
              bin_value,
              skip_zero,
              chr_list,
              out_fname):

    # reading data file and binning the data
    print >> sys.stderr, "reading data file"
    chr_Bdata = {}  
    First = True
    order = 2
    for line in open(data_fname):
        cols = line.strip().split()

        # find data type and range
        if First:
            if cols[2] == "PhysicalPosition":
                data_type = 'point'
                col_st = 3
            else:
                assert cols[2] == 'Start'
                assert cols[3] == 'End'
                data_type = 'binned'
                col_st = 4
            labels = cols[col_st:]
            First = False
            continue

        if data_type == 'point':
            chr, st, ed = cols[1], int(cols[2]), int(cols[2])
            ed +=1
        else:
            chr, st, ed = cols[1], int(cols[2]), int(cols[3])
        
        if chr not in chr_st:
            continue
        
        if math.log10(st+1) > order:
            print >> sys.stderr, chr + " st " + str(st) +" is reading"
            order += 1

        # limit the input end by chromosome size
        ed = min(ed, chr_size[chr])

        # find bins overlap with data range
        idx = st / bin_step
        bst, bed = bin_step*idx, bst + bin_size
        idx_st = idx
        while st >= bst and st < bed:
            if idx < idx_st:
                idx_st = idx
            idx -= 1
            if idx < 0:
                break
            bst, bed = bin_step*idx, bst + bin_size

        idx_ed = (ed - 1) / bin_step
        assert idx_st <= idx_ed
        
        # binning the data
        values = [float(value) for value in cols[col_st:]]
        for idx in range(idx_st, idx_ed+1):
            bst, bed = bin_step*idx, bst + bin_size
            a, b = max(st, bst), min(ed, bed)
            length = b - a

            if chr not in chr_Bdata:
                chr_Bdata[chr] = {}
            if idx not in chr_Bdata[chr]:
                chr_Bdata[chr][idx] = {}
                
            for name, value in zip(labels, values):
                if name not in chr_Bdata[chr][idx]:
                    chr_Bdata[chr][idx][name] = []

                chr_Bdata[chr][idx][name].append(value*length)
    
    # write bin data file
    print >> sys.stderr, "writing bin data file"
    
    f = open(out_fname + '_Bdata.cn', 'w')
    s = 'BinID\tChromosome\tStart\tEnd'
    s += '\t'.join(labels)
    print >> f, s

    ID = 0
    for chr in chr_list:
        try:
            chr_Bdata[chr]
        except:
            continue
        
        idx = 0        
        last_idx = chr_size[chr] / bin_step
        while idx <= last_idx:
            try:
                name_Bdata = chr_Bdata[chr][idx]
            except:
                idx +=1
                continue
            
            total = sum([sum(Bdata) for Bdata in name_Bdata.values()])
            if skip_zero and total <=0:
                idx +=1
                continue

            Bvalues = []
            for name in labels:
                Bdata = name_Bdata[name]
                if len(Bdata) <= 0:
                    Bvalue = 'NA'
                else:
                    if bin_value == 'sum':
                        Bvalue = sum(Bdata)
                    elif bin_value == 'mean':
                        Bvalue = float(sum(Bdata))/len(Bdata)
                    elif bin_value == 'median':
                        Bdata = sorted(Bdata)
                        if len(Bdata) % 2 !=0:
                            Bvalue = Bdata[len(Bdata)/2]
                        else:
                            Bvalue = 0.5*(Bdata[len(Bdata)/2] + Bdata[len(Bdata)/2 - 1])
                Bvalues.append(Bvalue)

            bst = bin_step*idx
            bed = min(bst + bin_size, chr_size[chr])
            s = '\t'.join([ID, chr, bst, bed] + Bvalues)
            print >> f, s

            ID +=1
            idx +=1
                        
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

    parser = ArgumentParser(description='Binning data cn file')
    parser.add_argument(metavar = '-f',
                        dest='data_fname',
                        type=str,
                        help='data cn file')
    parser.add_argument(metavar='-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence filename')
    parser.add_argument('--Bsize',
                        dest="bin_size",
                        type=int,
                        default=167,
                        help='bin window in bp')
    parser.add_argument('--Bstep',
                        dest="bin_step",
                        type=int,
                        default=25,
                        help='bin moving step in bp')
    parser.add_argument('--Bvalue',
                        dest="bin_value",
                        type=str,
                        default='sum',
                        help='binning value choice (sum/mean/median)')    
    parser.add_argument('--skip',
                        dest="skip_zero",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='skip the zero coverage positions')
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

    # get length for each chromosome
    chr_size = {}
    for line in open(args.ref_fname):
        line = line.strip()
        if line.startswith('>'):
            key = line[1:]
            assert key not in chr_size
            chr_size[key] = 0
            continue
        chr_size[key] += len(line)

    chr_list = []
    if not args.chr_list:
        chr_list = sorted(chr_size.keys(), cmp=chr_cmp)
    else:
        chr_list = sorted(args.chr_list, cmp=chr_cmp)

    Bin_data (args.data_fname,
              chr_size,
              args.bin_size,
              args.bin_step,
              args.bin_value,
              args.skip_zero,
              chr_list,
              args.out_fname)
