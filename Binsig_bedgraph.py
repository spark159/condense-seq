import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import copy
import time

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

# simple hash function when binning size/step is constant
class bin_hash:
    def __init__(self,
                 bin_size,
                 bin_step,
                 max_pos):

        self.idx_value = {}
        self.idx_count = {}
        self.bin_size = bin_size
        self.bin_step = bin_step
        self.max_pos = max_pos
        #print >> sys.stderr,"hash function is built"
        
    def find(self, pos):
        find_idxs = []
        idx = pos / self.bin_step
        st = self.bin_step*idx
        ed = min(st + self.bin_size, self.max_pos+1)
        while pos >= st and pos < ed:
            find_idxs.append(idx)
            idx -= 1
            if idx < 0:
                break
            st = self.bin_step*idx
            ed = min(st + self.bin_size, self.max_pos+1)
        return find_idxs

    def insert (self, pos, value):
        find_idxs = self.find(pos)
        for idx in find_idxs:
            if idx not in self.idx_value:
                self.idx_value[idx] = 0.0
            self.idx_value[idx] += value
            if idx not in self.idx_count:
                self.idx_count[idx] = 0.0
            self.idx_count[idx] += 1
        return find_idxs
        
    def find_range(self, rst, red):
        find_idxs = []
        min_idxs = self.find(rst)
        if len(min_idxs) <= 0:
            return find_idxs
        min_idx = min(min_idxs)
        max_idx = min(red - 1, self.max_pos) / self.bin_step
        find_idxs = range(min_idx, max_idx+1)
        return find_idxs

    def insert_range (self, rst, red, value):
        find_idxs = self.find_range(rst, red)
        for idx in find_idxs:
            st = self.bin_step*idx
            ed = min(st + self.bin_size, self.max_pos+1)
            a, b = max(st, rst), min(ed, red)
            length = b - a
            if idx not in self.idx_value:
                self.idx_value[idx] = 0.0
            self.idx_value[idx] += value*length
            if idx not in self.idx_count:
                self.idx_count[idx] = 0.0
            self.idx_count[idx] += length
        return find_idxs

    def get_value (self):
        return self.idx_value

    def get_count (self):
        return self.idx_count

def Bin_sig (sig_fnames,
             genome_size,
             bin_size,
             bin_step,
             skip_zero,
             mean_choice,
             chr_list,
             sample_choice,
             out_fnames,
             out_exten):

    # build interval dictionary for each chromosome
    print >> sys.stderr, "binning the genome"
    chr_intdic = {}
    for chr in chr_list:
        Int_dict = bin_hash(bin_size, bin_step, genome_size[chr])
        chr_intdic[chr] = Int_dict

    for k in range(len(sig_fnames)):
        sig_fname = sig_fnames[k]
        prefix, in_exten = sig_fname.rsplit('.', 1)

        print >> sys.stderr, "reading %s file" % (sig_fname)
        if in_exten == 'bedgraph':
            # reading signal file and get binned signal (bedgraph format)
            target_labels = ['Binsig']
            label_chr_intdic = {}
            for label in target_labels:
                label_chr_intdic[label] = copy.deepcopy(chr_intdic)

            linecount, order = 0, -1
            for line in open(sig_fname):
                linecount +=1
                if int(math.log10(linecount)) > order:
                    print >> sys.stderr, "line %d is reading" % (linecount)
                    order +=1

                cols = line.strip().split()
                chr, st, ed, value = cols

                if chr not in chr_list:
                    continue

                st, ed = int(st), int(ed)
                value = float(value)
                label_chr_intdic['Binsig'][chr].insert_range(st, ed, value)
                

        elif in_exten == 'cn':
            # reading signal file and get binned signal (cn format)
            linecount, order = 0, -1
            First = True
            for line in open(sig_fname):
                linecount +=1
                if int(math.log10(linecount)) > order:
                    print >> sys.stderr, "line %d is reading" % (linecount)
                    order +=1

                cols = line.strip().split()

                if First:
                    labels = cols[3:]
                    if sample_choice == 'all':
                        target_idxs = range(len(labels))
                    else:
                        target_idxs = [-1]
                    target_labels = [labels[idx] for idx in target_idxs]
                    label_chr_intdic = {}
                    for label in target_labels:
                        label_chr_intdic[label] = copy.deepcopy(chr_intdic)
                    First = False
                    continue

                _, chr, pos = cols[0], cols[1], int(cols[2])

                if chr not in chr_list:
                    continue

                values = cols[3:]
                for idx in target_idxs:
                    label = target_labels[idx]
                    value = float(values[idx])
                    label_chr_intdic[label][chr].insert(pos, value)


        print >> sys.stderr, "writing bin signal file"
        out_fname = out_fnames[k]
        
        if out_exten == 'bedgraph':                    
            # write bin signal file (bedgraph format)
            f = open(out_fname + '_Bsig.bedgraph', 'w')

            for chr in sorted(chr_list):
                bin_num = int(math.ceil(float(genome_size[chr])/bin_step))
                for idx in range(bin_num):
                    st = bin_step*idx
                    ed = min(st + bin_size, genome_size[chr])
                    string = [chr, str(st), str(ed)]
                    temp = []
                    for label in target_labels:
                        try:
                            value = label_chr_intdic[label][chr].idx_value[idx]
                            count = label_chr_intdic[label][chr].idx_count[idx]
                            bin_sig = float(value) / count
                        except:
                            bin_sig = 0.0
                        temp.append(bin_sig)
                    if skip_zero and sum(temp) == 0:
                        continue
                    string += [str(bin_sig) for bin_sig in temp]
                    print >> f, '\t'.join(string)
            f.close()

        elif out_exten == 'cn':
            # write bin signal file (cn format)
            f = open(out_fname + '_Bsig.cn', 'w')
            s = 'SNP\tChromosome\tPhysicalPosition'

            for label in target_labels:
                s += '\t' + label
                print >> f, s

            ID = 0
            for chr in sorted(chr_list):
                bin_num = int(math.ceil(float(genome_size[chr])/bin_step))
                for idx in range(bin_num):
                    st = bin_step*idx
                    ed = min(st + bin_size, genome_size[chr])
                    pos = (st + ed)/2
                    string = [str(ID), chr, str(pos)]
                    temp = []
                    for label in target_labels:
                        try:
                            value = label_chr_intdic[label][chr].idx_value[idx]
                            count = label_chr_intdic[label][chr].idx_count[idx]
                            bin_sig = float(value) / count
                        except:
                            bin_sig = 0.0
                        temp.append(bin_sig)
                    if skip_zero and sum(temp) == 0:
                        continue
                    string += [str(bin_sig) for bin_sig in temp]
                    print >> f, '\t'.join(string)
                    ID +=1

            f.close()

        print >> sys.stderr, "Done"
        print >> sys.stderr, ""


if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = ArgumentParser(description='Binning the data')
    parser.add_argument(metavar = '--sig',
                        dest='sig_fnames',
                        type=str,
                        nargs='+',
                        help='signal files')
    parser.add_argument(metavar='-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence filename')
    parser.add_argument('--Bsize',
                        dest="bin_size",
                        type=int,
                        default=167,
                        help='Bin window in bp')
    parser.add_argument('--Bstep',
                        dest="bin_step",
                        type=int,
                        default=25,
                        help='Bin moving step in bp')
    parser.add_argument('--skip',
                        dest="skip_zero",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='skip the zero coverage positions')
    parser.add_argument('--mean',
                        dest="mean_choice",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='Average the binned signals for each bin')
    parser.add_argument('--chr',
                        dest="chr_list",
                        type=str,
                        nargs='+',
                        help='tagert chromosome list')
    parser.add_argument('--sample-choice',
                        dest="sample_choice",
                        type=str,
                        default='all',
                        help='sample data choice (default:all)')
    parser.add_argument('-o',
                        dest='out_fnames',
                        type=str,
                        nargs='+',
                        help='output prefix filename')
    parser.add_argument('--oe',
                        dest='out_exten',
                        default='cn',
                        type=str,
                        help='output extension')
    
    args = parser.parse_args()

    # get length for each chromosome
    genome_size = {}

    for line in open(args.ref_fname):
        line = line.strip()
        if line.startswith('>'):
            key = line[1:]
            assert key not in genome_size
            genome_size[key] = 0
            continue
        genome_size[key] += len(line)

    chr_list = []
    if not args.chr_list:
        chr_list = sorted(genome_size.keys())
    else:
        chr_list = sorted(args.chr_list)

    out_fnames = []
    if not args.out_fnames:
        for sig_fname in args.sig_fnames:
            prefix, in_exten = sig_fname.rsplit('.', 1)
            out_fnames.append(prefix)
    else:
        out_fnames = args.out_fnames

    assert len(out_fnames) == len(args.sig_fnames)

    Bin_sig (args.sig_fnames,
             genome_size,
             args.bin_size,
             args.bin_step,
             args.skip_zero,
             args.mean_choice,
             chr_list,
             args.sample_choice,
             out_fnames,
             args.out_exten)
