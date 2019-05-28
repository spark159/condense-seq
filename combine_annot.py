import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import copy
import time

def binary_search (sortlist, target):
    st, ed = 0, len(sortlist)-1
    while st <= ed:
        mid = (st+ed) / 2
        if sortlist[mid] == target:
            return mid
        elif sortlist[mid] > target:
            ed = mid - 1 
        elif sortlist[mid] < target:
            st = mid + 1
    return st

def tuple_cmp (a,b):
    if a[0] < b[0]:
        return -1
    elif a[0] > b[0]:
        return 1
    else:
        if a[1] < b[1]:
            return -1
        elif a[1] > b[1]:
            return 1
        else:
            return 0

# simple hash function when binning size/step is constant
class bin_hash:

    def __init__(self,
                 ID_interval,
                 bin_size,
                 bin_step,
                 max_pos):

        self.ID_value = {}    
        self.bin_size = bin_size
        self.bin_step = bin_step
        self.max_pos = max_pos

        # map bin idx to bin ID
        self.idx_ID = {} 
        for ID in ID_interval:
            st, ed = ID_interval[ID]
            assert st % self.bin_step == 0
            assert ed == st + self.bin_size
            idx = st / self.bin_step
            assert idx not in self.idx_ID
            self.idx_ID[idx] = ID
        
        print >> sys.stderr, "hash function is built"
        
    def find(self, pos):
        find_IDs = []
        idx = pos / self.bin_step
        st = self.bin_step*idx
        ed = st + self.bin_size
        while pos >= st and pos < ed:
            try:
                ID = self.idx_ID[idx]
                find_IDs.append(ID)
            except:
                None
            idx -= 1
            if idx < 0:
                break
            st = self.bin_step*idx
            ed = st + self.bin_size
        return find_IDs

    def insert (self, pos, value):
        find_IDs = self.find(pos)
        for ID in find_IDs:
            if ID not in self.ID_value:
                self.ID_value[ID] = 0.0
            self.ID_value[ID] += value
        return find_IDs
        
    def find_range(self, rst, red):
        find_IDs = []

        idx = rst / self.bin_step
        min_idx = idx
        st = self.bin_step*idx
        ed = st + self.bin_size
        while rst >= st and rst < ed:
            if idx < min_idx:
                min_idx = idx
            idx -= 1
            if idx < 0:
                break
            st = self.bin_step*idx
            ed = st + self.bin_size
        max_idx = min((red - 1) / self.bin_step, self.max_pos / self.bin_step)

        for idx in range(min_idx, max_idx):
            try:
                ID = self.idx_ID[idx]
                find_IDs.append(ID)
            except:
                None        
        return find_IDs

    def insert_range (self, rst, red, value):
        find_IDs = self.find_range(rst, red)
        for ID in find_IDs:
            st = self.bin_step*ID
            ed = st + self.bin_size
            a, b = max(st, rst), min(ed, red)
            length = b - a
            if ID not in self.ID_value:
                self.ID_value[ID] = 0.0
            self.ID_value[ID] += value*length
        return find_IDs

    def keys (self):
        return self.ID_value.keys()

    def values (self):
        return self.ID_value.values()

    def ID (self, id):
        return self.ID_value[id]
        
    def get (self):
        return self.ID_value

                
# build interval dictionary by using double hashing
class double_hash:
    def __init__(self,
                 ID_interval,
                 domain_size,
                 max_pos):

        self.ID_value = {}
        self.ID_interval = ID_interval

        edID = []
        for ID, interval in ID_interval.items():
            st, ed = interval
            edID.append([ed, ID])
        edID = sorted(edID, cmp=tuple_cmp)

        edlist, IDlist = [], []
        for ed, ID in edID:
            edlist.append(ed)
            IDlist.append(ID)
        
        self.domain_size = domain_size
        self.max_pos = max_pos
        self.domain_IDs = {}
        self.domain_num = max_pos // domain_size + 1
        
        for i in xrange(self.domain_num):
            self.domain_IDs[i] = []
            dst = i * self.domain_size
            ded = min(dst + self.domain_size, max_pos+1)
            idx1 = binary_search(edlist, dst)
            if idx1 == len(edlist):
                continue
            for j in xrange(idx1, len(edlist)):
                ID = IDlist[j]
                st, ed = self.ID_interval[ID]
                if st < ded:
                    self.domain_IDs[i].append(ID)
                
        print >> sys.stderr, "hash function is built"

    def __str__ (self):
        print "%s\t%s\t%s\t%s" % ("ID", "st", "ed", "value")
        for ID, value in self.ID_value.items():
            st, ed = ID_interval[ID]
            print "%d\t%d\t%d\t%f" % (ID, st, ed, value)
        return

    def find (self, pos):
        find_IDs = []
        
        domain = pos // self.domain_size
        IDs = self.domain_IDs[domain]

        edID = []
        for ID in IDs:
            st, ed  = self.ID_interval[ID]
            edID.append([ed,ID])
        edID = sorted(edID, cmp=tuple_cmp)

        edlist, IDlist = [], []
        for ed, ID in edID:
            edlist.append(ed)
            IDlist.append(ID)
            
        idx = binary_search(edlist, pos)
        if idx == len(edlist):
            return find_IDs

        for i in xrange(idx, len(edlist)):
            ID = IDlist[i]
            st, ed = self.ID_interval[ID]
            if pos >= st and pos < ed:
                find_IDs.append(ID)
            
        return find_IDs

    def insert (self, pos, value):
        find_IDs = self.find(pos)
        for ID in find_IDs:
            if ID not in self.ID_value:
                self.ID_value[ID] = 0.0
            self.ID_value[ID] += value
        return find_IDs

    def find_range (self, rst, red):
        find_IDs = []
        domain1 = rst // self.domain_size
        domain2 = red // self.domain_size

        IDs = set([])
        for i in xrange(domain1, domain2 + 1):
            IDs |= set(self.domain_IDs[i])
        IDs = list(IDs)

        edID = []
        for ID in IDs:
            st, ed = self.ID_interval[ID]
            edID.append([ed, ID])
        edID = sorted(edID, cmp=tuple_cmp)

        edlist, IDlist = [], []
        for ed, ID in edID:
            edlist.append(ed)
            IDlist.append(ID)

        idx1 = binary_search(edlist, rst)
        if idx1 == len(edlist):
            return find_IDs
        
        for j in xrange(idx1, len(edlist)):
            ID = IDlist[j]
            st, ed = self.ID_interval[ID]
            if st < red:
                find_IDs.append(ID)

        return find_IDs
    
    def insert_range (self, rst, red, value):
        find_IDs = self.find_range(rst, red)
        for ID in find_IDs:
            st, ed = self.ID_interval[ID]
            a, b = max(st, rst), min(ed, red)
            length = b - a
            if ID not in self.ID_value:
                self.ID_value[ID] = 0.0
            self.ID_value[ID] += value*length
        return find_IDs

    def keys (self):
        return self.ID_value.keys()

    def values (self):
        return self.ID_value.values()

    def ID (self, id):
        return self.ID_value[id]
        
    def get (self):
        return self.ID_value


def combine_all(Ncov_fname,
                ref_fname,
                NCP_len,
                bin_step,
                bs_fname,
                chip_fname,
                chr_list,
                genome_size,
                out_fname):

    # metric for condensabiltiy
    def metric (test, control):
        if test <=0:
            test = 1.0
        if control <= 0:
            control = 1.0
        rcount = float(test)/float(control)
        return -math.log(rcount)

    # read Ncov file and get NCP start position and condensability
    def read_Ncov_file (Ncov_fname, chr_list, NCP_len):    
        ID_chrst, chrst_ID = {}, {}
        ID_metric_list = []
        names = []
        First = True
        for line in open(Ncov_fname):
            cols = line.strip().split()
            if First:
                names = cols[3:]
                ID_metric_list = [{} for i in range(len(names)-1)]
                First = False
                continue
            ID, chr, dyad = int(cols[0]), cols[1], int(cols[2])
            if chr not in chr_list:
                continue
            counts = [float(count) for count in cols[3:]]
            if sum(counts) <= 0:
                continue
            if chr not in chrst_ID:
                chrst_ID[chr] = {}
            st = dyad - NCP_len/2
            ID_chrst[ID] = (chr, st)
            chrst_ID[chr][st] = ID
            for i in range(len(ID_metric_list)):
                ID_metric = ID_metric_list[i]
                assert ID not in ID_metric
                ID_metric[ID] = metric(counts[i], counts[-1])
        return ID_chrst, chrst_ID, ID_metric_list, names
    
    ID_chrst, chrst_ID, ID_metric_list, names = read_Ncov_file(Ncov_fname, chr_list, NCP_len)
    chr_list = chrst_ID.keys()
    print >> sys.stderr, "Ncov reading is done"
    
    # get AT content of sequence
    def AT_content(seq):
        seq = seq.upper()
        output = 0.0
        for nt in seq:
            if nt in "AT":
                output += 1.0
        return output/len(seq)

    # get CpG motif number
    def CpG_num (seq):
        seq = seq.upper()
        num = 0
        for i in range(len(seq)-1):
            if seq[i:i+2] == 'CG':
                num += 1
        num = 2*num # both strands
        return num
        
    # read reference file and get AT content for each NCP
    def get_AT(ref_fname, chr, st_ID, win):
        seq = ""
        pt = -1
        k = 0
        left = []
        Find = False
        stID = [[st, st_ID[st]] for st in sorted(st_ID.keys())]
        pos, ID = stID[k]
        ID_AT = {}
        for line in open(ref_fname):
            line = line.strip()
            if line.startswith(">"):
                if Find:
                    break
                if line[1:] == chr:
                    Find = True
                continue
            if Find:
                if len(left) == 0 and pt + len(line) < pos:
                    pt += len(line)
                    continue
                for i in range(len(left)):
                    leftID, seq = left.pop(0)
                    ed = min(len(line), win-len(seq))
                    seq += line[:ed]
                    if len(seq) == win:
                        AT = AT_content(seq)
                        ID_AT[leftID] = AT
                    else:
                        left.append([leftID, seq])
                while pt + len(line) >= pos and k < len(stID):
                    loc = pos - pt - 1
                    seq = line[loc:min(loc+win,len(line))]
                    if len(seq) == win:
                        AT = AT_content(seq)
                        ID_AT[ID] = AT
                    else:
                        left.append([ID, seq])
                    k += 1
                    try:
                        pos, ID = stID[k]
                    except:
                        None
                if len(left) == 0 and len(ID_AT) == len(stID):
                    break
                pt += len(line)
        while len(left) > 0:
            leftID, seq = left.pop(0)
            AT = AT_content(seq)
            ID_AT[leftID] = AT
        assert len(ID_AT) == len(stID)
        return ID_AT
    
    ID_AT = {}
    for chr in chr_list:
        st_ID = chrst_ID[chr]
        temp = get_AT(ref_fname, chr, st_ID, NCP_len)
        ID_AT.update(temp)
    print >> sys.stderr, "reference reading is done"

    # build interval dictionary for each chromosome
    chr_intdic = {}
    if bs_fname or chip_fname:
        for chr in chr_list:
            st_ID = chrst_ID[chr]
            ID_interval = {}
            for ID in st_ID.values():
                chr, st = ID_chrst[ID]
                ed = st + NCP_len
                ID_interval[ID] = [st, ed]
            if bin_step:
                Int_dict = bin_hash(ID_interval, NCP_len, bin_step, genome_size[chr])
                chr_intdic[chr] = Int_dict
            else:
                #Int_dict = double_hash(ID_interval, 10000, genome_size[chr])
                #Int_dict = double_hash(ID_interval, 1000, genome_size[chr])
                Int_dict = double_hash(ID_interval, 100000, genome_size[chr])
                chr_intdic[chr] = Int_dict

    # read bisulfite-seq file and get the number of methylated GC for each NCP
    def read_BS_file (fname, Int_dict, chr):
        ID_CpG = {} # num of CpG site
        ID_me = {}  # num of methylated CpG site
        First = True
        line_count = -1
        for line in open(fname):
            line_count +=1
            cols = line.strip().split()
            if First:
                First = False
                continue
            chrname, pos, strand, count, total = cols
            if chrname != chr:
                continue
            pos = int(pos) - 1
            findIDs = Int_dict.insert(pos, 1)
            count, total = int(count), int(total)
            if total <= 0:
                sig = 0.5
                #continue
            else:
                sig = float(count) / total
            for ID in findIDs:
                if ID not in ID_me:
                    ID_me[ID] = 0.0
                ID_me[ID] += sig
        ID_CpG = Int_dict.get()
        return ID_CpG, ID_me

    ID_CpG = {}
    ID_me = {}
    if bs_fname:
        for chr in chr_list:
            intdict = copy.deepcopy(chr_intdic[chr])
            temp1, temp2 = read_BS_file(bs_fname, intdict, chr)
            ID_CpG.update(temp1)
            ID_me.update(temp2)
        del intdict
        print >> sys.stderr, "BS reading is done"

    # read chip-seq file and get signal for each NCP 
    def read_chip_file (fname, Int_dict, chr):
        count = -1
        for line in open(fname):
            count += 1
            cols = line.strip().split()
            chrname, st, ed = cols[0], int(cols[1]), int(cols[2])
            if chrname != chr:
                continue
            score = float(cols[4])
            Int_dict.insert_range(st, ed+1, score)
        return Int_dict.get()

    chip_ID_signal = {}
    if chip_fname:
        for chr in chr_list:
            for chip in chip_fname:
                intdict = copy.deepcopy(chr_intdic[chr])
                fname = chip_fname[chip]
                ID_signal = read_chip_file(fname, intdict, chr)
                if chip not in chip_ID_signal:
                    chip_ID_signal[chip] = ID_signal
                else:
                    chip_ID_signal[chip].update(ID_signal)
        del intdict
        print >> sys.stderr, "Chip reading is done"

    del chr_intdic

    #start = time.time()
    
    # write annotation file
    print >> sys.stderr, "writing annotation file"
    
    f = open(out_fname + '_anot.cn', 'w')
    s = 'SNP\tChromosome\tPhysicalPosition'
    for i in range(len(ID_metric_list)):
        s += '\t' + names[i]
    s += '\t' + 'ATcontent'
    if bs_fname:
        s += '\t' + 'CpGNumber'
        s += '\t' + 'meGCNumber'
    chip_names = sorted(chip_ID_signal.keys())
    for chip in chip_names:
        s += '\t' + chip
    print >> f, s

    IDs = sorted(ID_metric_list[0].keys())
    for ID in IDs:
        chr, st = ID_chrst[ID]
        pos = st + NCP_len/2
        s = str(ID) + "\t" + chr + "\t" + str(pos)
        for i in range(len(ID_metric_list)):
            ID_metric = ID_metric_list[i]
            metric = round(ID_metric[ID], 5)
            s += '\t' + str(metric)
        s += '\t' + str(round(ID_AT[ID], 5))
        if bs_fname:
            try:
                CpG = int(ID_CpG[ID])
            except:
                CpG = 0
            s += '\t' + str(CpG)
            try:
                me = round(ID_me[ID],5)
            except:
                me = 0.0
            s += '\t' + str(me)
        for chip in chip_names:
            ID_signal = chip_ID_signal[chip]
            try:
                signal = round(ID_signal[ID],5)
            except:
                signal = 0.0
            s += '\t' + str(signal)
        print >> f, s

    f.close()

    print >> sys.stderr, "Done"

    #end = time.time()
    #print end-start
    
if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = ArgumentParser(description='combine all annotations of mono-nucleosome data')
    parser.add_argument(metavar='--Ncov',
                        dest="Ncov_fname",
                        type=str,
                        help='Ncov cn file')
    parser.add_argument(metavar = '--ref',
                        dest='ref_fname',
                        type=str,
                        help='reference genome file')
    parser.add_argument('--Nlen',
                        dest="NCP_len",
                        type=int,
                        default=171,
                        help='Mono-nucleosomal length in bp')
    parser.add_argument('--Bin',
                        dest="Bin",
                        type=int,
                        nargs='+',
                        help='bin size and step in bp')
    parser.add_argument('--bs',
                        dest="bs_fname",
                        type=str,
                        help='bisulfite sequencing file')
    parser.add_argument('--chip',
                        dest="chip_fname",
                        type=str,
                        nargs='+',
                        help='chip sequencing files (order: name, file)')
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

    chip_fname = {}
    if args.chip_fname:
        flen = len(args.chip_fname)
        if flen % 2 != 0:
            print >> sys.stderr, "not right order of names or miss name"
            sys.exit(1)
        for i in range(flen):
            if i % 2 == 0:
                chip = args.chip_fname[i]
            else:
                fname = args.chip_fname[i]
                chip_fname[chip] = fname

    if args.Bin:
        if len(args.Bin) != 2:
            print >> sys.stderr, "provide bin 'size' and 'step' information"
            sys.exit(1)
        bin_size, bin_step = args.Bin
        NCP_len = bin_size
    else:
        bin_step = None
        NCP_len = args.NCP_len
                    
    combine_all(args.Ncov_fname,
                args.ref_fname,
                NCP_len,
                bin_step,
                args.bs_fname,
                chip_fname,
                chr_list,
                genome_size,
                args.out_fname
               )
