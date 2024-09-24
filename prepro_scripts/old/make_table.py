import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import copy
import gzip

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

# get the reverse complementary of sequence
def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

# get AT content of sequence
def AT_content(seq):
    seq = seq.upper()
    output = 0.0
    for nt in seq:
        if nt in "AT":
            output += 1.0
    return output/len(seq)

# get the number of C in the target motif in both strand
def C_motif (seq, motif='CG', both=True):
    seq = seq.upper()
    num = 0
    for i in range(len(seq)-1):
        if seq[i:i+2] == motif:
            num += 1
    if both:
        rev_seq = rev_cmp(seq)
        for i in range(len(rev_seq)-1):
            if rev_seq[i:i+2] == motif:
                num +=1
    return num

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
        self.ID_idx = {}
        for ID in ID_interval:
            st, ed = ID_interval[ID]
            assert st % self.bin_step == 0
            assert ed == st + self.bin_size
            idx = st / self.bin_step
            assert idx not in self.idx_ID
            self.idx_ID[idx] = ID
            self.ID_idx[ID] = idx
            
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

        red = min(red, self.max_pos + 1)
        max_idx = (red - 1) / self.bin_step

        for idx in range(min_idx, max_idx+1):
            try:
                ID = self.idx_ID[idx]
                find_IDs.append(ID)
            except:
                None        
        return find_IDs

    def insert_range (self, rst, red, value):
        find_IDs = self.find_range(rst, red)
        for ID in find_IDs:
            idx = self.ID_idx[ID]
            st = self.bin_step*idx
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


# read cn file
def read_cn_file (data_fname, chr_list, bin_size, col_st=3, offset=None):    
    ID_chrst, chrst_ID = {}, {}
    ID_data_list = []
    names = []
    data_type = None
    First = True

    if data_fname.endswith('.gz'):
        reading_file = gzip.open(data_fname, 'rb')
    else:
        reading_file = open(data_fname, 'r')
    
    for line in reading_file:
        cols = line.strip().split()
        if First:
            # cn file: point data
            if cols[2] == "PhysicalPosition":
                data_type = 'point'
                col_st = 3
                if offset == None:
                    offset = -bin_size/2 # point is middle of bin
            # otherwise, assume binned data
            else:
                assert cols[2] == 'Start'
                assert cols[3] == 'End'
                data_type = 'binned'
                col_st = 4
                offset = 0
            names = [name.rsplit('/', 1)[-1].rsplit('.', 1)[0] for name in cols[col_st:]]
            ID_data_list = [{} for i in range(len(names))]
            First = False
            continue
        ID, chr, pos = int(cols[0]), cols[1], int(cols[2])
        if chr not in chr_list:
            continue
        datas = [float(data) for data in cols[col_st:]]
        if chr not in chrst_ID:
            chrst_ID[chr] = {}
        st = pos + offset
        ID_chrst[ID] = (chr, st)
        chrst_ID[chr][st] = ID
        for i in range(len(ID_data_list)):
            ID_data = ID_data_list[i]
            assert ID not in ID_data
            ID_data[ID] = datas[i]

    reading_file.close()
    return ID_chrst, chrst_ID, ID_data_list, names, data_type


# read reference file and get sequence for each NCP
def get_seq(ref_fname, chr, st_ID, win):
    seq = ""
    pt = -1
    k = 0
    left = []
    Find = False
    stID = [[st, st_ID[st]] for st in sorted(st_ID.keys())]
    pos, ID = stID[k]
    ID_seq = {}

    if ref_fname.endswith('.gz'):
        reading_file = gzip.open(ref_fname, 'rb')
    else:
        reading_file = open(ref_fname, 'r')
    
    for line in reading_file:
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
                    #AT = AT_content(seq)
                    #ID_AT[leftID] = AT
                    ID_seq[leftID] = seq
                else:
                    left.append([leftID, seq])
            while pt + len(line) >= pos and k < len(stID):
                loc = pos - pt - 1
                seq = line[loc:min(loc+win,len(line))]
                if len(seq) == win:
                    #AT = AT_content(seq)
                    #ID_AT[ID] = AT
                    ID_seq[ID] = seq
                else:
                    left.append([ID, seq])
                k += 1
                try:
                    pos, ID = stID[k]
                except:
                    None
            if len(left) == 0 and len(ID_seq) == len(stID):
                break
            pt += len(line)
    while len(left) > 0:
        leftID, seq = left.pop(0)
        #AT = AT_content(seq)
        #ID_AT[leftID] = AT
        ID_seq[leftID] = seq
    assert len(ID_seq) == len(stID)

    reading_file.close()
    return ID_seq


# read bisulfite-seq data (ENCODE bethylbed file) and count methylated C for each bin
def read_BS_file (fname, Int_dict, chr):
    ID_C = {} # num of "detected" C in the target motif
    ID_meC = {}  # num of "detected" methylated C in the target motif

    if fname.endswith('.gz'):
        reading_file = gzip.open(fname, 'rb')
    else:
        reading_file = open(fname, 'r')

    for line in reading_file:
        cols = line.strip().split()
        chrname, st, ed, _, _, strand, _, _, _, reads, frac = cols[:11]
        if chrname != chr:
            continue
        pos = int(st)
        reads, frac = int(reads), 0.01*float(frac)
        if reads <= 0: # skip "undetected" C in the target motif
            continue
        findIDs = Int_dict.insert(pos, 1)
        for ID in findIDs:
            if ID not in ID_meC:
                ID_meC[ID] = 0.0
            ID_meC[ID] += frac
    ID_C = Int_dict.get()
    
    reading_file.close()
    return ID_C, ID_meC


# read chip-seq data (ENCODE peak-bed file) and get signal for each bin
def read_chip_file (fname, Int_dict, chr, unit='signal'):
    if fname.endswith('.gz'):
        reading_file = gzip.open(fname, 'rb')
    else:
        reading_file = open(fname, 'r')
    
    for line in reading_file:
        cols = line.strip().split()
        chrname, st, ed, peakname, _, strand, signal, pvalue, qvalue = cols[:9]
        if chrname != chr:
            continue
        if unit == 'signal':
            score = float(signal)
        elif unit == 'pvalue':
            score = float(pvalue)
        elif unit == 'qvalue':
            score = float(qvalue)
        st, ed = int(st), int(ed)
        Int_dict.insert_range(st, ed, score)

    reading_file.close()
    return Int_dict.get()


# read bedgraph file and get value for each bin
def read_bedgraph_file (fname, Int_dict, chr):
    if fname.endswith('.gz'):
        reading_file = gzip.open(fname, 'rb')
    else:
        reading_file = open(fname, 'r')
    
    for line in reading_file:
        if not line.startswith('chr'):
            continue
        cols = line.strip().split()
        chrname, st, ed, value = cols
        if chrname != chr:
            continue
        st, ed = int(st), int(ed)
        value = float(value)
        Int_dict.insert_range(st, ed, value)

    reading_file.close()
    return Int_dict.get()

def combine_all(data_fname,
                ref_fname,
                bin_size,
                bin_step,
                bs_fname,
                chip_fname,
                bedgraph_fname,
                full_seq,
                chr_list,
                genome_size,
                out_fname):

    # read data cn file
    ID_chrst, chrst_ID, ID_data_list, names, data_type = read_cn_file(data_fname, chr_list, bin_size)
    chr_list = chrst_ID.keys()
    print >> sys.stderr, "data cn file reading is done"

    
    # extract sequence for each bin
    ID_seq = {}
    for chr in chr_list:
        st_ID = chrst_ID[chr]
        temp = get_seq(ref_fname, chr, st_ID, bin_size)
        ID_seq.update(temp)
    print >> sys.stderr, "reference reading is done"
    
    # build interval dictionary for each chromosome
    chr_intdic = {}
    if bs_fname or chip_fname or bg_fname:
        for chr in chr_list:
            st_ID = chrst_ID[chr]
            ID_interval = {}
            for ID in st_ID.values():
                chr, st = ID_chrst[ID]
                ed = st + bin_size
                ID_interval[ID] = [st, ed]
            if bin_step:
                Int_dict = bin_hash(ID_interval, bin_size, bin_step, genome_size[chr])
                chr_intdic[chr] = Int_dict
            else:
                #Int_dict = double_hash(ID_interval, 10000, genome_size[chr])
                #Int_dict = double_hash(ID_interval, 1000, genome_size[chr])
                Int_dict = double_hash(ID_interval, 100000, genome_size[chr])
                chr_intdic[chr] = Int_dict

    # read bisulfite-seq data (ENCODE bethylbed file) and count methylated C for each bin
    bs_ID_C = {}
    bs_ID_meC = {}
    if bs_fname:
        for chr in chr_list:
            for bs in bs_fname:
                intdict = copy.deepcopy(chr_intdic[chr])
                fname = bs_fname[bs]
                ID_C, ID_meC = read_BS_file(fname, intdict, chr)
                if bs not in bs_ID_C:
                    bs_ID_C[bs] = {}
                if bs not in bs_ID_meC:
                    bs_ID_meC[bs] = {}
                bs_ID_C[bs].update(ID_C)
                bs_ID_meC[bs].update(ID_meC)
        del intdict, ID_C, ID_meC
        print >> sys.stderr, "BS reading is done"

    # read chip-seq data (ENCODE peak-bed file) and get singal for each bin
    chip_ID_signal = {}
    if chip_fname:
        for chr in chr_list:
            for chip in chip_fname:
                intdict = copy.deepcopy(chr_intdic[chr])
                fname = chip_fname[chip]
                ID_signal = read_chip_file(fname, intdict, chr)
                if chip not in chip_ID_signal:
                    chip_ID_signal[chip] = ID_signal
                chip_ID_signal[chip].update(ID_signal)
        del intdict, ID_signal
        print >> sys.stderr, "Chip reading is done"


    # read bedgraph file and get value for each bin
    bg_ID_value = {}
    if bg_fname:
        for chr in chr_list:
            for bg in bg_fname:
                intdict = copy.deepcopy(chr_intdic[chr])
                fname = bg_fname[bg]
                ID_value = read_bedgraph_file(fname, intdict, chr)
                if bg not in bg_ID_value:
                    bg_ID_value[bg] = ID_value
                bg_ID_value[bg].update(ID_value)
        del intdict, ID_value
        print >> sys.stderr, "bedgraph file reading is done"

    del chr_intdic

    #start = time.time()

    # write annotation file
    print >> sys.stderr, "writing annotation table file"

    if data_type == 'point':
        f = open(out_fname + '_table.cn', 'w')
        s = 'SNP\tChromosome\tPhysicalPosition'
    else:
        f = open(out_fname + '_table.cn', 'w')
        s = 'BinID\tChromosome\tStart\tEnd'
        
    for i in range(len(ID_data_list)):
        s += '\t' + names[i]
    if full_seq:
        s += '\t' + 'Sequence'
    s += '\t' + 'ATcontent'
    bs_names = sorted(bs_ID_C.keys())
    for bs in bs_names:
        s += '\t' + 'CNumber(%s)' % (bs)
        s += '\t' + 'meCNumber(%s)' % (bs)
    chip_names = sorted(chip_ID_signal.keys())
    for chip in chip_names:
        s += '\t' + chip
    bg_names = sorted(bg_ID_value.keys())
    for bg in bg_names:
        s += '\t' + bg
    print >> f, s


    IDs = sorted(ID_data_list[0].keys())
    for ID in IDs:
        chr, st = ID_chrst[ID]
        
        if data_type == 'point':
            pos = st + bin_size/2
            s = str(ID) + "\t" + chr + "\t" + str(pos)
        else:
            ed = st + bin_size
            s = str(ID) + "\t" + chr + "\t" + str(st) + "\t" + str(ed)

        for i in range(len(ID_data_list)):
            ID_data = ID_data_list[i]
            try:
                data = round(ID_data[ID], 5)
            except:
                data = ID_data[ID]
            s += '\t' + str(data)
            
        seq = ID_seq[ID]
        if full_seq:
            s += '\t' + seq
        s += '\t' + str(round(AT_content(seq), 5))

        for bs in bs_names:
            try:
                C = int(bs_ID_C[bs][ID])
            except:
                C = 0
            s += '\t' + str(C)
            try:
                meC = round(bs_ID_meC[bs][ID],5)
            except:
                meC = 0.0
            s += '\t' + str(meC)

        for chip in chip_names:
            ID_signal = chip_ID_signal[chip]
            try:
                signal = round(ID_signal[ID],5)
            except:
                signal = 0.0
            s += '\t' + str(signal)

        for bg in bg_names:
            ID_value = bg_ID_value[bg]
            try:
                value = round(ID_value[ID],5)
            except:
                value = 0.0
            s += '\t' + str(value)

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

    parser = ArgumentParser(description='Make annotation table by putting together cn file with other epigenetic data')
    parser.add_argument(metavar='-f',
                        dest="data_fname",
                        type=str,
                        help='data cn files')
    parser.add_argument(metavar = '-x',
                        dest='ref_fname',
                        type=str,
                        help='reference genome file')
    parser.add_argument('--binsize',
                        dest="bin_size",
                        type=int,
                        default=171,
                        help='bin window size of cn file in bp (default: 171)')
    parser.add_argument('--binstep',
                        dest="bin_step",
                        type=int,
                        default=0,
                        help='bin step size for cn file in bp (regular binning case)')
    parser.add_argument('--bs',
                        dest="bs_fname",
                        type=str,
                        nargs='+',
                        help='bisulfite sequencing files (order: name, ENCODE bedMethyl file)')
    parser.add_argument('--chip',
                        dest="chip_fname",
                        type=str,
                        nargs='+',
                        help='chip sequencing files (order: name, ENCODE peak-bed file)')
    parser.add_argument('--bedgraph',
                        dest="bg_fname",
                        type=str,
                        nargs='+',
                        help='bedgraph files (order: name, bedgraph file)')
    parser.add_argument('--full-seq',
                        dest="full_seq",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='write full sequence information')
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
        chr_list = sorted(genome_size.keys(), cmp=chr_cmp)
    else:
        chr_list = sorted(args.chr_list, cmp=chr_cmp)

    bs_fname = {}
    if args.bs_fname:
        flen = len(args.bs_fname)
        if flen % 2 != 0:
            print >> sys.stderr, "not right order of names or miss name"
            sys.exit(1)
        for i in range(flen):
            if i % 2 == 0:
                bs = args.bs_fname[i]
            else:
                fname = args.bs_fname[i]
                bs_fname[bs] = fname

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

    bg_fname = {}
    if args.bg_fname:
        flen = len(args.bg_fname)
        if flen % 2 != 0:
            print >> sys.stderr, "not right order of names or miss name"
            sys.exit(1)
        for i in range(flen):
            if i % 2 == 0:
                bg = args.bg_fname[i]
            else:
                fname = args.bg_fname[i]
                bg_fname[bg] = fname

    if args.bin_step:
        bin_step = args.bin_step
    else:
        bin_step = None
                    
    combine_all(args.data_fname,
                args.ref_fname,
                args.bin_size,
                bin_step,
                bs_fname,
                chip_fname,
                bg_fname,
                args.full_seq,
                chr_list,
                genome_size,
                args.out_fname)
