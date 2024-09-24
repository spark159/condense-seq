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

def gzopen (fname):
    if fname.endswith('.gz'):
        reading_file = gzip.open(fname, 'rb')
    else:
        reading_file = open(fname, 'r')
    return reading_file

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

# read gtab file
def read_gtab (fname,
               mode='row',
               chr_choices=None):

    if mode == 'row':
        ID_field_value = {}
    elif mode == 'col':
        field_ID_value = {}
    elif mode == 'both':
        ID_field_value, field_ID_value = {}, {}

    First = True
    data_type = None
    for line in gzopen(fname):
        line = line.strip()

        if not line:
            continue

        cols = line.split('\t')
        if First:
            if cols[1] == 'Position':
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

            field_names = cols[col_st:col_ed]
            field_idxs = range(col_st, col_ed)
                    
            First = False
            continue
        
        if data_type == 'point':
            chr, pos = cols[:col_st]
            ID = (chr, int(pos))
        elif data_type == 'binned':
            chr, start, end = cols[:col_st]
            ID = (chr, int(start), int(end))

        if chr_choices != None and chr not in chr_choices:
            continue
                    
        for field, k in zip(field_names, field_idxs):
            value = cols[k]
            try:
                value = float(value)
            except:
                pass

            if mode in ['row', 'both']:
                try:
                    ID_field_value[ID]
                except:
                    ID_field_value[ID] = {}
                ID_field_value[ID][field] = value

            if mode in ['col', 'both']:
                try:
                    field_ID_value[field]
                except:
                    field_ID_value[field] = {}
                field_ID_value[field][ID] = value
 
    if mode == 'row':
        return ID_field_value, field_names, data_type
    if mode == 'col':
        return field_ID_value, field_names, data_type
    if mode == 'both':
        return ID_field_value, field_ID_value, field_names, data_type


# read fasta and get sequence for each windows
def get_seq_from_FASTA (fasta_fname,
                        chr_ID_win,
                        mode='fullseq'):
    
    chrs = sorted(chr_ID_win.keys(), cmp=chr_cmp)
    #print chrs
    Find = False
    chr_ID_seq = {}
    for line in gzopen(fasta_fname):
        line = line.strip()
        if line.startswith(">"):
            if Find:
                while len(left) > 0:
                    leftID, win_size, seq = left.pop(0)
                    if mode == 'fullseq':
                        ID_seq[leftID] = seq
                    elif mode == 'ATcontent':
                        ID_seq[leftID] = AT_content(seq)
                assert len(ID_seq) == len(winID)
                assert chr not in chr_ID_seq
                chr_ID_seq[chr] = ID_seq
                #print chr
                #print len(ID_seq), len(winID)

            if len(chr_ID_seq) >= len(chrs):
                Find = False
                continue
                
            chr_name = line.split()[0][1:]
            if chr_name in chrs:
                chr = chr_name
                ID_win = chr_ID_win[chr]
                winID = []
                for ID, win in ID_win.items():
                    wst, wed = win
                    winID.append((wst, wed, ID))
                winID = sorted(winID)
                
                k = 0
                wst, wed, ID = winID[k]
                pt = -1
                ID_seq = {}
                left = []

                Find = True
            else:
                Find = False
            continue

        if Find:                    
            if len(left) == 0 and len(ID_seq) == len(winID):
                continue
            if len(left) == 0 and pt + len(line) < wst:
                pt += len(line)
                continue
            
            for i in range(len(left)):
                leftID, win_size, seq = left.pop(0)
                ed = min(len(line), win_size-len(seq))
                seq += line[:ed]
                if len(seq) == win_size:
                    if mode == 'fullseq':
                        ID_seq[leftID] = seq
                    elif mode == 'ATcontent':
                        ID_seq[leftID] = AT_content(seq)
                else:
                    left.append([leftID, win_size, seq])

            while pt + len(line) >= wst and k < len(winID):
                loc = wst - pt - 1
                win_size = wed - wst
                seq = line[loc:min(loc+win_size, len(line))]
                if len(seq) == win_size:
                    if mode == 'fullseq':
                        ID_seq[ID] = seq
                    elif mode == 'ATcontent':
                        ID_seq[ID] = AT_content(seq)
                else:
                    left.append([ID, win_size, seq])
                k += 1
                if k >= len(winID):
                    break
                wst, wed, ID =  winID[k]

            pt += len(line)

    if Find:
        while len(left) > 0:
            leftID, win_size, seq = left.pop(0)
            if mode == 'fullseq':
                ID_seq[leftID] = seq
            elif mode == 'ATcontent':
                ID_seq[leftID] = AT_content(seq)
        assert len(ID_seq) == len(winID)
        assert chr not in chr_ID_seq
        chr_ID_seq[chr] = ID_seq

    #print sorted(chr_ID_seq.keys(), cmp=chr_cmp)
    assert len(chrs) == len(chr_ID_seq)
    return chr_ID_seq


# read bisulfite-seq data (ENCODE bethylbed file) and count methylated C for each bin
def read_BS_file (fnames, chr_intdict):
    chr_ID_C = {} # num of "detected" C in the target motif

    for fname in fnames:
        for line in gzopen(fname):
            cols = line.strip().split()
            chrname, st, ed, _, _, strand, _, _, _, reads, frac = cols[:11]

            try:
                chr_intdict[chrname]
            except:
                continue
            
            pos = int(st)
            reads, frac = int(reads), 0.01*float(frac)
            if reads <= 0: # skip "undetected" C in the target motif
                continue

            findIDs = chr_intdict[chrname].insert(pos, frac)
            for ID in findIDs:
                if chrname not in chr_ID_C:
                    chr_ID_C[chrname] = {}
                if ID not in chr_ID_C[chrname]:
                    chr_ID_C[chrname][ID] = 0
                chr_ID_C[chrname][ID] += 1

    chr_ID_meC = {}  # num of "detected" methylated C in the target motif
    for chr in chr_intdict:
        chr_ID_meC[chr] = chr_intdict[chr].get()
        
    return chr_ID_C, chr_ID_meC


# read chip-seq data (ENCODE peak-bed file) and get signal for each bin
def read_chip_file (fnames, chr_intdict, unit='signal'):
    for fname in fnames:        
        for line in gzopen(fname):
            cols = line.strip().split()
            chrname, st, ed, peakname, _, strand, signal, pvalue, qvalue = cols[:9]

            try:
                chr_intdict[chrname]
            except:
                continue
            
            if unit == 'signal':
                score = float(signal)
            elif unit == 'pvalue':
                score = float(pvalue)
            elif unit == 'qvalue':
                score = float(qvalue)
                
            st, ed = int(st), int(ed)
            chr_intdict[chrname].insert_range(st, ed, score)

    chr_ID_value = {}
    for chr in chr_intdict:
        chr_ID_value[chr] = chr_intdict[chr].get()
        
    return chr_ID_value


# read bedgraph file and get value for each bin
def read_bedgraph_file (fnames, chr_intdict):
    for fname in fnames:
        for line in gzopen(fname):
            if not line.startswith('chr'):
                continue
            cols = line.strip().split()
            chrname, st, ed, value = cols

            try:
                chr_intdict[chrname]
            except:
                continue
            
            st, ed = int(st), int(ed)
            value = float(value)
            chr_intdict[chrname].insert_range(st, ed, value)

    chr_ID_value = {}
    for chr in chr_intdict:
        chr_ID_value[chr] = chr_intdict[chr].get()
    
    return chr_ID_value

def make_table(data_fname,
               ref_fname,
               bin_size,
               bin_step,
               bs_fnames,
               chip_fnames,
               bedgraph_fnames,
               full_seq,
               chr_list,
               genome_size,
               out_fname):    

    # read data gtab file
    field_ID_data, field_names, data_type = read_gtab(data_fname,
                                                      mode='col',
                                                      chr_choices=chr_list)

    chr_ID_win = {}
    for ID in field_ID_data[field_names[0]]:
        try:
            chr, st, ed = ID
        except:
            chr, pos = ID
            st = pos - bin_size/2
            ed = pos + bin_size/2

            if bin_size % 2 != 0:
                ed +=1
        
        if chr not in chr_ID_win:
            chr_ID_win[chr] = {}
        chr_ID_win[chr][ID] = (st, ed)

    chr_list = sorted(chr_ID_win.keys(), cmp=chr_cmp)
    print >> sys.stderr, "data gtab file reading is done"

    # extract sequence information for each bin
    if full_seq:
        chr_ID_seq = get_seq_from_FASTA(ref_fname,
                                        chr_ID_win,
                                        mode='fullseq')
    else:
        chr_ID_AT = get_seq_from_FASTA(ref_fname,
                                       chr_ID_win,
                                       mode='ATcontent')
        
    print >> sys.stderr, "reference reading is done"
    
    # build interval dictionary for each chromosome
    chr_intdict = {}
    if len(bs_fnames) + len(chip_fnames) + len(bg_fnames) > 0:
        print >> sys.stderr, "building interval dictionary"
        for chr in chr_list:
            print >> sys.stderr, "%s" % (chr)
            if bin_step:
                Int_dict = bin_hash(chr_ID_win[chr],
                                    bin_size,
                                    bin_step,
                                    genome_size[chr])
            else:
                domain_size = 10**(int(math.log10(genome_size[chr]))/2)
                Int_dict = double_hash(chr_ID_win[chr],
                                       domain_size,
                                       genome_size[chr])
            chr_intdict[chr] = Int_dict

    # read bisulfite-seq data (ENCODE bethylbed file) and count methylated C for each bin
    bs_chr_ID_C = {}
    bs_chr_ID_meC = {}
    if bs_fnames:
        for bs in bs_fnames:
            chr_intdict_copy = copy.deepcopy(chr_intdict)
            fnames = bs_fnames[bs]
            chr_ID_C, chr_ID_meC = read_BS_file(fnames, chr_intdict_copy)
            bs_chr_ID_C[bs] = chr_ID_C
            bs_chr_ID_meC[bs] = chr_ID_meC
            del chr_intdict_copy, chr_ID_C, chr_ID_meC
            print >> sys.stderr, "%s BS reading is done" % (bs)

    # read chip-seq data (ENCODE peak-bed file) and get value for each bin
    chip_chr_ID_value = {}
    if chip_fnames:
        for chip in chip_fnames:
            chr_intdict_copy = copy.deepcopy(chr_intdict)
            fnames = chip_fnames[chip]
            chr_ID_value = read_chip_file(fnames, chr_intdict_copy)
            chip_chr_ID_value[chip] = chr_ID_value
            del chr_intdict_copy, chr_ID_value
            print >> sys.stderr, "%s Chip reading is done" % (chip)

    # read bedgraph file and get value for each bin
    bg_chr_ID_value = {}
    if bg_fnames:
        for bg in bg_fnames:
            chr_intdict_copy = copy.deepcopy(chr_intdict)
            fnames = bg_fnames[bg]
            chr_ID_value = read_bedgraph_file(fnames, chr_intdict_copy)
            bg_chr_ID_value[bg] = chr_ID_value
            del chr_intdict_copy, chr_ID_value
            print >> sys.stderr, "%s Bedgraph reading is done" % (bg)

    del chr_intdict

    #start = time.time()

    # write annotation file
    print >> sys.stderr, "writing annotation table file"

    f = gzip.open(out_fname + '_table.gtab.gz', 'wb')
    if data_type == 'point':
        s = 'Chromosome\tPosition'
    else:
        s = 'Chromosome\tStart\tEnd'
        
    for field_name in field_names:
        s += '\t' + field_name

    if full_seq:
        s += '\t' + 'Sequence'
    s += '\t' + 'ATcontent'

    bs_names = sorted(bs_chr_ID_C.keys())
    for bs in bs_names:
        s += '\t' + 'CNumber(%s)' % (bs)
        s += '\t' + 'meCNumber(%s)' % (bs)

    chip_names = sorted(chip_chr_ID_value.keys())
    for chip in chip_names:
        s += '\t' + chip

    bg_names = sorted(bg_chr_ID_value.keys())
    for bg in bg_names:
        s += '\t' + bg
    print >> f, s

    for chr in chr_list:
        ID_win = chr_ID_win[chr]
        for ID in sorted(ID_win.keys()):
            try:
                chr, pos = ID
                s = chr + "\t" + str(pos)
            except:
                chr, st, ed = ID
                s = chr + "\t" + str(st) + "\t" + str(ed)            

            for field_name in field_names:
                data = field_ID_data[field_name][ID]
                try:
                    data = round(data, 5)
                except:
                    pass
                s += '\t' + str(data)
            
            if full_seq:
                seq = chr_ID_seq[chr][ID]
                s += '\t' + seq
                s += '\t' + str(round(AT_content(seq), 5))
            else:
                s += '\t' + str(round(chr_ID_AT[chr][ID], 5))

            for bs in bs_names:
                try:
                    C = int(bs_chr_ID_C[bs][chr][ID])
                except:
                    C = 0
                s += '\t' + str(C)
                try:
                    meC = round(bs_chr_ID_meC[bs][chr][ID],5)
                except:
                    meC = 0.0
                s += '\t' + str(meC)

            for chip in chip_names:
                ID_value = chip_chr_ID_value[chip][chr]
                try:
                    value = round(ID_value[ID],5)
                except:
                    value = 0.0
                s += '\t' + str(value)

            for bg in bg_names:
                ID_value = bg_chr_ID_value[bg][chr]
                try:
                    value = round(ID_value[ID],5)
                except:
                    value = 0.0
                s += '\t' + str(value)

            print >> f, s

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

    parser = ArgumentParser(description='Make annotation table by putting together gtab file with other epigenetic data')
    parser.add_argument(metavar='-f',
                        dest="data_fname",
                        type=str,
                        help='data gtab files')
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
                        dest="bs_fnames",
                        type=str,
                        nargs='+',
                        action='append',
                        help='bisulfite sequencing files (order: name, ENCODE bedMethyl files)')
    parser.add_argument('--chip',
                        dest="chip_fnames",
                        type=str,
                        nargs='+',
                        action='append',
                        help='chip sequencing files (order: name, ENCODE peak-bed files)')
    parser.add_argument('--bedgraph',
                        dest="bg_fnames",
                        type=str,
                        nargs='+',
                        action='append',
                        help='bedgraph files (order: name, bedgraph files)')
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
    for line in gzopen(args.ref_fname):
        line = line.strip()
        if line.startswith('>'):
            key = line.split()[0][1:]
            assert key not in genome_size
            genome_size[key] = 0
            continue
        genome_size[key] += len(line)

    chr_list = []
    if not args.chr_list:
        chr_list = sorted(genome_size.keys(), cmp=chr_cmp)
    else:
        chr_list = sorted(args.chr_list, cmp=chr_cmp)

    bs_fnames = {}
    if args.bs_fnames:
        for bs_fname in args.bs_fnames:
            if len(bs_fname) <= 1:
                print >> sys.stderr, "need name and file information"
                sys.exit(1)
            bs, fnames = bs_fname[0], bs_fname[1:]
            assert bs not in bs_fnames
            bs_fnames[bs] = fnames

    chip_fnames = {}
    if args.chip_fnames:
        for chip_fname in args.chip_fnames:
            if len(chip_fname) <= 1:
                print >> sys.stderr, "need name and file information"
                sys.exit(1)
            chip, fnames = chip_fname[0], chip_fname[1:]
            assert chip not in chip_fnames
            chip_fnames[chip] = fnames

    bg_fnames = {}
    if args.bg_fnames:
        for bg_fname in args.bg_fnames:
            if len(bg_fname) <= 1:
                print >> sys.stderr, "need name and file information"
                sys.exit(1)
            bg, fnames = bg_fname[0], bg_fname[1:]
            assert bg not in bg_fnames
            bg_fnames[bg] = fnames

    if args.bin_step:
        bin_step = args.bin_step
    else:
        bin_step = None
                    
    make_table(args.data_fname,
               args.ref_fname,
               args.bin_size,
               bin_step,
               bs_fnames,
               chip_fnames,
               bg_fnames,
               args.full_seq,
               chr_list,
               genome_size,
               args.out_fname)
