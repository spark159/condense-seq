import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import copy
import time
#import random

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

def remap_pos (pos, st, ed, new_st, new_ed):
    new_pos = new_st + (new_ed - new_st) * float(pos - st)/(ed - st) 
    new_pos = int(round(new_pos))
    if new_pos >= new_ed:
        new_pos -= 1
    return new_pos

def is_empty (data):
    if not data:
        return True
    if type(data) == list:
        for element in data:
            if not is_empty(element):
                return False
        return True
    return False

def mean (data_list):
    if len(data_list) <= 0:
        return 'NA'
    return float(sum(data_list))/len(data_list)

# nearest-neighbor interpolation
def NN_interpolate (raw_data_list):
    data_list = copy.deepcopy(raw_data_list)
    has_data = False
    # handle first one
    for i in range(len(data_list)):
        if data_list[i] != 'NA':
            data_list[:i] = [data_list[i]]*i
            has_data = True
            break
    # no data
    if not has_data:
        return data_list
    # handle extra
    while i < len(data_list)-1:
        j = i + 1
        while j < len(data_list)-1:
            if data_list[j] != 'NA':
                break
            j += 1
        mid = int(math.ceil((i+j)/2.0))
        #print i, mid, j
        data_list[i+1:mid] = [data_list[i]]*(mid-i-1)
        if data_list[j] != 'NA':
            data_list[mid:j] = [data_list[j]]*(j-mid)
        else:
            data_list[mid:j+1] = [data_list[i]]*(j-mid+1)
        i = j
    return data_list
                
# read UCSC gene annotation table file
def read_hgtable(fname, chr_list, mode='gene'):
    ID_field_values = {}   
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            First = False
        _, geneID, chr, strand, TSS, TTS, CSS, CTS = cols[:8]
        if chr not in chr_list:
            continue
        exon_num, exon_starts, exon_ends = cols[8:11]
        exon_num = int(exon_num)
        exon_temp1 = [int(value) for value in  exon_starts.split(',')[:-1]]
        exon_temp2 = [int(value) for value in exon_ends.split(',')[:-1]]
        assert exon_num == len(exon_temp1) == len(exon_temp2)
        if strand == '+':
            TSS, TTS, CSS, CTS = int(TSS), int(TTS), int(CSS), int(CTS)
            exon_starts = sorted(exon_temp1)
            exon_ends = sorted(exon_temp2)
        else:
            TSS, TTS, CSS, CTS = int(TTS), int(TSS), int(CTS), int(CSS)
            exon_starts = sorted(exon_temp2, reverse=True)
            exon_ends = sorted(exon_temp1, reverse=True)
        exons = [[exon_starts[i], exon_ends[i]] for i in range(exon_num)]
        if geneID not in ID_field_values:
            ID_field_values[geneID] = {}
        ID_field_values[geneID]['chr'] = chr
        ID_field_values[geneID]['strand'] = strand
        ID_field_values[geneID]['TSS'] = TSS
        ID_field_values[geneID]['TTS'] = TTS
        ID_field_values[geneID]['CSS'] = CSS
        ID_field_values[geneID]['CTS'] = CTS
        ID_field_values[geneID]['exons'] = exons

    if mode in ['field', 'both']:
        field_ID_values = {}
        for ID in ID_field_values:
            field_values = ID_field_values[ID]
            for field in field_values:
                values = field_values[field]
                if field not in field_ID_values:
                    field_ID_values[field] = {}
                if ID not in field_ID_values[field]:
                    field_ID_values[field][ID] = values

    if mode == 'gene':
        return ID_field_values
    if mode == 'field':
        return field_ID_values
    if mode == 'both':
        return ID_field_values, field_ID_values

# read GTF file (updated for new GTF format)
def read_GTF (fname, chr_list=None, mode="gene", strip_ver=True):
    ID_field_values = {}
    for line in open(fname):
        if line.startswith("#"):
            continue
        cols = line.strip().split('\t')
        chr, source, feature, start, end, score, strand, frame, attribute = cols[:9]
        if not chr.startswith('chr'):
            chr = "chr" + chr
        if chr_list and chr not in chr_list:
            continue
        if feature not in ["gene", "exon", "start_codon", "stop_codon"]:
            continue
        start = int(start) - 1
        end = int(end) - 1
        attcols = attribute.strip(';').split('; ')
        tag_value = {}
        for item in attcols:
            tag, value = item.strip().split(' ')
            value = value.strip('"')
            tag_value[tag] = value
        geneID = tag_value["gene_id"]
        if strip_ver:
            geneID = geneID.split('.')[0]
        geneType = tag_value["gene_type"]
        geneName = tag_value["gene_name"]
        if geneID not in ID_field_values:
            ID_field_values[geneID] = {}
        if "chr" not in ID_field_values[geneID]:
            ID_field_values[geneID]["chr"] = chr
        if "strand" not in ID_field_values[geneID]:
            ID_field_values[geneID]["strand"] = strand
        if "geneType" not in ID_field_values[geneID]:
            ID_field_values[geneID]["geneType"] = geneType
        if "geneName" not in ID_field_values[geneID]:
            ID_field_values[geneID]["geneName"] = geneName            
        if feature == "gene":
            if strand == "+":
                TSS, TTS = start, end
            else:
                TTS, TSS = start, end
            ID_field_values[geneID]["TSS"] = TSS
            ID_field_values[geneID]["TTS"] = TTS
        if feature == "exon":
            interval = [start, end]
            if "exons" not in ID_field_values[geneID]:
                ID_field_values[geneID]["exons"] = []
            ID_field_values[geneID]["exons"].append(interval)
        if feature == "start_codon":
            if strand == "+":
                CSS = start
            else:
                CSS = end
            if "CSS" not in ID_field_values[geneID]:
                ID_field_values[geneID]["CSS"] = CSS
            else:
                prev = ID_field_values[geneID]["CSS"]
                if strand == "+":
                    CSS = min(prev, CSS)
                else:
                    CSS = max(prev, CSS)
                ID_field_values[geneID]["CSS"] = CSS
        if feature == "stop_codon":
            if strand == "+":
                CTS = end
            else:
                CTS = start
            if "CTS" not in ID_field_values[geneID]:
                ID_field_values[geneID]["CTS"] = CTS
            else:
                prev = ID_field_values[geneID]["CTS"]
                if strand == "+":
                    CTS = max(prev, CTS)
                else:
                    CTS = min(prev, CTS)
                ID_field_values[geneID]["CTS"] = CTS

    for ID in ID_field_values:
        try:
            exons = ID_field_values[ID]["exons"]
        except:
            continue
        new = []
        for start, end in sorted(exons):
            if new and new[-1][1] >= start:
                new[-1][1] = max(new[-1][1], end)
            else:
                new.append([start, end])
        ID_field_values[ID]["exons"] = new

    if mode == "gene":
        return ID_field_values
    if mode == "field" or mode == "both":
        field_ID_values = {}
        for ID in ID_field_values:
            field_values = ID_field_values[ID]
            for field in field_values:
                values = field_values[field]
                if field not in field_ID_values:
                    field_ID_values[field] = {}
                field_ID_values[field][ID] = values

    if mode == "field":
        return field_ID_values
    if mode == "both":
        return ID_field_values, field_ID_values

def Profiling (data_fname,
               region_fname,
               feature_choice,
               genome_size,
               chr_list,
               profile_len,
               domain_frac,
               up_win,
               down_win,
               interpolate,
               skip_zero,
               out_fname):

    # reading genomic region file
    print >> sys.stderr, "reading genomic region file"
    file_type = region_fname.split('.')[-1]
    if file_type == 'table':
        ID_field_values = read_hgtable(region_fname, chr_list, mode='gene')
    elif file_type == 'gtf':
        ID_field_values = read_GTF(region_fname, chr_list, mode='gene')
    else:
        None #need to work on

    chr_win = {}
    chr_start, chr_end = {}, {}
    chr_ID_intervals = {}
    chr_ID_info = {}

    chr_stpt_ID = {}
    chr_edpt_ID = {}

    feature_count = 0

    features = feature_choice.split('-')
    if len(features) == 2:
        mark1, mark2 = features
    else:
        assert len(features) == 1
        mark1, mark2 = features[0], features[0]

    for ID in ID_field_values:
        field_values = ID_field_values[ID]
        chr = field_values['chr']
        try:
            pos1, pos2 = field_values[mark1], field_values[mark2]
        except:
            continue
        strand = field_values['strand']
        if strand == '+':
            start = max(0, pos1 - up_win)
            end = min(genome_size[chr] - 1, pos2 + down_win)
            intervals = [(start, pos1), (pos1, pos2+1), (pos2+1, end+1)] 
        else:
            start = max(0, pos2 - down_win)
            end = min(genome_size[chr] - 1, pos1 + up_win)
            intervals = [(start, pos2), (pos2, pos1+1), (pos1+1, end+1)]
        win = (start, end)
        
        if chr not in chr_win:
            chr_win[chr] = []
        if chr not in chr_start:
            chr_start[chr] = []
        if chr not in chr_end:
            chr_end[chr] = []
        if chr not in chr_ID_intervals:
            chr_ID_intervals[chr] = {}
        if chr not in chr_ID_info:
            chr_ID_info[chr] = {}
        if chr not in chr_stpt_ID:
            chr_stpt_ID[chr] = []
        if chr not in chr_edpt_ID:
            chr_edpt_ID[chr] = []

        if win not in chr_win[chr]:
            chr_win[chr].append(win)
            chr_start[chr].append(start)
            chr_end[chr].append(end)
            chr_stpt_ID[chr].append([start, ID])
            chr_edpt_ID[chr].append([end, ID])
            if ID not in chr_ID_intervals[chr]:
                chr_ID_intervals[chr][ID] = intervals
            if ID not in chr_ID_info[chr]:
                chr_ID_info[chr][ID] = (ID, feature_choice, chr, str(win[0])+'-'+str(win[1]), strand)
            feature_count += 1

    for chr in chr_start:
        chr_start[chr] = sorted(chr_start[chr])
        chr_end[chr] = sorted(chr_end[chr])
        chr_stpt_ID[chr] = sorted(chr_stpt_ID[chr])
        chr_edpt_ID[chr] = sorted(chr_edpt_ID[chr])

    for chr in chr_stpt_ID:
        chr_stpt_ID[chr] = [value[1] for value in chr_stpt_ID[chr]]
        chr_edpt_ID[chr] = [value[1] for value in chr_edpt_ID[chr]]

    # make profile intervals
    if len(features) == 2:
        margin_frac = 1.0 - domain_frac
        left_len = int(round(up_win*float(profile_len*margin_frac)/(up_win+down_win)))
        right_len = int(round(down_win*float(profile_len*margin_frac)/(up_win+down_win)))
    else:
        assert len(features) == 1        
        left_len = int(round(up_win*float(profile_len)/(up_win+down_win+1)))
        right_len = profile_len - left_len - 1

    new_intervals = {}
    new_intervals['+'] = [(0, left_len), (left_len, profile_len-right_len), (profile_len-right_len, profile_len)]
    new_intervals['-'] = [(0, right_len), (right_len, profile_len-left_len), (profile_len-left_len, profile_len)]            

    # start writing profile file
    print >> sys.stderr, "writing profile file"
    f = open(out_fname + '_profile.txt', 'w')
    s = 'Sample\tID\tFeature\tChromosome\tPhysicalPosition\tStrand'
    a, b, c, d = 0, left_len, profile_len-right_len, profile_len
    for i in range(b-a):
        s += '\t' + str(i-b)
    for i in range(b, c):
        s += '\t' + str(i-b) + 'D'
    for i in range(d-c):
        s += '\t' + str(i+1)
    print >> f, s
    
    # read input data file
    print >> sys.stderr, "reading input data cn file"
    data_type = None
    col_st = None
    stpt, edpt = 0, 0
    prev_chr = None
    ID_domain = {}
    First = True
    order = None
    for line in open(data_fname):
        cols = line.strip().split()
        if First:
            # find data type and range
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

        if chr not in chr_start:
            continue

        if prev_chr != chr:
            while len(ID_domain) > 0:
                ID = ID_domain.keys().pop()
                domain = ID_domain[ID]
                del ID_domain[ID]
                ID, feature_choice, chr, position, strand = chr_ID_info[chr][ID]
                if skip_zero and is_empty(domain):
                    continue
                for i in range(len(labels)):
                    name = labels[i]
                    profile = [ mean(values) for values in domain[i] ]
                    if strand == '-':
                        profile = profile[::-1]
                    if interpolate:
                        for ist, ied in new_intervals[strand]:
                            profile[ist:ied] = NN_interpolate(profile[ist:ied])
                    s = name + "\t" + ID + "\t" + feature_choice + "\t" + chr + "\t" + position + "\t" + strand + "\t"
                    for u in range(len(profile)-1):
                        value = profile[u]
                        s += str(value) + "\t"
                    s += str(profile[-1])
                    print >> f, s
            stpt = 0 
            edpt = 0
            prev_chr = chr
            order = None

        if order == None:
            print >> sys.stderr, chr + " pos " + str(st) +" is reading"
            order = int(math.log10(max(st, 1)))
        elif int(math.log10(max(st, 1))) > order:
            print >> sys.stderr, chr + " pos " + str(st) +" is reading"
            order += 1

        while stpt < len(chr_start[chr]) and ed > chr_start[chr][stpt]:
            ID = chr_stpt_ID[chr][stpt]
            domain = []
            for i in range(len(labels)):
                domain.append([ [] for k in range(profile_len) ])
            ID_domain[ID] = domain
            stpt += 1

        while edpt < len(chr_end[chr]) and st > chr_end[chr][edpt]:
            ID = chr_edpt_ID[chr][edpt]
            domain = ID_domain[ID]
            del ID_domain[ID]
            ID, feature_choice, chr, position, strand = chr_ID_info[chr][ID]
            edpt += 1
            if skip_zero and is_empty(domain):
                continue
            for i in range(len(labels)):
                name = labels[i]
                profile = [ mean(values) for values in domain[i] ]
                if strand == '-':
                    profile = profile[::-1]
                #print profile
                if interpolate:
                    for ist, ied in new_intervals[strand]:
                        profile[ist:ied] = NN_interpolate(profile[ist:ied])
                s = name + "\t" + ID + "\t" + feature_choice + "\t" + chr + "\t" + position + "\t" + strand + "\t"
                for u in range(len(profile)-1):
                    value = profile[u]
                    s += str(value) + "\t"
                s += str(profile[-1])
                print >> f, s

        if len(ID_domain) <= 0:
            continue
        
        counts = cols[col_st:]
        for ID in ID_domain:
            intervals = chr_ID_intervals[chr][ID]
            strand = chr_ID_info[chr][ID][-1]

            if data_type == 'point':
                for k in range(len(intervals)):
                    ist, ied = intervals[k]
                    if st < ied:
                        if k > 0:
                            intervals = intervals[k:]
                            chr_ID_intervals[chr][ID] = intervals
                        break
                    
                assert st >= ist and st < ied

                new_ist, new_ied = new_intervals[strand][3 - len(intervals)]
                idx = remap_pos(st, ist, ied, new_ist, new_ied)

                for i in range(len(labels)):
                    count = float(counts[i])
                    ID_domain[ID][i][idx].append(count)
                
            else:
                wst, wed = intervals[0][0], intervals[-1][1]
                
                if st <= wst:
                    idx_st = 0
                    k = 0
                else:
                    for k in range(len(intervals)):
                        if st < intervals[k][1]:
                            break
                    assert st >= intervals[k][0] and st < intervals[k][1]    
                    new_ist1, new_ied1 = new_intervals[strand][3 - len(intervals) + k]
                    idx_st = remap_pos(st, intervals[k][0], intervals[k][1], new_ist1, new_ied1)
                    
                if ed >= wed:
                    idx_ed = profile_len
                else:
                    for v in range(k, len(intervals)):
                        if ed < intervals[v][1]:
                            break
                    assert ed >= intervals[v][0] and ed < intervals[v][1]
                    new_ist2, new_ied2 = new_intervals[strand][3 - len(intervals) + v]
                    idx_ed = remap_pos(ed, intervals[v][0], intervals[v][1], new_ist2, new_ied2)

                if k > 0:
                    chr_ID_intervals[chr][ID] = intervals[k:]

                for i in range(len(labels)):
                    count = float(counts[i])
                    for idx in range(idx_st, idx_ed):
                        ID_domain[ID][i][idx].append(count)

    while len(ID_domain) > 0:
        ID = ID_domain.keys().pop()
        domain = ID_domain[ID]
        del ID_domain[ID]
        ID, feature_choice, chr, position, strand = chr_ID_info[chr][ID]
        if skip_zero and is_empty(domain):
            continue
        for i in range(len(labels)):
            name = labels[i]
            profile = [ mean(values) for values in domain[i] ]
            if strand == '-':
                profile = profile[::-1]
            if interpolate:
                for ist, ied in new_intervals[strand]:
                    profile[ist:ied] = NN_interpolate(profile[ist:ied])
            s = name + "\t" + ID + "\t" + feature_choice + "\t" + chr + "\t" + position + "\t" + strand + "\t"
            for u in range(len(profile)-1):
                value = profile[u]
                s += str(value) + "\t"
            s += str(profile[-1])
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

    parser = ArgumentParser(description='Make data profiles projected on the gene annotation')
    parser.add_argument(metavar='--data',
                        dest="data_fname",
                        type=str,
                        help='input data cn file')
    parser.add_argument(metavar='--region',
                        dest='region_fname',
                        type=str,
                        help='Genomic region file (UCSC table/gtf/custom file)')
    parser.add_argument(metavar='-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence filename')
    parser.add_argument('--feature',
                        dest="feature_choice",
                        type=str,
                        help='TSS: tx start site \nTTS: tx terminal site\nCSS: coding start site\nCTS: conding terminal site\nESS: exon start site\nETS: exon terminal site')
    parser.add_argument('--chr',
                        dest="chr_list",
                        type=str,
                        nargs='+',
                        help='tagert chromosome list')
    parser.add_argument('--profile-len',
                        dest="profile_len",
                        type=int,
                        default=1000,
                        help='set output profile length (default: 1000 bins)')
    parser.add_argument('--domain-frac',
                        dest="domain_frac",
                        type=float,
                        default=0.6,
                        help='length fraction of profile when feature is domain (default: 0.6)')
    parser.add_argument('--up',
                        dest="up_win",
                        type=int,
                        default=1000,
                        help='up stream window size (default: 1kb)')
    parser.add_argument('--down',
                        dest="down_win",
                        type=int,
                        default=2000,
                        help='down stream window size (default: 2kb)')
    parser.add_argument('--interpolate',
                        dest="interpolate",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='smoothing the profile by NN-interpolation')
    parser.add_argument('--skip',
                        dest="skip_zero",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='skip the empty profile')
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()

    # get length for each chromosome
    genome_size = {}

    if not args.ref_fname:
        print >> sys.stderr, "Error: there is no reference file input."
        sys.exit(1)

    for line in open(args.ref_fname):
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

    profile_len = args.profile_len
    
    if not args.feature_choice:
        feature_choice = 'TSS'
    else:
        feature_choice = args.feature_choice

    Profiling (args.data_fname,
               args.region_fname,
               feature_choice,
               genome_size,
               chr_list,
               profile_len,
               args.domain_frac,
               args.up_win,
               args.down_win,
               args.interpolate,
               args.skip_zero,
               args.out_fname)
