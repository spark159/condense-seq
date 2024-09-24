import sys
import gzip
import glob
import numpy as np
import statis_edit as statis

# open file w/o gzip compression
def gzopen (fname):
    if fname.endswith('.gz'):
        reading_file = gzip.open(fname, 'rb')
    else:
        reading_file = open(fname, 'r')
    return reading_file

# read csv file
def read_csv_file (fname,
                   mode='row',
                   delim=',',
                   header=True,
                   rowID=True,
                   jump=None):
    if rowID:
        col_st = 1
    else:
        col_st = 0
        
    ID_field_value = {}
    First = True
    counter = -1

    #for cols in csv.reader(open(fname), delimiter=delim):
    for cols in csv.reader(codecs.EncodedFile(gzopen(fname),
                                              'utf-8',
                                              'utf-8-sig'),
                           delimiter=delim):

        if First and header:
            field_names = cols[col_st:]
            First = False
            continue
        elif First and not header:
            field_names = range(len(cols[col_st:]))
            First = False
            pass

        counter += 1
        if jump and counter % jump != 0:
            continue

        if rowID:
            ID = cols[0]
        else:
            ID = counter

        if ID not in ID_field_value:
            ID_field_value[ID] = {}

        cols = cols[col_st:]
        #print cols
        for i in range(len(cols)):
            field = field_names[i]
            try:
                value = float(cols[i])
            except:
                value = cols[i]
            if field not in ID_field_value[ID]:
                ID_field_value[ID][field] = value
            else:
                if type(ID_field_value[ID][field]) != list:
                    ID_field_value[ID][field] = [ID_field_value[ID][field]]
                ID_field_value[ID][field].append(value)

    if mode == 'row':
        return ID_field_value

    if mode == 'col' or mode == 'both':
        field_ID_value = {}
        for ID in ID_field_value:
            field_value = ID_field_value[ID]
            for field in field_value:
                value = field_value[field]
                if field not in field_ID_value:
                    field_ID_value[field] = {}
                field_ID_value[field][ID] = value

    if mode == 'col':
        return field_ID_value

    if mode == 'both':
        return ID_field_value, field_ID_value

# read tabular text file
def read_tabular_file (fname,
                       mode='row',
                       delim='\t',
                       header=True,
                       rowID=True,
                       rowID_choices=None,
                       field_choices=None,
                       skip_first=0,
                       skip_nan=False):
    if rowID:
        col_st = 1
    else:
        col_st = 0

    if field_choices:
        assert header == True
        
    ID_field_value = {}
    First = True
    counter = 0

    for line in gzopen(fname):
        if skip_first > counter:
            counter +=1
            continue
        
        cols = line.strip().split(delim)
        
        if First:
            if header:            
                if field_choices == None:
                    field_names = cols[col_st:]
                    field_idxs = range(col_st, len(cols))
                else:
                    field_names, field_idxs = [], []
                    for k in range(col_st, len(cols)):
                        field_name = cols[k]
                        if field_name in field_choices:
                            field_names.append(field_name)
                            field_idxs.append(k)
            else:
                field_names = range(len(cols[col_st:]))
                field_idxs = range(col_st, len(cols))
            First = False
            counter +=1
            continue

        if rowID:
            try:
                ID = int(cols[0])
            except:
                ID = cols[0]
        else:
            ID = counter

        if rowID_choices and ID not in rowID_choices:
            counter +=1
            continue

        if ID not in ID_field_value:
            ID_field_value[ID] = {}

        for field, idx in zip(field_names, field_idxs):
            try:
                value = float(cols[idx])
            except:
                value = cols[idx]
                
            if value == 'NA':
                if skip_nan:
                    continue
                else:
                    value = np.NaN

            if field not in ID_field_value[ID]:
                ID_field_value[ID][field] = value
            else:
                if type(ID_field_value[ID][field]) != list:
                    ID_field_value[ID][field] = [ID_field_value[ID][field]]
                ID_field_value[ID][field].append(value)
        counter +=1

    if mode == 'row':
        return ID_field_value

    if mode == 'col' or mode == 'both':
        field_ID_value = {}
        for ID in ID_field_value:
            field_value = ID_field_value[ID]
            for field in field_value:
                value = field_value[field]
                if field not in field_ID_value:
                    field_ID_value[field] = {}
                field_ID_value[field][ID] = value

    if mode == 'col':
        return field_ID_value

    if mode == 'both':
        return ID_field_value, field_ID_value

# read gtab file
def read_gtab (fname,
               mode='row',
               field_choices=None,
               chr_choices=None,
               skip_nan=False,
               by_chr=False):

    # sort by chromosomes
    if by_chr:
        if mode == 'row':
            chr_ID_field_value = {}
        elif mode == 'col':
            field_chr_ID_value = {}
        elif mode == 'both':
            chr_ID_field_value, field_chr_ID_value = {}, {}
    else:
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
            else:
                assert cols[1] == 'Start'
                assert cols[2] == 'End'
                data_type = 'binned'
                col_st = 3

            if field_choices == None:
                field_names = cols[col_st:]
                field_idxs = range(col_st, len(cols))
            else:
                field_names, field_idxs = [], []
                for k in range(col_st, len(cols)):
                    field_name = cols[k]
                    if field_name not in field_choices:
                        continue
                    field_names.append(field_name)
                    field_idxs.append(k)
                    
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
            
            if value == 'NA':
                if skip_nan:
                    continue
                else:
                    value = np.NaN

            if by_chr:
                if mode in ['row', 'both']:
                    try:
                        chr_ID_field_value[chr]
                    except:
                        chr_ID_field_value[chr] = {}
                    try:
                        chr_ID_field_value[chr][ID]
                    except:
                        chr_ID_field_value[chr][ID] = {}
                    chr_ID_field_value[chr][ID][field] = value

                if mode in ['col', 'both']:
                    try:
                        field_chr_ID_value[field]
                    except:
                        field_chr_ID_value[field] = {}
                    try:
                        field_chr_ID_value[field][chr]
                    except:
                        field_chr_ID_value[field][chr] = {}
                    field_chr_ID_value[field][chr][ID] = value
            else:
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

    if by_chr:
        if mode == 'row':
            return chr_ID_field_value
        if mode == 'col':
            return field_chr_ID_value
        if mode == 'both':
            return chr_ID_field_value, field_chr_ID_value
    else:    
        if mode == 'row':
            return ID_field_value
        if mode == 'col':
            return field_ID_value
        if mode == 'both':
            return ID_field_value, field_ID_value

# read referenec fasta file and get each chromosome length
def read_genome_size(fname,
                     chr_choices=None):
    genome_size = {}
    for line in gzopen(fname):
        line = line.strip()
        if line.startswith('>'):
            chr = line.split()[0][1:]
            if chr_choices!=None and chr not in chr_choices:
                continue
            assert chr not in genome_size
            genome_size[chr] = 0
            continue
        if chr_choices!=None and chr not in chr_choices:
            continue
        genome_size[chr] += len(line)
    return genome_size

# read chromosome G-banding file
def read_Gband (fname,
                chr_choices=None):

    chr_ID_Gband = {}
    ID = 0
    for line in gzopen(fname):
        if line.startswith("#"):
            continue
        cols = line.strip().split()
        chr, start, end, name, stain = cols
        if chr_choices !=None and chr not in chr_choices:
            continue
        start, end = int(start), int(end)
        if stain.startswith('g'):
            type = stain[1:4]
            if type == 'neg':
                value = 0
            elif type == 'pos':
                value = int(stain[4:])
            else:
                assert type == 'var'
                value = np.nan
        else:
            type = stain
            value = np.nan
        if chr not in chr_ID_Gband:
            chr_ID_Gband[chr] = {}
        assert ID not in chr_ID_Gband[chr]
        chr_ID_Gband[chr][ID] = {}
        chr_ID_Gband[chr][ID]['interval'] = (start, end)
        chr_ID_Gband[chr][ID]['name'] = name
        chr_ID_Gband[chr][ID]['type'] = type
        chr_ID_Gband[chr][ID]['value'] = value
        ID +=1
    return chr_ID_Gband

# read ChromHMM file
def read_chromHMM(fname,
                  chr_choices=None,
                  state_name=None):
    
    chr_state_intervals = {}
    for line in gzopen(fname):
        cols = line.strip().split()
        chr, st, ed, state = cols[:4]
        if chr_choices and chr not in chr_choices:
            continue
        st, ed = int(st), int(ed)
        if state_name:
            state = state_name[state]
        if chr not in chr_state_intervals:
            chr_state_intervals[chr] = {}
        if state not in chr_state_intervals[chr]:
            chr_state_intervals[chr][state] = []
        chr_state_intervals[chr][state].append((st,ed))
    return chr_state_intervals

# read profile file
def read_profile(fname,
                 moving_win=None,
                 name_choice=None,
                 ID_choice=None,
                 strip_ver=True,
                 average=False):

    name_ID_profile = {}
    First = True
    for line in gzopen(fname):
        line = line.strip()
        if not line:
            continue
        
        cols = line.split() 
        if First:
            First = False
            continue

        name, ID = cols[:2]
        if strip_ver:
            ID = ID.split('.')[0]

        if name_choice and name not in name_choice:
            continue
        if ID_choice and ID not in ID_choice:
            continue

        if name not in name_ID_profile:
            name_ID_profile[name] = {}
            
        if ID not in name_ID_profile[name]:
            name_ID_profile[name][ID]= {}

        profile = []
        for value in cols[6:]:
            try:
                value = float(value)
            except:
                value = np.NaN
            profile.append(value)
        profile = np.asarray(profile)
        
        if moving_win != None:
            profile = statis.moving_average(profile, moving_win)
        name_ID_profile[name][ID] = profile

    if not average:
        return name_ID_profile

    name_mean_profile = {}
    for name in name_ID_profile:
        profiles = name_ID_profile[name].values()
        mean_profile = np.nanmean(profiles, axis=0)
        assert name not in name_mean_profile
        name_mean_profile[name] = mean_profile
    return name_mean_profile, name_ID_profile

# read GTF file
def read_GTF (fname,
              chr_list=None,
              mode="gene",
              strip_ver=True):

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

# read ENCODE RNA-seq gene quantification tsv file
def read_ENCODE_RNA_seq (fname,
                         strip_ver=True,
                         unit='FPKM'):    
    First = True
    geneID_value = {}
    for line in open(fname):
        if First:
            First = False
            continue
        
        cols = line.strip().split()

        geneID = cols[0]
        if strip_ver:
            geneID = geneID.split('.')[0]

        if unit == 'TPM':
            value = float(cols[5])
        elif unit == 'FPKM':
            value = float(cols[6])

        geneID_value[geneID] = value
    return geneID_value

# read GSEA rank file
def read_rank (fname):
    gene_value = {}
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        gene, value = cols
        value = float(value)
        gene_value[gene] = value
    return gene_value

# read GSEA output
def read_GSEA (path):
    def read_report (fname, cutoff=100):
        gs_list = []
        First = True
        for line in open(fname):
            if First:
                First = False
                continue
            cols = line.strip().split('\t')
            #name, size, nes = cols[0], cols[3], cols[5]
            name, size, nes, pvalue = cols[0], cols[3], cols[5], cols[8]
            size = int(size)
            #nes = float(nes)
            try:
                nes = float(nes)
            except:
                nes = 3
            try:
                pvalue = float(pvalue) + 10**-10 # add dummy value
            except:
                pvalue = 10**-10
            #gs_list.append({'name':name, 'size':size, 'nes':nes})
            gs_list.append({'name':name, 'size':size, 'nes':nes, 'pvalue':pvalue})
            if len(gs_list) >= cutoff:
                break
        return gs_list

    def read_gs (fname):
        gene_info = {}
        First = True
        for line in open(fname):
            if First:
                First = False
                continue
            cols = line.strip().split('\t')
            gene, rank, core, es = cols[1], cols[2], cols[-1], cols[-2]
            rank = int(rank)
            es = float(es)
            if core == 'Yes':
                core = True
            else:
                core = False
            if gene not in gene_info:
                gene_info[gene] = {}
            gene_info[gene]['rank'] = rank
            gene_info[gene]['core'] = core
            gene_info[gene]['es'] = es
        return gene_info

    pos_fname = glob.glob(path + "/gsea_report_for_na_pos_*.tsv")[0]
    pos_gs_list = []
    for gs in read_report(pos_fname):
        try:
            gs_name = gs['name']
            gene_info = read_gs(path + '/' + gs_name + '.tsv')
            gs['genes'] = gene_info
            pos_gs_list.append(gs)
        except:
            break

    neg_fname = glob.glob(path + "/gsea_report_for_na_neg_*.tsv")[0]
    neg_gs_list = []
    for gs in read_report(neg_fname):
        try:
            gs_name = gs['name']
            gene_info = read_gs(path + '/' + gs_name + '.tsv')
            gs['genes'] = gene_info
            neg_gs_list.append(gs)
        except:
            break

    return pos_gs_list, neg_gs_list

