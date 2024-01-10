import sys
import numpy as np
import statis

def norm(L):
    total = sum(L)
    return [L[i]/float(total) for i in range(len(L))]

def read_genome_size(fname):
    genome_size = {}
    for line in open(fname):
        line = line.strip()
        if line.startswith('>'):
            key = line[1:]
            assert key not in genome_size
            genome_size[key] = 0
            continue
        genome_size[key] += len(line)
    return genome_size

def read_tabular_file (fname,
                       mode='row',
                       delim='\t',
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
    for line in open(fname):
        cols = line.strip().split(delim)
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
            if value == 'NA':
                value = np.NaN
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
            
            
def read_anot_file(fname,
                   chr_choice=None,
                   target_names=None,
                   jump=None,
                   num_max=sys.maxsize):
    ID_chr = {}
    name_ID_value = {}
    First = True
    count = 0
    counter = -1
    for line in open(fname):
        if count > num_max:
            break
        cols = line.strip().split()
        if First:
            if cols[2] == 'PhysicalPosition':
                dtype = 'point'
                ID_pos = {}
                col_st = 3
                names = cols[col_st:]
            elif cols[2] == 'Start' and cols[3] == 'End':
                dtype = 'bin'
                ID_start, ID_end = {}, {}
                col_st = 4
                names = cols[col_st:]
            First = False
            continue
        counter += 1
        if jump and counter % jump != 0:
            continue
        
        if dtype == 'point':
            ID, chr, pos = int(cols[0]), cols[1], int(cols[2])
        elif dtype == 'bin':
            ID, chr, start, end = int(cols[0]), cols[1], int(cols[2]), int(cols[3])
            ID_start[ID] = start
            ID_end[ID] = end

        if chr_choice!=None and chr not in chr_choice:
            continue
        ID_chr[ID] = chr

        if dtype == 'point':
            ID_pos[ID] = pos
        elif dtype == 'bin':
            ID_start[ID] = start
            ID_end[ID] = end
        
        cols = cols[col_st:]
        for i in range(len(cols)):
            name = names[i]
            if target_names and name not in target_names:
                continue
            if name not in name_ID_value:
                name_ID_value[name] = {}
            assert ID not in name_ID_value[name]
            try:
                value = float(cols[i])
            except:
                value = cols[i]
            if value == 'NA':
                value = np.NaN
            name_ID_value[name][ID] = value
        count += 1
    if dtype == 'point':
        return ID_chr, ID_pos, name_ID_value
    elif dtype == 'bin':
        return ID_chr, ID_start, ID_end, name_ID_value


def read_profile(fname,
                 moving_win=None,
                 name_choice=None,
                 ID_choice=None,
                 strip_ver=True,
                 average=True):

    name_ID_profile = {}
    First = True
    for line in open(fname):
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
        #profile = [float(value) for value in cols[6:]]
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

def read_tsv (fname):
    First = True
    geneID_FPKM = {}
    for line in open(fname):
        if First:
            First = False
            continue
        cols = line.strip().split()
        geneID, FPKM = cols[0], float(cols[6])
        geneID = geneID.split('.')[0]
        geneID_FPKM[geneID] = FPKM
    return geneID_FPKM


def read_hgtable(fname, chr_list=None, mode='gene'):
    ID_field_values = {}   
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            First = False
        _, geneID, chr, strand, TSS, TTS, CSS, CTS = cols[:8]
        if chr_list and chr not in chr_list:
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

"""
def read_hgtable(fname, chr_target, TSS_range=1000, TTS_range=1000, Prom_range=5000):
    feature_ID_interval = {}
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            First = False
        _, geneID, chr, strand = cols[:4]
        
        #if geneID.startswith("NR"):
            #continue
        #if geneID.startswith("XR"):
            #continue
        if chr != chr_target:
            continue
        if strand == "+":
            TSS_st = max(int(cols[4])-TSS_range/2, 0)
            TSS_ed = min(int(cols[4])+TSS_range/2 + 1, genome_size[chr])
            TTS_st = max(int(cols[5])-TTS_range/2, 0)
            TTS_ed = min(int(cols[5])+TTS_range/2 + 1, genome_size[chr])
            Body_st = min(TSS_ed+1, genome_size[chr])
            Body_ed = max(TTS_st, 0)
            Prom_st = max(TSS_st-Prom_range, 0)
            Prom_ed = max(TSS_st, 0)
        else:
            TSS_st = max(int(cols[5])-TSS_range/2, 0)
            TSS_ed = min(int(cols[5])+TSS_range/2 + 1, genome_size[chr])
            TTS_st = max(int(cols[4])-TTS_range/2, 0)
            TTS_ed = min(int(cols[4])+TTS_range/2 + 1, genome_size[chr])
            Body_st = min(TTS_ed+1, genome_size[chr])
            Body_ed = max(TSS_st, 0)
            Prom_st = min(TSS_ed+1, genome_size[chr])
            Prom_ed = min(TSS_ed+Prom_range + 1, 0)
        if TSS_ed - TSS_st > 0:
            if "TSS" not in feature_ID_interval:
                feature_ID_interval['TSS'] = {}
            if geneID not in feature_ID_interval['TSS']:
                feature_ID_interval['TSS'][geneID] = {}
            feature_ID_interval['TSS'][geneID] = [TSS_st, TSS_ed]
        if TTS_ed - TTS_st > 0:
            if "TTS" not in feature_ID_interval:
                feature_ID_interval['TTS'] = {}
            if geneID not in feature_ID_interval['TTS']:
                feature_ID_interval['TTS'][geneID] = {}
            feature_ID_interval['TTS'][geneID] = [TTS_st, TTS_ed]
        if Body_ed - Body_st > 0:
            if "Body" not in feature_ID_interval:
                feature_ID_interval['Body'] = {}
            if geneID not in feature_ID_interval['Body']:
                feature_ID_interval['Body'][geneID] = {}
            feature_ID_interval['Body'][geneID] = [Body_st, Body_ed]
        if Prom_ed - Prom_st > 0:
            if "Prom" not in feature_ID_interval:
                feature_ID_interval['Prom'] = {}
            if geneID not in feature_ID_interval['Prom']:
                feature_ID_interval['Prom'][geneID] = {}
            feature_ID_interval['Prom'][geneID] = [Prom_st, Prom_ed]
    return feature_ID_interval
"""

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


def read_GTF_old (fname, chr_list=None, mode="gene"):
    ID_field_values = {}
    for line in open(fname):
        if line.startswith("#"):
            continue
        cols = line.strip().split('\t')
        chrnum, source, feature, start, end, score, strand, frame, attribute = cols
        if not chrnum.startswith('chr'):
            chr = "chr" + chrnum
        else:
            chr = chrnum
        if chr_list and chr not in chr_list:
            continue
        if feature not in ["gene", "exon", "start_codon", "stop_codon"]:
            continue
        start = int(start) - 1
        end = int(end) - 1
        attcols = attribute.strip(' ;').split(';')
        tag_value = {}
        for item in attcols:
            tag, value = item.strip().split(' ')
            value = value.strip('"')
            tag_value[tag] = value
        geneID = tag_value["gene_id"]
        try:
            geneType = tag_value["gene_biotype"]
        except:
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

def read_RPKM (fname, gtf_fname, chr_list=None):
    gID_field_values, field_gID_values = read_GTF_old (gtf_fname, chr_list=chr_list, mode="both")
    #gID_field_values, field_gID_values = read_GTF (gtf_fname, chr_list=chr_list, mode="both")
    
    gID_exons = field_gID_values['exons']
    gID_exonlen = {}
    for gID in gID_exons:
        length = 0
        for start, end in gID_exons[gID]:
            length +=  end - start + 1
        gID_exonlen[gID] = length

    gID_exp_counts, exp_gID_counts = read_tabular_file (fname, mode="both")
    gID_counts1 = exp_gID_counts['38-Per_rep1']
    gID_counts2 = exp_gID_counts['38-Per_rep2']
    #gID_counts1 = exp_gID_counts['group4Stim_1']
    #gID_counts2 = exp_gID_counts['group4Stim_3']


    total_counts = 0.0
    gID_counts = {}
    for gID in gID_counts1:
        counts = (gID_counts1[gID] + gID_counts2[gID])*0.5
        counts += 1  # exclude the case of zero counts
        gID_counts[gID] = counts
        total_counts += counts

    gID_RPKM = {}
    for gID in gID_exonlen:
        try:
            RPM = (gID_counts[gID] / total_counts)*(10**6)
            RPKM = float(RPM)/(gID_exonlen[gID]/1000.0)
        except:
            continue
        gID_RPKM[gID] = RPKM

    return gID_RPKM

def read_RPKM_new (fname, gtf_fname, chr_list=None):
    gID_field_values, field_gID_values = read_GTF (gtf_fname, chr_list=chr_list, mode="both")
    
    gID_exons = field_gID_values['exons']
    gID_exonlen = {}
    for gID in gID_exons:
        length = 0
        for start, end in gID_exons[gID]:
            length +=  end - start + 1
        gID_exonlen[gID] = length

    gname_exp_counts, exp_gname_counts = read_tabular_file (fname, mode="both")
    gname_counts1 = exp_gname_counts['group4Stim_1']
    gname_counts2 = exp_gname_counts['group4Stim_3']

    total_counts = 0.0
    gname_counts = {}
    for gname in gname_counts1:
        counts = (gname_counts1[gname] + gname_counts2[gname])*0.5
        counts += 1  # exclude the case of zero counts
        gname_counts[gname] = counts
        total_counts += counts

    gID_RPKM = {}
    for gID in gID_exonlen:
        try:
            gname = gID_field_values[gID]['geneName']
            RPM = (gname_counts[gname] / total_counts)*(10**6)
            RPKM = float(RPM)/(gID_exonlen[gID]/1000.0)
        except:
            continue
        gID_RPKM[gID] = RPKM

    return gID_RPKM

# read score file
def read_score (fname, chr_choice=None):
    ID_pos = {}
    name_ID_score = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        if First:
            names = cols[3:]
            First = False
            continue
        ID, chr, pos = cols[:3]
        if chr_choice != None and chr not in chr_choice:
            continue
        ID = chr + ':' + ID
        pos = int(pos)
        ID_pos[ID] = pos
        scores = [float(score) for score in cols[3:]]
        for name, score in zip(names, scores):
            if name not in name_ID_score:
                name_ID_score[name] = {}
            name_ID_score[name][ID] = score
    return ID_pos, name_ID_score

# read score file
def read_bin_score (fname, chr_choice=None):
    name_chr_range_score = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        if First:
            names = cols[4:]
            First = False
            continue
        ID, chr, st, ed = cols[:4]
        if chr_choice != None and chr not in chr_choice:
            continue
        st, ed = int(st), int(ed)
        scores = [float(score) for score in cols[4:]]
        for name, score in zip(names, scores):
            if name not in name_chr_range_score:
                name_chr_range_score[name] = {}
            if chr not in name_chr_range_score[name]:
                name_chr_range_score[name][chr] = {}
            name_chr_range_score[name][chr][(st,ed)] = score
        return name_chr_range_score

# read score file
def read_bin_score_new (fname, chr_choice=None):
    ID_chr_range = {}
    name_ID_score = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        if First:
            names = cols[4:]
            First = False
            continue
        ID, chr, st, ed = cols[:4]
        if chr_choice != None and chr not in chr_choice:
            continue
        st, ed = int(st), int(ed)
        assert ID not in ID_chr_range
        ID_chr_range[ID] = [chr, (st, ed)]
        scores = [float(score) for score in cols[4:]]
        for name, score in zip(names, scores):
            if name not in name_ID_score:
                name_ID_score[name] = {}
            assert ID not in name_ID_score[name]
            name_ID_score[name][ID] = score
    return ID_chr_range, name_ID_score

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
