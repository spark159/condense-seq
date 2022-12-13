import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
from sklearn import linear_model
import random
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from scipy.optimize import curve_fit
from sklearn import linear_model

# temporal sigmoid function
def sigmoid(x, L ,x0, k):
    y = L / (1 + np.exp(k*(x-x0)))
    return (y)


# temporal function for reading RNA-seq data 
def read_txtRNA_seq (fname):
    gname_exonlen = {}
    gname_count = {}
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        line = line.strip()
        if not line:
            continue
        if line.startswith('Geneid'):
            continue
        cols = line.split()
        gname, exonlen, count = cols[0], cols[-2], cols[-1]
        exonlen, count = int(exonlen), int(count)
        assert gname not in gname_exonlen
        assert gname not in gname_count
        gname_exonlen[gname] = exonlen
        gname_count[gname] = count

    totalcount = sum(gname_count.values())
    gname_RPKM = {}
    for gname in gname_count:
        gname_RPKM[gname] = (10**9)*float(gname_count[gname])/(totalcount*gname_exonlen[gname])
    return gname_RPKM

def read_Gband (fname, chr_choices=None):
    chr_ID_Gband = {}
    ID = 0
    for line in open(fname):
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

def read_bin (fname, tlen=False, chr_choices=None):
    First = True
    ID_pos = {}
    name_ID_count = {}
    for line in open(fname):
        cols = line.strip().split()
        if First:
            if tlen:
                names = cols[4:-2]
            else:
                names = cols[4:-1]
            First = False
            continue
        ID, chr, start, end = cols[:4]
        start, end = int(start), int(end)
        if chr_choices!=None and chr not in chr_choices:
            continue
        pos = int(float(start + end)/2)
        assert ID not in ID_pos
        ID_pos[ID] = pos
        for i in range(len(names)):
            name = names[i]
            count = float(cols[4+i]) + 1
            if name not in name_ID_count:
                name_ID_count[name] = {}
            assert ID not in name_ID_count[name]
            name_ID_count[name][ID] = count

    name_total = {name:float(sum(name_ID_count[name].values())) for name in name_ID_count}
    name_ID_score = {}
    for name in names[:-1]:
        for ID in name_ID_count[name]:
            control_count = name_ID_count[names[-1]][ID]
            if control_count <= 1:
                continue
            control_rcount = control_count / name_total[names[-1]]
            rcount = name_ID_count[name][ID] / name_total[name]
            score = -np.log(rcount/control_rcount)
            if name not in name_ID_score:
                name_ID_score[name] = {}
            assert ID not in name_ID_score[name]
            name_ID_score[name][ID] = score
        
    return ID_pos, name_ID_score

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

def read_bin_num (fname, chr_choices=None):
    First = True
    ID_pos = {}
    name_ID_score = {}
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        if First:
            names = []
            for name in cols[4:-1]:
                agent = name.rsplit('.', 1)[0].split('-')[-2]
                tnum = int(name.rsplit('.', 1)[0].split('-')[-1])
                names.append(tnum) 
            First = False
            continue
        if cols[-1] == '*':
            continue
        ID, chr, start, end = cols[:4]
        start, end = int(start), int(end)
        if chr_choices!=None and chr not in chr_choices:
            continue
        pos = int(float(start + end)/2)
        assert ID not in ID_pos
        ID_pos[ID] = pos
        counts = [int(count) for count in cols[4:]]
        if sum(counts) <= 0:
            continue
        control = float(counts[-1]) + 1
        if control <=1:
            continue
        for i in range(len(counts)-1):
            count = counts[i] + 1
            score = float(count)/control
            #score = -np.log(float(count)/control)
            name = names[i]
            if name not in name_ID_score:
                name_ID_score[name] = {}
            name_ID_score[name][ID] = score
        
    return ID_pos, name_ID_score


def read_bin_Chalf (fname, chr_choices=None):
    ID_pos = {}
    #name_ID_value = {'Chalf':{}}
    ID_value = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        if First:
            First = False
            continue
        cols = line.split('\t')
        binID, chr, st, ed, Chalf = cols
        if chr_choices and chr not in chr_choices:
            continue
        st, ed = int(st), int(ed)
        ID = (st, ed)
        st, ed = int(st), int(ed)
        Chalf = float(Chalf)
        pos = int(float(st + ed)/2)
        ID_pos[ID] = pos
        #name_ID_value['Chalf'][ID] = Chalf
        ID_value[ID] = Chalf
    #return ID_pos, name_ID_value
    return ID_pos, ID_value
        
def read_ATACbin (fname, chr_choices=None):
    First = True
    ID_pos, ID_score = {}, {}
    for line in open(fname):
        if First:
            First = False
            continue
        cols = line.strip().split()
        ID, chr, pos, value = cols
        if chr_choices and chr not in chr_choices:
            continue
        ID = int(ID)
        pos = int(pos)
        value = float(value)
        score = np.log2(1+value)
        assert ID not in ID_pos
        assert ID not in ID_score
        ID_pos[ID] = pos
        ID_score[ID] = score
    return ID_pos, ID_score
        

# read A/B compartment annotation
def read_eigenfile (fname, bin_size=1000000):
    eigenBinID_value = []
    eigenBinID_interval = []
    i = 0
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        try:
            value = float(line)
        except:
            value = np.nan
        st = i*bin_size
        ed = (i+1)*bin_size
        eigenBinID_value.append(value)
        eigenBinID_interval.append((st,ed))
        i +=1
    return eigenBinID_value, eigenBinID_interval

# read A/B compartment annotation (bedgraph)
def read_eigenbedgraph (fname, chr_choice):
    eigenBinID_value = []
    eigenBinID_interval = []
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        try:
            cols = line.split()
            assert len(cols) == 4
        except:
            continue
        chr, st, ed, value = cols
        if chr != chr_choice:
            continue
        st, ed = int(st), int(ed)
        value = float(value)
        eigenBinID_value.append(value)
        eigenBinID_interval.append((st, ed))
    return eigenBinID_value, eigenBinID_interval


#### temporal
tnum_conc, tnum_frac = read_titration("mCD8T_WT-NCP_sp_titration.csv")
#tnum_conc, tnum_frac = read_titration("mCD8T_inht-NCP_sp_titration.csv")
#tnum_conc, tnum_frac = read_titration("mCD8T_KO-NCP_sp_titration.csv")


#### set parameters
# set binning resolution

i = 20
#i = 10
#i = 5
bin_size = int(0.5*(10**6) / i) # binsize (unit of bp)
#bin_size = int(0.5*(10**5) / i)
#blur_win = 1
#blur_win = int(2*i + 1)
blur_win = int(4*i + 1) # sliding window (unit of bin)
#blur_win = int(6*i + 1)
#blur_win = int(10*i + 1)

# set chromosomes
chr_choices = ['chr1']
#chr_choices = ['chr10']
#chr_choices = ['chr%d' % (i) for i in range(1, 13)]
#chr_choices = ['chr1']
#chr_choices = ['chrX', 'chrY']
#chr_choices = ['chr%d' % (i) for i in range(1, 23)] + ['chrX', 'chrY']
#chr_choices = ['chr%d' % (i) for i in range(1, 4)]

# set target names and feature names
target_names = [k for k in range(1,10)]
#target_names = ['WT-Chalf', 'inht-Chalf', 'KO-Chalf']
#target_names = ['WT-Chalf']
#target_names = ["H1-NCP-HP1a-%s.bam" % (k) for k in range(1,6)]
#target_names = ["H1-DNA-HP1a-%d.bam" % (k) for k in range(1,6)]
#target_names = ["GM-NCP-sp-%d.bam" % (k) for k in range(1,10)]
#target_names = ["GM-NCP-sp-4.bam", "GM-NCP-sp-8.bam"]
#target_names = ["mCD8T-WT-NCP-sp-%d.bam" % (k) for k in range(1,10)]
#target_names = ["mCD8T-inht-NCP-sp-%d.bam" % (k) for k in range(1,10)]
#target_names = ["mCD8T-KO-NCP-sp-%d.bam" % (k) for k in range(1,10)]
#target_names = ["HGPS-NCP-sp-%d.bam" % (k) for k in range(1,10)]
#feature_names = ['Gene activity']
feature_names = ['eigen']
#target_names = ["H1-NCP-sp-8.bam"]
#feature_names = ['ATAC score']
#feature_names = []
names = target_names + feature_names

# plot mode
#plot_mode = 'targetVSfeature'
#plot_mode = 'targetVSfeature_nospace'
#plot_mode = 'targetonly_nospace'

# chromosome cartoon aspect ratio
aspect = (0.05*bin_size) / (10**6)


#### Binning the data and draw the genome-wide plot by chromosome
for chr_choice in chr_choices:

    #### read files
    path = ""
    # read G-band information
    #chr_gID_Gband = read_Gband(path+"Gband_information.txt", chr_choices=[chr_choice]) # human
    chr_gID_Gband = read_Gband(path+"Gbanding_mouse.txt", chr_choices=[chr_choice]) # mouse

    # read annotation file
    #field_ID_value = load_file.read_tabular_file (path+"H1_NCP_sp_10kb_anot.cn", mode='col')
    #ID_pos = field_ID_value['PhysicalPosition']
    #ID_pos, field_ID_value = read_bin("mCD8T_WT-NCP_sp_10kb_bin.cn", chr_choices=[chr_choice])
    #ID_pos, field_ID_value = read_bin("mCD8T_inht-NCP_sp_10kb_bin.cn", chr_choices=[chr_choice])
    #ID_pos, field_ID_value = read_bin("mCD8T_KO-NCP_sp_bin.cn", chr_choices=[chr_choice], tlen=True)
    #ID_pos, field_ID_value = read_bin("GM_NCP_sp_10kb_bin.cn", chr_choices=[chr_choice])
    #ID_pos, field_ID_value = read_bin("H1_NCP_HP1a_10kb_bin.cn", chr_choices=[chr_choice])
    #ID_pos, field_ID_value = read_bin("H1_DNA_HP1a_10kb_bin.cn", chr_choices=[chr_choice])
    #ID_pos, field_ID_value = read_bin("HGPS_NCP_sp_bin.cn", chr_choices=[chr_choice], tlen=True)
    #ID_pos, field_ID_value = read_bin_Chalf("H1_NCP_sp_10kb_Chalf.cn", chr_choices=[chr_choice])
    #ID_pos, field_ID_value = read_bin_Chalf("GM_NCP_sp_10kb_Chalf.cn", chr_choices=[chr_choice])

    # check Chalf
    #ID_pos1, ID_Chalf1 = read_bin_Chalf("mCD8T_WT-NCP_sp_10kb_Chalf.cn", chr_choices=[chr_choice])
    #ID_pos2, ID_Chalf2 = read_bin_Chalf("mCD8T_inht-NCP_sp_10kb_Chalf.cn", chr_choices=[chr_choice])
    #ID_pos3, ID_Chalf3 = read_bin_Chalf("mCD8T_KO-NCP_sp_10kb_Chalf.cn", chr_choices=[chr_choice])

    #common_IDs = list(set(ID_pos1) & set(ID_pos2) & set(ID_pos3))
    #ID_pos = {ID:ID_pos1[ID] for ID in common_IDs}
    #field_ID_value = {'WT-Chalf':{}, 'inht-Chalf':{}, 'KO-Chalf':{}}
    #for ID in common_IDs:
    #    Chalf1 = ID_Chalf1[ID]
    #    Chalf2 = ID_Chalf2[ID]
    #    Chalf3 = ID_Chalf3[ID]
    #    field_ID_value['WT-Chalf'][ID] = Chalf1
    #    field_ID_value['inht-Chalf'][ID] = Chalf2
    #    field_ID_value['KO-Chalf'][ID] = Chalf3

    # bin num file
    #ID_pos, field_ID_value = read_bin_num("/home/spark159/../../media/spark159/sw/" +
    #                                      "mCD8T_WT-NCP_sp_10kb_num.cn", chr_choices=[chr_choice])
    #ID_pos, field_ID_value = read_bin_num("/home/spark159/../../media/spark159/sw/" +
    #                                      "mCD8T_inht-NCP_sp_10kb_num.cn", chr_choices=[chr_choice])
    ID_pos, field_ID_value = read_bin_num("/home/spark159/../../media/spark159/sw/" +
                                          "mCD8T_KO-NCP_sp_10kb_num.cn", chr_choices=[chr_choice])

        
    
    # read GTF file
    #geneID_field_values, field_geneID_values = load_file.read_GTF (path+"ENCFF159KBI.gtf", mode="both", chr_list=[chr_choice]) # human
    geneID_field_values, field_geneID_values = load_file.read_GTF (path+"gencodeM21pri-UCSC-tRNAs-ERCC-phiX.gtf", mode="both", chr_list=[chr_choice]) # mouse
    
    geneID_pos = {}
    for geneID in geneID_field_values:
        try:
            pos = geneID_field_values[geneID]['TSS']
            geneID_pos[geneID] = pos
        except:
            continue

    # read H1 tsv file
    #geneID_FPKM = load_file.read_tsv(path+"ENCFF174OMR.tsv")

    # read GM tsv file
    #geneID_FPKM = load_file.read_tsv(path + "ENCFF345SHY.tsv")


    ## read txt RNA_seq file and get geneID_FPKM (mouse CD8 T cell)
    gname_FPKM = read_txtRNA_seq("GSM3721901_W_effector_batch_1_1.txt")
    geneID_FPKM = {}
    for geneID in geneID_field_values:
        gname = geneID_field_values[geneID]["geneName"]
        try:
            FPKM = gname_FPKM[gname]
            geneID_FPKM[geneID] = FPKM
        except:
            continue

    #geneID_FPKM = load_file.read_RPKM_new ("GSE136898_rawCounts.txt", "gencodeM21pri-UCSC-tRNAs-ERCC-phiX.gtf", chr_list=[chr_choice]) 
        
    # read compartment score file
    #eigenbinID_value, eigenbinID_interval = read_eigenfile(path+"eigen_H1_100kb.txt", bin_size=100000)
    #eigenbinID_value, eigenbinID_interval = read_eigenbedgraph(path+"eigen_H1_100kb.bedgraph", chr_choice=chr_choice) # H1
    #eigenbinID_value, eigenbinID_interval = read_eigenbedgraph(path+"eigen_GM_compart.bedgraph", chr_choice=chr_choice) # GM
    eigenbinID_value, eigenbinID_interval = read_eigenbedgraph(path+"eigen_mouseCD8Tcell_100kb.bedgraph", chr_choice=chr_choice) # mouse CD8 Tcell


    #### binning the data
    # binning the annotation data
    name_binID_mean = {}
    name_binID_count = {}
    for ID in ID_pos:
        binID = int(ID_pos[ID]) / int(bin_size)
        for name in names:
            if name not in name_binID_mean:
                name_binID_mean[name] = {}
            if name not in name_binID_count:
                name_binID_count[name] = {}
            try:
                value = field_ID_value[name][ID]
            except:
                continue
            if np.isnan(value):
                continue
            if binID not in name_binID_mean[name]:
                name_binID_mean[name][binID] = 0.0
            name_binID_mean[name][binID] += value
            if binID not in name_binID_count[name]:
                name_binID_count[name][binID] = 0
            name_binID_count[name][binID] += 1

    for name in name_binID_mean:
        for binID in name_binID_mean[name]:
            name_binID_mean[name][binID] = float(name_binID_mean[name][binID]) / name_binID_count[name][binID]

    # binning the RNA-seq data
    if 'Gene activity' in names:
        temp = {}
        for ID in ID_pos:
            binID = int(ID_pos[ID]) / int(bin_size)
            temp[binID] = 0

        name_binID_mean['Gene density'] = copy.deepcopy(temp)
        name_binID_count['Gene density'] = copy.deepcopy(temp)
        name_binID_mean['Gene activity'] = copy.deepcopy(temp)
        name_binID_count['Gene activity'] = copy.deepcopy(temp)

        min_FPKM = min(geneID_FPKM.values())
        for geneID in geneID_pos:
            binID = int(geneID_pos[geneID]) / int(bin_size)
            try:
                name_binID_mean['Gene density'][binID] += 1.0
            except:
                continue
            try:
                #name_binID_mean['Gene activity'][binID] += geneID_FPKM[geneID]
                #name_binID_mean['Gene activity'][binID] += np.log2(geneID_FPKM[geneID])
                name_binID_mean['Gene activity'][binID] += np.log2(geneID_FPKM[geneID] - min_FPKM + 1)
                name_binID_count['Gene activity'][binID] += 1
            except:
                #name_binID_mean['Gene density'][binID] += 0.0
                #name_binID_mean['Gene activity'][binID] += np.nan
                continue
            #name_binID_count['Gene activity'][binID] += 1

        for binID in name_binID_mean['Gene activity']:
            if name_binID_count['Gene activity'][binID] <= 0:
                name_binID_mean['Gene activity'][binID] = np.nan
                continue
            name_binID_mean['Gene activity'][binID] /= name_binID_count['Gene activity'][binID] 

        #names.append('Gene density')
        #names.append('Gene activity')


    # binning ATAC-seq score data
    if 'ATAC score' in names:
        ATAC_fname = "H1_ATAC_%s_Bsig.cn" % (chr_choice)
        aID_pos, aID_score = read_ATACbin(ATAC_fname, chr_choices=[chr_choice])

        binID_ATAC_mean = {}
        binID_ATAC_count = {}
        for aID in aID_pos:
            binID = int(aID_pos[aID]) / int(bin_size)
            try:
                score = aID_score[aID]
            except:
                continue
            if np.isnan(score):
                continue
            if binID not in binID_ATAC_mean:
                binID_ATAC_mean[binID] = 0.0
            binID_ATAC_mean[binID] += score
            if binID not in binID_ATAC_count:
                binID_ATAC_count[binID] = 0
            binID_ATAC_count[binID] += 1

        for binID in binID_ATAC_mean:
            binID_ATAC_mean[binID] = float(binID_ATAC_mean[binID]) / binID_ATAC_count[binID]

        name_binID_mean['ATAC score'] = binID_ATAC_mean
        name_binID_count['ATAC score'] = binID_ATAC_count

        #names.append('ATAC score')

    ## binning the compartment score data
    if 'eigen' in names:    
        binID_eigens = {}
        for value, interval in zip(eigenbinID_value, eigenbinID_interval):
            st, ed = interval
            st_binID, ed_binID = st / bin_size, ed / bin_size
            for binID in range(st_binID, ed_binID):
                if binID not in binID_eigens:
                    binID_eigens[binID] = []
                binID_eigens[binID].append(value)

        binID_eigen = {}
        for binID in binID_eigens:
            try:
                binID_eigen[binID] = np.mean(binID_eigens[binID])
            except:
                binID_eigen[binID] = np.nan
        name_binID_mean['eigen'] = binID_eigen


    # binning the G-banding data and make ideogram
    gID_Gband = chr_gID_Gband[chr_choice]
    binID_Gvalue, binID_Gtype = {}, {}
    gID_binterval = {}
    gband_img, gband_cenimg, gband_varimg = [], [], []
    gtick_locs, gtick_labels = [], []
    for gID in sorted(gID_Gband.keys()):
        st, ed = gID_Gband[gID]['interval']
        st_binID, ed_binID = st / bin_size, ed / bin_size
        gtype = gID_Gband[gID]['type']
        gvalue = gID_Gband[gID]['value']
        for binID in range(st_binID, ed_binID):
            binID_Gvalue[binID] = gID_Gband[gID]['value']
            binID_Gtype[binID] = gID_Gband[gID]['type']
        gID_binterval[gID] = (st_binID, ed_binID)
        bsize = ed_binID - st_binID
        gband_img += [[gvalue] for k in range(bsize)] 
        if gtype == 'acen':
            gband_cenimg += [[10] for k in range(bsize)]
        else:
            gband_cenimg += [[np.nan] for k in range(bsize)]
        if gtype == 'var':
            gband_varimg += [[10] for k in range(bsize)]
        else:
            gband_varimg += [[np.nan] for k in range(bsize)]
        mid = (st_binID + ed_binID)/2
        gname = gID_Gband[gID]['name']
        gtick_locs.append(mid)
        gtick_labels.append(gname)


    # set xtick labels along chromosome
    xtick_locs, xtick_labels = [], []
    for binID in sorted(binID_Gvalue.keys()):
        pos = bin_size*binID + bin_size/2
        Mb_pos = int(round(float(pos)/(10**6)))
        label = str(Mb_pos)
        if label not in xtick_labels:
            xtick_locs.append(binID)
            xtick_labels.append(label)

            
    # smoothing the binned data by sliding window average
    name_sig = {}
    for name in names:
        binID_mean = name_binID_mean[name]
        sig = []
        for binID in range(len(gband_img)):
            #size = binID_size[binID]
            size = 1
            try:
                sig += [binID_mean[binID] for k in range(size)]
            except:
                sig += [np.nan for k in range(size)]
        if name == 'eigen':
            sig = statis.slow_moving_average2(sig, int(4*i/5.0 + 1))
        else:
            sig = statis.slow_moving_average2(sig, blur_win)
        name_sig[name] = sig


    # temporal titration plot
    if True:
        Y_list = []
        sig_len = len(name_sig[1])
        for i in range(sig_len):
            Y = [1.0]
            valid = True
            for j in range(1, 10):        
                value = name_sig[j][i]
                if np.isnan(value):
                    valid = False
                    break
                Y.append(value)
            if valid:
                Y_list.append(Y)

        X = [tnum_conc[i] for i in range(10)]

        fig = plt.figure(figsize=(2.4,1.8))
        for Y in Y_list:
            try:
                p0 = [max(Y), np.median(X), 1]
                bounds = ([0.0, 0.0, 0.0], [max(Y)+max(Y)*0.1, np.inf, np.inf])
                popt, pcov = curve_fit(sigmoid, X, Y, p0, bounds = bounds,  method='dogbox')
                residuals = np.asarray(Y)- sigmoid(X, *popt)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((np.asarray(Y)-np.mean(Y))**2)
                r_squared = 1 - (ss_res / ss_tot)
                pred_X = np.linspace(min(X), max(X), 1000)
                pred_Y = sigmoid(pred_X, *popt)
                plt.plot(pred_X, pred_Y, 'k-', alpha=0.1)
                #plt.plot(X, Y, '.-', color='k', alpha=0.1)
            except:
                pass
        plt.xscale('log', basex=10)
        plt.ylim([0, 1])
        plt.xlabel("Spermine concentration (mM)", fontsize=8)
        plt.ylabel("Soluble fraction", fontsize=8)
        plt.tick_params(axis='both', which='major', labelsize=5)
        plt.tick_params(axis='both', which='minor', labelsize=5)
        plt.savefig("smoothed_10kb_curve.png", dpi=500, bbox_inches='tight')
        #plt.show()
        plt.close()
            

    # plot genome cartoon and genome-wide data (target VS feature)
    if True:
        for i in range(len(feature_names)):
            fig = plt.figure(figsize=(15,5))
            ax1 = plt.subplot(211)
            ax2 = plt.subplot(212)

            cmap = mpl.cm.get_cmap("jet")
            color_list = np.linspace(0.01, 0.99, num=len(target_names))
            for j in range(len(target_names)):
                target_name = target_names[j]
                #ax1.plot(name_sig[target_name], '#1f77b4', alpha=1)
                #ax1.plot(name_sig[target_name], alpha=0.5, lw=2, label=str(j+1))
                #ax1.plot(name_sig[target_name], alpha=1, lw=2, label=target_name)
                ax1.plot(name_sig[target_name], alpha=1, lw=2, color=cmap(color_list[j]),
                         label='[sp]=%.2fmM' % (tnum_conc[target_name]))

            for gID in sorted(gID_binterval):
                st, ed = gID_binterval[gID]
                if gID_Gband[gID]['type'] == 'pos':
                    ax1.axvspan(st, ed-1, alpha=0.15, color='grey')

            ax1.set_xticks(xtick_locs[::10])
            ax1.set_xticklabels(xtick_labels[::10])
            ax1.set_xlabel("Position (Mb)")
            #ax1.set_ylabel(target_name, color='blue')
            ax1.set_ylabel('Survival Probability', color='blue')
            ax1.tick_params('y', colors='blue')
            ax1.set_xlim([0, len(gband_img)+1])
            #ax1.set_ylim([-0.4, 0.6])
            #ax1.set_ylim([-1.2, 0])
            #ax1.set_ylim([-0.8, -0.2])
            #ax1.set_ylim([-0.5, 0.5])
            #ax1.set_ylim([-0.5, 1.5])
            #ax1.set_ylim([-0.5, 2.0])
            #ax1.set_ylim([-0.5, 1.0])
            #ax1.set_ylim([0.14, 0.22])
            #ax1.set_ylim([0.2, 3.5])
            #ax1.set_ylim([0, 1.0])
            ax1.set_ylim([2**-4.5, 2**-0.5])
            ax1.set_yscale('log', basey=2)
            #ax1.legend()
            ax1.legend(facecolor='white', framealpha=1) # temporal


            ax1p = ax1.twinx()
            feature_name = feature_names[i]
            #ax1p.plot(name_sig[feature_name], '#d62728', alpha=0.8, lw=2)
            #ax1p.plot(name_sig[feature_name], 'k', alpha=0.8, lw=2)
            ax1p.plot(name_sig[feature_name], 'k--', alpha=0.5, lw=1.5)
            #ax1p.plot(name_sig[feature_name], 'tab:orange', alpha=0.8)
            #ax1p.plot(np.log(-1*np.asarray(name_sig[feature_name])), '#d62728', alpha=0.5)
            #ax1p.set_ylabel(feature_name, color='r')
            ax1p.set_ylabel(feature_name, color='k')
            #ax1p.set_ylabel('Eigenvector', color='orangered')
            #ax1p.tick_params('y', colors='#d62728')
            #ax1p.tick_params('y', colors='orangered')
            #ax1p.set_ylim([-0.1, 2.0])
            #ax1p.set_ylim([-0.2, 3.0])
            #ax1p.set_ylim([0.3, 0.85])

            ax2.imshow(np.transpose(gband_img), cmap='Greys', aspect=0.3/aspect)
            ax2.imshow(np.transpose(gband_cenimg), cmap ='Reds', vmin=0, vmax=20, aspect=0.3/aspect)
            ax2.imshow(np.transpose(gband_varimg), cmap ='Purples', vmin=0, vmax=20, aspect=0.3/aspect)
            ax2.set_yticks([])
            ax2.set_xticks(gtick_locs)
            ax2.set_xticklabels(gtick_labels)
            ax2.tick_params(axis="x", labelsize=5, rotation=90)
            ax2.set_xlim([0, len(gband_img)+1])
            
            plt.tight_layout()
            #plt.savefig("Gwide_" + chr_choice + '_' + target_name + '_' + feature_name + ".png", bbox_inches='tight', dpi=1000)
            plt.savefig("Gwide_" + chr_choice + '_' + str(target_name) + '_' + feature_name + ".svg", format='svg', bbox_inches='tight')
            #plt.show()
            plt.close()


    # plot genome cartoon and genome-wide data (no space, horizontal)
    if False:
        figwidth = round(10*float(len(gband_img))/9970, 1)
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(figwidth,1.5), sharex=True, gridspec_kw = {'hspace':0.02, 'wspace':None, 'height_ratios':[6,1]})

        for j in range(len(target_names)):
            target_name = target_names[j]
            axes[0].plot(name_sig[target_name], '#1f77b4', alpha=1)

        ax0p = axes[0].twinx()
        for j in range(len(feature_names)):
            feature_name = feature_names[j]
            ax0p.plot(name_sig[feature_name], '#d62728', alpha=0.25)

        for gID in sorted(gID_binterval):
            gtype = gID_Gband[gID]['type']
            if gtype == 'acen' or gtype == 'var':
                st, ed = gID_binterval[gID]
                axes[0].axvspan(st, ed-1, alpha=1, color='white', zorder=10)
                ax0p.axvspan(st, ed-1, alpha=1, color='white', zorder=10)


        axes[0].set_ylabel(chr_choice, rotation=90+40, fontsize=12, labelpad=10)
        axes[0].set_ylim([-0.5, 0.6])
        axes[0].spines['top'].set_visible(False)
        axes[0].spines['bottom'].set_visible(False)
        axes[0].spines['left'].set_visible(False)
        axes[0].spines['right'].set_visible(False)
        axes[0].tick_params(top='off', bottom='off', left='off', right='off', labelleft='off', labelbottom='off')

        ax0p.set_ylim([-0.1, 2.0])
        ax0p.spines['top'].set_visible(False)
        ax0p.spines['bottom'].set_visible(False)
        ax0p.spines['left'].set_visible(False)
        ax0p.spines['right'].set_visible(False)
        ax0p.tick_params(top='off', bottom='off', left='off', right='off', labelright='off', labelbottom='off')

        axes[1].imshow(np.transpose(gband_img), cmap='Greys', aspect='auto')
        axes[1].imshow(np.transpose(gband_cenimg), cmap ='Reds', vmin=0, vmax=20, aspect='auto')
        axes[1].imshow(np.transpose(gband_varimg), cmap ='Purples', vmin=0, vmax=20, aspect='auto')
        axes[1].set_xticks([])
        axes[1].set_yticks([])
        
        #plt.tight_layout()
        #plt.savefig("Gwide_" + chr_choice + '_' + target_name + '_' + feature_name + ".png", bbox_inches='tight', dpi=500)
        #plt.xticks(rotation=90)
        plt.show()
        plt.close()

    # plot genome cartoon and genome-wide data (no space, vertical)
    if False:
        figheight = round(8*float(len(gband_img))/9970, 1)
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(2, figheight), sharey=True, gridspec_kw = {'hspace':0.02, 'wspace':0.02, 'width_ratios':[1,8]})

        for j in range(len(target_names)):
            target_name = target_names[j]
            axes[1].plot(name_sig[target_name], range(len(name_sig[target_name])), '#1f77b4', alpha=1)

        ax1p = axes[1].twiny()
        for j in range(len(feature_names)):
            feature_name = feature_names[j]
            ax1p.plot(name_sig[feature_name], range(len(name_sig[feature_name])), '#d62728', alpha=0.25)

        for gID in sorted(gID_binterval):
            gtype = gID_Gband[gID]['type']
            if gtype == 'acen' or gtype == 'var':
                st, ed = gID_binterval[gID]
                axes[1].axhspan(st, ed-1, alpha=1, color='white', zorder=10)
                ax1p.axhspan(st, ed-1, alpha=1, color='white', zorder=10)

        axes[1].set_xlim([-0.5, 0.6])
        axes[1].spines['top'].set_visible(False)
        axes[1].spines['bottom'].set_visible(False)
        axes[1].spines['left'].set_visible(False)
        axes[1].spines['right'].set_visible(False)
        axes[1].tick_params(top='off', bottom='off', left='off', right='off', labelbottom='off')

        ax1p.set_xlabel(chr_choice, rotation=40, fontsize=15, labelpad=23, ha='right')
        ax1p.set_xlim([-0.1, 2.0])
        ax1p.spines['top'].set_visible(False)
        ax1p.spines['bottom'].set_visible(False)
        ax1p.spines['left'].set_visible(False)
        ax1p.spines['right'].set_visible(False)
        ax1p.tick_params(top='off', bottom='off', left='off', right='off', labeltop='off')

        axes[0].imshow(gband_img, cmap='Greys', aspect='auto')
        axes[0].imshow(gband_cenimg, cmap ='Reds', vmin=0, vmax=20, aspect='auto')
        axes[0].imshow(gband_varimg, cmap ='Purples', vmin=0, vmax=20, aspect='auto')
        axes[0].set_xticks([])
        axes[0].set_yticks([])
        
        #plt.tight_layout()
        #plt.savefig("Gwide_" + chr_choice + '_' + target_name + '_' + feature_name + ".png", bbox_inches='tight')
        plt.savefig("Gwide_" + chr_choice + '_' + target_name + '_' + feature_name + ".svg", format='svg', bbox_inches='tight')

        #plt.show()
        plt.close()
