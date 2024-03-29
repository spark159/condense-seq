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

def read_bin_num (fname, skip_star=False, chr_choices=None):
    First = True
    ID_pos = {}
    name_ID_count = {}
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
            if skip_star:
                continue
            cols = cols[:-1]
        ID, chr, start, end = cols[:4]
        start, end = int(start), int(end)
        if chr_choices!=None and chr not in chr_choices:
            continue
        pos = int(float(start + end)/2)
        assert ID not in ID_pos
        ID_pos[ID] = pos
        counts = [int(count) for count in cols[4:]]
        for name, count in zip(names, counts):
            if name not in name_ID_count:
                name_ID_count[name] = {}
            assert ID not in name_ID_count[name]
            name_ID_count[name][ID] = count
    return ID_pos, name_ID_count


def read_bin_file (fname, chr_choices=None):
    First = True
    ID_pos = {}
    name_ID_value = {}
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        if First:
            names = []
            for name in cols[4:]:
                #agent = name.rsplit('.', 1)[0].split('-')[-2]
                #tnum = int(name.rsplit('.', 1)[0].split('-')[-1])
                #names.append(tnum)
                names.append(name.rsplit('.', 1)[0])
            First = False
            continue
        ID, chr, start, end = cols[:4]
        start, end = int(start), int(end)
        if chr_choices!=None and chr not in chr_choices:
            continue
        pos = int(float(start + end)/2)
        assert ID not in ID_pos
        ID_pos[ID] = pos
        values = [float(value) for value in cols[4:]]
        for name, value in zip(names, values):
            if name not in name_ID_value:
                name_ID_value[name] = {}
            assert ID not in name_ID_value[name]
            name_ID_value[name][ID] = value
    return ID_pos, name_ID_value

def read_bin_Chalf (fname, chr_choices=None):
    ID_pos = {}
    ID_value = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        if First:
            First = False
            continue
        binID, chr, st, ed, Chalf = cols
        if chr_choices != None and chr not in chr_choices:
            continue
        st, ed = int(st), int(ed)
        ID = (st, ed)
        st, ed = int(st), int(ed)
        Chalf = float(Chalf)
        pos = int(float(st + ed)/2)
        ID_pos[ID] = pos
        ID_value[ID] = Chalf
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


# parameters
# set path
#path = ""
path = "/home/spark159/../../media/spark159/sw/"

agent_fullname = {'sp':'Spermine(4+)',
                  'spd':'Spermidine(3+)',
                  'CoH':'Cobalt Hexammine(3+)',
                  'PEG':'PEG 8000',
                  'HP1a':'HP1$\\alpha$',
                  'HP1bSUV':'HP1$\\beta$+tSUV',
                  'LKH':'Linker histone1',
                  'Ki67':'Ki67',
                  'FUS':'FUS',
                  'Mg':'Magnesium',
                  'Ca':'Calcium'}

cell = 'H1'
# experiment list (cell, sample, agent)
# it should be same cell line
#exp_list = [(cell, 'NCP', 'sp'),
#            (cell, 'NCP', 'spd'),
#            (cell, 'NCP', 'CoH'),
#            (cell, 'NCP', 'PEG'),
#            (cell, 'NCP', 'Ca'),
#            (cell, 'NCP', 'Mg'),
#            (cell, 'NCP', 'HP1a'),
#            (cell, 'NCP', 'HP1bSUV'),
#            (cell, 'NCP', 'LKH'),
#            (cell, 'NCP', 'Ki67'),
#            (cell, 'NCP', 'FUS')]

#exp_list = [(cell, 'NCP', 'sp'),
#            (cell, 'NCP', 'spd'),
#            (cell, 'NCP', 'CoH'),
#            (cell, 'NCP', 'PEG'),
#            (cell, 'NCP', 'Mg'),
#            (cell, 'NCP', 'Ca'),
#            (cell, 'NCP', 'HP1a')]

exp_list = [(cell, 'NCP', 'sp', 8),
            (cell, 'NCP', 'spd', 6),
            (cell, 'NCP', 'CoH', 5),
            (cell, 'NCP', 'PEG', 6),
            (cell, 'NCP', 'Ca', 5),
            (cell, 'NCP', 'Mg', 5),
            (cell, 'NCP', 'HP1a', 3),
            (cell, 'NCP', 'HP1bSUV', 4),
            (cell, 'NCP', 'LKH', 3),
            (cell, 'NCP', 'Ki67', 4),
            (cell, 'NCP', 'FUS', 5)]

#exp_list = [(cell, 'NCP', 'HP1a')]

# set species and gender
if cell in ['H1', 'GM']:
    species = 'human'
elif cell in ['mCD8T']:
    species = 'mouse'

if cell in ['H1']:
    gender = 'male'
elif cell in ['GM', 'mCD8T']:
    gender = 'female'

#### read G-band information
if species == 'human':
    Gband_fname = "Gband_information.txt"
elif species == 'mouse':
    Gband_fname = "Gbanding_mouse.txt"
chr_gID_Gband = read_Gband(path+"Gband_information.txt")

#### read titration file
for cell, sample, agent, tnum in exp_list:
    tnum_conc, tnum_frac = load_file.read_titration("%s_%s_%s_titration.csv" % (cell, sample, agent))

#### read RNA-seq file
if cell == 'H1':
    tsv_fname = "ENCFF174OMR.tsv"
    geneID_FPKM = load_file.read_tsv(path+tsv_fname)
elif cell == 'GM':
    tsv_fname = "ENCFF345SHY.tsv"
    geneID_FPKM = load_file.read_tsv(path+tsv_fname)
elif cell == 'mCD8T':
    ## read txt RNA_seq file and get geneID_FPKM (mouse CD8 T cell)
    txt_fname = "GSM3721901_W_effector_batch_1_1.txt"
    gname_FPKM = read_txtRNA_seq(path + txt_fname)
    geneID_FPKM = {}
    for geneID in geneID_field_values:
        gname = geneID_field_values[geneID]["geneName"]
        try:
            FPKM = gname_FPKM[gname]
            geneID_FPKM[geneID] = FPKM
        except:
            continue


#### set parameters
# set binning resolution
i = 20
bin_size = int(0.5*(10**6) / i) # binsize (unit of bp)
blur_win = int(4*i + 1) # sliding window (unit of bin)

# set chromosomes
chr_choices = ['chr1']

# set target names and feature names
#target_names = [k for k in range(1,10)]
#target_names = ['Chalf']

target_names = ['-'.join([cell, sample, agent, str(tnum)]) for cell, sample, agent, tnum in exp_list]

#target_names = []
#for cell, sample, agent in exp_list:
#    target_names.append("%s-%s-%s-Chalf" % (cell, sample, agent))
feature_names = ['Gene activity']
#feature_names = ['eigen']
names = target_names + feature_names

# plot mode
#plot_mode = 'targetVSfeature'
#plot_mode = 'targetVSfeature_nospace'
#plot_mode = 'targetonly_nospace'

#plot_mode = 'targetVSfeature'
#plot_mode = 'targetByfeature'
plot_mode = 'featureBytarget'

# target file information
target_ftype = "zscore"
#target_ftype = "Chalf"
#target_ftype = "zChalf"
target_binsize = 10000

# chromosome cartoon aspect ratio
aspect = (0.05*bin_size) / (10**6)

# graphics parameters
#cmap = mpl.cm.get_cmap("Blues")
#color_list = np.linspace(0.1, 0.9, num=len(target_names))
#target_colors = [cmap(color_list[k]) for k in range(len(target_names))]
#target_labels = ['[%s]=%.2fmM' % (agent, tnum_conc[tnum])
#                 for cell, sample, agent, tnum in target_names]
#target_ylims = [-1.2, 1.5]
#target_ylims = [None, None]
#target_ylabel = "Condensability"
#target_ycolor = "blue"

#target_colors = ['tab:blue']
#target_labels = [None]
#target_ylims = [None, None]
#target_ylabel = "C-half"

target_colors = [None for k in range(len(target_names))]
target_labels = [agent_fullname[agent] for cell, sample, agent, tnum in exp_list]
target_ylims = [-1.2, 1.2]
#target_ylims = [None, None]
target_ylabel = "zscore"
target_ycolor = "blue"

target_alpha = 0.8
target_lw = 1.5
target_linestyle = None


#target_ycolor = "black"
legend_loc = 'upper right'

feature_colors = ["tab:red"]
feature_labels = [None]
feature_ylims = [None, None]
feature_alpha = 0.7
feature_lw = 1
feature_linestyle = None
feature_ylabel = "Gene activity"
feature_ycolor = "red"

#feature_colors = ["black"]
#feature_labels = [None]
#feature_ylims = [None, None]
#feature_alpha = 0.7
#feature_lw = 1
#feature_linestyle = '--'
#feature_ylabel = "A/B compartment score"
#feature_ycolor = "black"

note =""

#### Binning the data and draw the genome-wide plot by chromosome
for chr_choice in chr_choices:

    #### read and binning the data
    name_binID_mean = {}
    name_binID_count = {}

    ## read and binning target data
    # read target data
    name_ID_pos = {}
    name_ID_value = {}
    
    for cell, sample, agent, tnum in exp_list:

        target_fname = '_'.join([cell, sample, agent,
                                 str(int(target_binsize/1000.0)) + 'kb',
                                 target_ftype]) + '.cn'

        if target_ftype == 'anot':
           field_ID_value = load_file.read_tabular_file (path + target_fname, mode='col')
           ID_pos = field_ID_value['PhysicalPosition']

        elif target_ftype in ['num', 'score', 'zscore']:
            ID_pos, field_ID_value = read_bin_file(path + target_fname, chr_choices=[chr_choice])

        elif target_ftype in ['Chalf', 'zChalf']:
            ID_pos, ID_value = read_bin_Chalf(path + target_fname, chr_choices=[chr_choice])
            field_ID_value = {}
            field_ID_value['%s-%s-%s-Chalf' % (cell, sample, agent)] = ID_value

        for name in names:
            try:
                ID_value = field_ID_value[name]
            except:
                continue
            if name not in name_ID_pos:
                name_ID_pos[name] = {}
            name_ID_pos[name].update(ID_pos)
            if name not in name_ID_value:
                name_ID_value[name] = {}
            name_ID_value[name].update(ID_value)
    
    # binning the target data
    for name in name_ID_pos:
        ID_pos = name_ID_pos[name]
        ID_value = name_ID_value[name]
        for ID in ID_pos:
            binID = int(ID_pos[ID]) / int(bin_size)
            value = ID_value[ID]
            if np.isnan(value):
                continue
            if name not in name_binID_mean:
                name_binID_mean[name] = {}
            if binID not in name_binID_mean[name]:
                name_binID_mean[name][binID] = 0.0
            name_binID_mean[name][binID] += value
            if name not in name_binID_count:
                name_binID_count[name] = {}
            if binID not in name_binID_count[name]:
                name_binID_count[name][binID] = 0
            name_binID_count[name][binID] += 1

    for name in name_binID_mean:
        for binID in name_binID_mean[name]:
            name_binID_mean[name][binID] = float(name_binID_mean[name][binID]) / name_binID_count[name][binID]

    ## read and binning the RNA-seq data
    if 'Gene activity' in names:

        # read GTF file
        if species == 'human':
            gtf_fname = "ENCFF159KBI.gtf"
        elif species == 'mouse':
            gtf_fname = "gencodeM21pri-UCSC-tRNAs-ERCC-phiX.gtf"
                
        geneID_field_values, field_geneID_values = load_file.read_GTF (path + gtf_fname,
                                                                       mode="both",
                                                                       chr_list=[chr_choice])

        geneID_pos = {}
        for geneID in geneID_field_values:
            try:
                pos = geneID_field_values[geneID]['TSS']
                geneID_pos[geneID] = pos
            except:
                continue
        
        min_FPKM = min(geneID_FPKM.values())

        FPKM_binID_mean = {}
        FPKM_binID_count = {}
        for geneID in geneID_FPKM:
            try:
                binID = int(geneID_pos[geneID]) / int(bin_size)
            except:
                continue

            if binID not in FPKM_binID_mean:
                FPKM_binID_mean[binID] = 0.0

            ex_level = np.log2(geneID_FPKM[geneID] - min_FPKM + 1) #rescaling
            FPKM_binID_mean[binID] += ex_level

            if binID not in FPKM_binID_count:
                FPKM_binID_count[binID] = 0
            FPKM_binID_count[binID] += 1

        for binID in FPKM_binID_mean:
            FPKM_binID_mean[binID] = float(FPKM_binID_mean[binID]) / FPKM_binID_count[binID]

        name_binID_mean['Gene activity'] = FPKM_binID_mean
        name_binID_count['Gene activity'] = FPKM_binID_count



    ## read and binning ATAC-seq score data
    if 'ATAC score' in names:
        # read ATAC file
        ATAC_fname = "H1_ATAC_%s_Bsig.cn" % (chr_choice)
        aID_pos, aID_score = read_ATACbin(ATAC_fname, chr_choices=[chr_choice])

        ATAC_binID_mean = {}
        ATAC_binID_count = {}
        for aID in aID_pos:
            binID = int(aID_pos[aID]) / int(bin_size)
            try:
                score = aID_score[aID]
            except:
                continue
            if np.isnan(score):
                continue
            if binID not in ATAC_binID_mean:
                ATAC_binID_mean[binID] = 0.0
            ATAC_binID_mean[binID] += score
            if binID not in ATAC_binID_count:
                ATAC_binID_count[binID] = 0
            ATAC_binID_count[binID] += 1

        for binID in ATAC_binID_mean:
            ATAC_binID_mean[binID] = float(ATAC_binID_mean[binID]) / ATAC_binID_count[binID]

        name_binID_mean['ATAC score'] = ATAC_binID_mean
        name_binID_count['ATAC score'] = ATAC_binID_count


    ## read and binning the compartment score data
    if 'eigen' in names:
        # read compartment score file
        if cell == 'H1':
            eigen_fname = "eigen_H1_100kb.bedgraph"
        elif cell == 'GM':
            eigen_fname = "eigen_GM_compart.bedgraph"
        elif cell == 'mCD8T':
            eigen_fname = "eigen_mouseCD8Tcell_100kb.bedgraph"
        eigenbinID_value, eigenbinID_interval = read_eigenbedgraph(path+eigen_fname,
                                                                   chr_choice=chr_choice)
  
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


    ## binning the G-banding data and make ideogram
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


    #### set xtick labels along chromosome
    xtick_locs, xtick_labels = [], []
    for binID in sorted(binID_Gvalue.keys()):
        pos = bin_size*binID + bin_size/2
        Mb_pos = int(round(float(pos)/(10**6)))
        label = str(Mb_pos)
        if label not in xtick_labels:
            xtick_locs.append(binID)
            xtick_labels.append(label)

            
    #### smoothing the binned data by sliding window average
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
            

    #### plot genome cartoon and genome-wide data
    ## plot targets by feature
    if plot_mode == 'targetByfeature':
        for i in range(len(feature_names)):
            fig = plt.figure(figsize=(15,5))
            ax1 = plt.subplot(211)
            ax2 = plt.subplot(212)

            for j in range(len(target_names)):
                target_name = target_names[j]
                ax1.plot(name_sig[target_name],
                         color=target_colors[j],
                         alpha=target_alpha,
                         lw=target_lw,
                         linestyle=target_linestyle,
                         label=target_labels[j])

            for gID in sorted(gID_binterval):
                st, ed = gID_binterval[gID]
                if gID_Gband[gID]['type'] == 'pos':
                    ax1.axvspan(st, ed-1, alpha=0.15, color='grey')

            ax1.set_xticks(xtick_locs[::10])
            ax1.set_xticklabels(xtick_labels[::10])
            ax1.set_xlabel("Position (Mb)")
            ax1.set_ylabel(target_ylabel, color=target_ycolor)
            ax1.tick_params('y', colors=target_ycolor)
            ax1.set_xlim([0, len(gband_img)+1])
            ax1.set_ylim(target_ylims)
            leg = ax1.legend(framealpha=1, loc=legend_loc)
            for legobj in leg.legendHandles:
                legobj.set_alpha(1.0)
                legobj.set_linewidth(2.0)

            ax1p = ax1.twinx()
            feature_name = feature_names[i]
            ax1p.plot(name_sig[feature_name],
                      color=feature_colors[i],
                      alpha=feature_alpha,
                      lw=feature_lw,
                      linestyle=feature_linestyle)
            ax1p.set_ylabel(feature_ylabel, color=feature_ycolor)
            ax1p.tick_params('y', colors=feature_ycolor)
            ax1p.set_ylim(feature_ylims)

            ax1.set_zorder(1)  # default zorder is 0 for ax1 and ax2
            ax1.patch.set_visible(False)  # prevents ax1 from hiding ax2

            ax2.imshow(np.transpose(gband_img), cmap='Greys', aspect=0.3/aspect)
            ax2.imshow(np.transpose(gband_cenimg), cmap ='Reds', vmin=0, vmax=20, aspect=0.3/aspect)
            ax2.imshow(np.transpose(gband_varimg), cmap ='Purples', vmin=0, vmax=20, aspect=0.3/aspect)
            ax2.set_yticks([])
            ax2.set_xticks(gtick_locs)
            ax2.set_xticklabels(gtick_labels)
            ax2.tick_params(axis="x", labelsize=5, rotation=90)
            ax2.set_xlim([0, len(gband_img)+1])
            
            plt.tight_layout()
            plt.savefig("_".join(['Gwide', chr_choice, str(target_name), str(feature_name), note])
                        + ".svg",
                        format='svg',
                        bbox_inches='tight')
            #plt.show()
            plt.close()


    ## plot feature by target
    if plot_mode == 'featureBytarget':
        for i in range(len(target_names)):
            fig = plt.figure(figsize=(15,5))
            ax1 = plt.subplot(211)
            ax2 = plt.subplot(212)

            target_name = target_names[i]
            ax1.plot(name_sig[target_name],
                     color=target_colors[i],
                     alpha=target_alpha,
                     lw=target_lw,
                     linestyle=target_linestyle,
                     label=target_labels[i])
            ax1.set_ylabel(target_ylabel, color=target_ycolor)
            ax1.tick_params('y', colors=target_ycolor)
            ax1.set_ylim(target_ylims)

            for gID in sorted(gID_binterval):
                st, ed = gID_binterval[gID]
                if gID_Gband[gID]['type'] == 'pos':
                    ax1.axvspan(st, ed-1, alpha=0.15, color='grey')

            ax1.set_xticks(xtick_locs[::10])
            ax1.set_xticklabels(xtick_labels[::10])
            ax1.set_xlabel("Position (Mb)")
            ax1.set_ylabel(target_ylabel, color=target_ycolor)
            ax1.tick_params('y', colors=target_ycolor)
            ax1.set_xlim([0, len(gband_img)+1])
            ax1.set_ylim(target_ylims)
            leg = ax1.legend(framealpha=1, loc=legend_loc)
            for legobj in leg.legendHandles:
                legobj.set_alpha(1.0)
                legobj.set_linewidth(2.0)


            ax1p = ax1.twinx()
            for j in range(len(feature_names)):
                feature_name = feature_names[j]
                ax1p.plot(name_sig[feature_name],
                          color=feature_colors[j],
                          alpha=feature_alpha,
                          lw=feature_lw,
                          linestyle=feature_linestyle)

            ax1.set_zorder(1)  # default zorder is 0 for ax1 and ax2
            ax1.patch.set_visible(False)  # prevents ax1 from hiding ax2

            ax2.imshow(np.transpose(gband_img), cmap='Greys', aspect=0.3/aspect)
            ax2.imshow(np.transpose(gband_cenimg), cmap ='Reds', vmin=0, vmax=20, aspect=0.3/aspect)
            ax2.imshow(np.transpose(gband_varimg), cmap ='Purples', vmin=0, vmax=20, aspect=0.3/aspect)
            ax2.set_yticks([])
            ax2.set_xticks(gtick_locs)
            ax2.set_xticklabels(gtick_labels)
            ax2.tick_params(axis="x", labelsize=5, rotation=90)
            ax2.set_xlim([0, len(gband_img)+1])
            
            plt.tight_layout()
            plt.savefig("_".join(['Gwide', chr_choice, str(target_name), str(feature_name), note])
                        + ".svg",
                        format='svg',
                        bbox_inches='tight')
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
