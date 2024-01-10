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
                #agent = name.rsplit('.', 1)[0].split('-')[-2]
                #tnum = int(name.rsplit('.', 1)[0].split('-')[-1])
                names.append(name.rsplit('.', 1)[0]) 
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

def read_bin_fract (fname, skip_star=False, chr_choices=None):
    First = True
    ID_pos = {}
    name_ID_fract = {}
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        if First:
            names = []
            for name in cols[4:-2]:
                #agent = name.rsplit('.', 1)[0].split('-')[-2]
                #tnum = int(name.rsplit('.', 1)[0].split('-')[-1])
                names.append(name.rsplit('.', 1)[0]) 
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
        control = int(cols[-1])
        if control <=0:
            continue
        pos = int(float(start + end)/2)
        assert ID not in ID_pos
        ID_pos[ID] = pos
        fracts = [float(count)/control for count in cols[4:-1]]
        for name, fract in zip(names, fracts):
            if name not in name_ID_fract:
                name_ID_fract[name] = {}
            assert ID not in name_ID_fract[name]
            name_ID_fract[name][ID] = fract
    return ID_pos, name_ID_fract



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
#path = "/home/spark159/../../media/spark159/sw/"
path = "/home/spark159/../../storage/"
path2 = "/home/spark159/../../storage/replicates/"


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

exp_list = [(cell, 'NCP', 'sp', 8),
            (cell, 'NCP', 'spd', 6),
            (cell, 'NCP', 'CoH', 5),
            (cell, 'NCP', 'PEG', 6),
            (cell, 'NCP', 'Ca', 5),
            (cell, 'NCP', 'HP1a', 3),
            (cell, 'NCP', 'HP1bSUV', 4)]

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

#exp_list = [(cell, 'NCP', 'sp', 8),
#            (cell, 'NCP', 'spd', 6),
#            (cell, 'NCP', 'CoH', 5),
#            (cell, 'NCP', 'PEG', 6),
#            (cell, 'NCP', 'Ca', 5),
#            (cell, 'NCP', 'Mg', 5),
#            (cell, 'NCP', 'HP1a', 3),
#            (cell, 'NCP', 'HP1bSUV', 4),
#            (cell, 'NCP', 'LKH', 3),
#            (cell, 'NCP', 'Ki67', 4),
#            (cell, 'NCP', 'FUS', 5)]

cell = 'mCD8T'
#exp_list = [(cell, 'WT-NCP', 'sp', 8),
#            (cell, 'inht-NCP', 'sp', 8),
#            (cell, 'KO-NCP', 'sp', 8)]
#exp_list = [(cell, 'WT-NCP', 'sp', i) for i in range(1, 10)]
#exp_list = [(cell, 'inht-NCP', 'sp', i) for i in range(1, 10)]
exp_list = [(cell, 'KO-NCP', 'sp', i) for i in range(1, 10)]

#exp_list = [(cell, 'NCP', 'HP1bSUV', i) for i in range(1, 6)]


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
chr_gID_Gband = read_Gband(path+Gband_fname)

#### read titration file
for cell, sample, agent, tnum in exp_list:
    tnum_conc, tnum_frac = load_file.read_titration("%s_%s_%s_titration.csv" % (cell, sample, agent))

#### set parameters
# set binning resolution
i = 20
bin_size = int(0.5*(10**6) / i) # binsize (unit of bp)
blur_win = int(4*i + 1) # sliding window (unit of bin)

# set chromosomes
if species == 'human':
    chr_choices = ['chr' + str(i) for i in range(1, 23)]
elif species == 'mouse':
    chr_choices = ['chr' + str(i) for i in range(1, 20)]
chr_choices += ['chrX']
if gender == 'male':
    chr_choices += ['chrY']
    
chr_choices = ['chr1']


# set target names and feature names
target_names = ['-'.join([cell, sample, agent, str(tnum)]) for cell, sample, agent, tnum in exp_list]
feature_names = ['Gene activity']
names = target_names + feature_names

# target file information
target_ftype = "score"
#target_ftype = "Chalf"
#target_ftype = "zChalf"
#target_ftype = "fract"
target_binsize = 10000

# chromosome cartoon aspect ratio
aspect = (0.05*bin_size) / (10**6)

# graphics parameters
alphas = [0.8, 0.4]
lws = [2, 1.5]
linestyles = [None, None]
ylabels = ['Condensability', 'Gene activity']
ylabels = ["", ""]
ycolors = ['Blue', 'Red']
#ycolors = [None, None]
#ylims_list = [[-1.2, 1.2], [None, None]]
#ylims_list = [[0.0, 1.0], [None, None]]
ylims_list = [[None, None], [None, None]]
yscales = [None, None] #To dos
legend_loc = 'upper right'
shaded = False

# optional graphics parameters
#name_color = {'Gene activity':'tab:red'}
name_color = {'Gene activity':'black'}
#cmap = mpl.cm.get_cmap("jet")
#color_list = np.linspace(0.01, 0.99, num=len(target_names))
#name_color = {target_names[i]:cmap(color_list[i]) for i in range(len(target_names))}
name_label = {target_name:target_name.split('-')[-1] for target_name in target_names}
#name_label = {target:target.split('-')[2] for target in target_names}



# set the list of figure frames (each frame is a matrix of varables to be plotted)
#fig_list = [ [[[target], ['Gene activity']]] for target in target_names ]
#fig_list = [[[[target], ['Gene activity']] for target in target_names]]
fig_list = [ [[[target for target in target_names], []]] ]
#fig_list = [ [[[target for target in target_names], ['Gene activity']]] ]
note_list = [target for target in target_names]


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

        if target_ftype == 'fract':
            target_fname = '_'.join([cell, sample, agent,
                                     str(int(target_binsize/1000.0)) + 'kb',
                                     'num']) + '.cn'

        else:
            target_fname = '_'.join([cell, sample, agent,
                                     str(int(target_binsize/1000.0)) + 'kb',
                                     target_ftype]) + '.cn'

        if target_ftype == 'anot':
           field_ID_value = load_file.read_tabular_file (path + target_fname, mode='col')
           ID_pos = field_ID_value['PhysicalPosition']

        elif target_ftype in ['score', 'zscore']:
            ID_pos, field_ID_value = read_bin_file(path2 + target_fname, chr_choices=[chr_choice])

        elif target_ftype in ['num']:
            ID_pos, field_ID_value = read_bin_num(path + target_fname, chr_choices=[chr_choice])

        elif target_ftype == 'fract':
            ID_pos, field_ID_value = read_bin_fract(path2 + target_fname, chr_choices=[chr_choice])      

        elif target_ftype in ['Chalf', 'zChalf']:
            ID_pos, ID_value = read_bin_Chalf(path + target_fname, chr_choices=[chr_choice])
            field_ID_value = {}
            field_ID_value['%s-%s-%s-%s' % (cell, sample, agent, tnum)] = ID_value

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

        # read RNA-seq file
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
    for fig_frame, note in zip(fig_list, note_list):
        # the case of single plot
        if len(fig_frame) == 1:

            fig = plt.figure(figsize=(15,5))
            ax1 = plt.subplot(211)
            ax2 = plt.subplot(212)
            ax1p = ax1.twinx()

            axes = [ax1, ax1p]
            left_names, right_names = fig_frame[0]
            var_names_list = [left_names, right_names]

            for i in range(len(axes)):
                ax = axes[i]
                var_names = var_names_list[i]
                alpha = alphas[i]
                lw = lws[i]
                linestyle = linestyles[i]
                ylabel = ylabels[i]
                ycolor = ycolors[i]
                ylims = ylims_list[i]

                for name in var_names:
                    try:
                        color = name_color[name]
                    except:
                        color = None
                    try:
                        label = name_label[name]
                    except:
                        label = None

                    ax.plot(name_sig[name],
                            color=color,
                            alpha=alpha,
                            lw=lw,
                            linestyle=linestyle,
                            label=label)

                    ax.set_ylabel(ylabel, color=ycolor)
                    ax.tick_params('y', colors=ycolor)
                    ax.set_ylim(ylims)

            if shaded:
                for gID in sorted(gID_binterval):
                    st, ed = gID_binterval[gID]
                    if gID_Gband[gID]['type'] == 'pos':
                        ax1.axvspan(st, ed-1, alpha=0.15, color='grey')

            ax1.set_xlim([0, len(gband_img)+1])
            ax1.set_xticks(xtick_locs[::10])
            ax1.set_xticklabels(xtick_labels[::10])
            ax1.set_xlabel("Position (Mb)")

            if name_label:
                leg = ax1.legend(framealpha=1, loc=legend_loc)
                for legobj in leg.legendHandles:
                    legobj.set_alpha(1.0)
                    legobj.set_linewidth(2.0)

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
            plt.savefig("_".join(['Gwide', chr_choice, note]) + ".svg",
                        format='svg',
                        bbox_inches='tight')
            #plt.show()
            plt.close()

        else:
            # the case of multiple target-feature pairs
            figwidth = round(6.4*float(len(gband_img))/10000, 1)
            figheight = 1 + 1.4*len(fig_frame)
            
            fig, axes = plt.subplots(nrows=1+len(fig_frame),
                                     ncols=1,
                                     figsize=(figwidth, figheight),
                                     sharex=True,
                                     gridspec_kw = {'hspace':0.12,
                                                    'wspace':0.12,
                                                    'height_ratios':[7]*len(fig_frame)+[1]})

            for k in range(len(fig_frame)):
                k_axes = [axes[k], axes[k].twinx()]
                left_names, right_names = fig_frame[k]
                var_names_list = [left_names, right_names]

                for i in range(len(k_axes)):
                    ax = k_axes[i]
                    var_names = var_names_list[i]
                    alpha = alphas[i]
                    lw = lws[i]
                    linestyle = linestyles[i]
                    ylabel = ylabels[i]
                    ycolor = ycolors[i]
                    ylims = ylims_list[i]

                    for name in var_names:
                        try:
                            color = name_color[name]
                        except:
                            color = None
                        try:
                            label = name_label[name]
                        except:
                            label = None

                        ax.plot(name_sig[name],
                                color=color,
                                alpha=alpha,
                                lw=lw,
                                linestyle=linestyle,
                                label=label)

                        ax.set_ylabel(ylabel, color=ycolor)
                        ax.tick_params('y', colors=ycolor)
                        ax.set_ylim(ylims)

                        #minv, maxv = ylims
                        #if minv !=None and maxv!=None:
                        #    ax.set_yticks([minv + 0.2*(maxv-minv),
                        #                   maxv - 0.2*(maxv-minv)])
                        #    ax.set_yticklabels([str(minv + 0.2*(maxv-minv)),
                        #                        str(maxv - 0.2*(maxv-minv))])

                        ax.spines['top'].set_visible(False)
                        ax.spines['bottom'].set_visible(False)
                        ax.spines['left'].set_visible(False)
                        ax.spines['right'].set_visible(False)
                        ax.tick_params(top='off',
                                       bottom='off',
                                       left='on',
                                       right='on',
                                       labelbottom='off')

                if shaded:
                    for gID in sorted(gID_binterval):
                        st, ed = gID_binterval[gID]
                        if gID_Gband[gID]['type'] == 'pos':
                            axes[k].axvspan(st, ed-1, alpha=0.15, color='grey')

                if name_label:
                    leg = axes[k].legend(framealpha=1, loc=legend_loc)
                    for legobj in leg.legendHandles:
                        legobj.set_alpha(1.0)
                        legobj.set_linewidth(2.0)

                axes[k].set_zorder(1)  # default zorder is 0 for ax1 and ax2
                axes[k].patch.set_visible(False)  # prevents ax1 from hiding ax2


            axes[len(fig_frame)].imshow(np.transpose(gband_img), cmap='Greys', aspect='auto')
            axes[len(fig_frame)].imshow(np.transpose(gband_cenimg), cmap ='Reds', vmin=0, vmax=20, aspect='auto')
            axes[len(fig_frame)].imshow(np.transpose(gband_varimg), cmap ='Purples', vmin=0, vmax=20, aspect='auto')
            axes[len(fig_frame)].set_xticks([])
            axes[len(fig_frame)].set_yticks([])

            plt.savefig("_".join(['Gwide', chr_choice, note]) + ".svg",
                        format='svg',
                        bbox_inches='tight')
            #plt.show()
            plt.close()
