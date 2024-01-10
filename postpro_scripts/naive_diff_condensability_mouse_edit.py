import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
import random
import pickle
from sklearn import linear_model
from sklearn.neighbors import LocalOutlierFactor
from sklearn.covariance import EllipticEnvelope
import seaborn as sns
from scipy import stats
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn
from sklearn.decomposition import PCA


# custom diverging colormap with white background
pastel_jet_div = LinearSegmentedColormap.from_list('white_viridis',
                                             [(0, 'darkblue'),
                                              (0.1, 'blue'),
                                              (0.2, 'tab:blue'),
                                              (0.4, 'tab:cyan'),
                                              (0.5, 'ivory'),
                                              (0.6, 'tab:orange'),
                                              (0.8, 'tab:red'),
                                              (0.9, 'red'),
                                              (1, 'darkred')
                                             ], N=256)

# "jet-like" colormap with white background
pastel_jet = LinearSegmentedColormap.from_list('white_viridis',
                                             [(0, '#ffffff'),
                                              (0.03, 'tab:cyan'),
                                              (0.1, 'tab:blue'),
                                              (0.3, 'tab:green'),
                                              (0.5, 'yellow'),
                                              (0.7, 'tab:orange'),
                                              (0.9, 'tab:red'),
                                              (1, 'darkred')
                                             ], N=256)

def draw_GO (fname, num, color='gold', alpha=1):
    goID_field_value = load_file.read_tabular_file (fname)

    sig_goID = sorted([(goID_field_value[goID]['P-value'], goID) for goID in goID_field_value])
    #sig_goID = sig_goID[:min(20, len(sig_goID))]
    sig_goID = sig_goID[:min(num, len(sig_goID))]

    xlabels, sigs = [], []
    for sig, goID in sig_goID:
        xlabels.append(goID_field_value[goID]['Description'])
        sigs.append(-np.log10(sig))

    #fig = plt.figure(figsize=(4,7))
    fig = plt.figure()
    #plt.barh(range(len(sigs))[::-1], sigs, color='gold')
    plt.barh(range(len(sigs))[::-1], sigs, color=color, alpha=alpha)
    for i in range(len(xlabels)):
        plt.annotate(xlabels[i], (0.5, len(xlabels)-1-i), ha='left', va='center', fontsize=8, weight='bold')

    plt.ylim([-1+0.25, len(xlabels)-0.25])
    plt.gca().tick_params(left='off')
    plt.yticks([])
    plt.xlabel('-log10(p)')
    #plt.savefig(fname.rsplit('.', 1)[0] + '.png', bbox_inches='tight')
    plt.savefig(fname.rsplit('.', 1)[0] + '.svg', format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()


def density_scatter(x , y, ax = None, sort = True, bins = 20, density = False, **kwargs )   :
    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d(x, y, bins = bins, density=density )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T ,
                 method = "splinef2d", bounds_error = False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        #print idx
        x, y, z = x[idx], y[idx], z[idx]

    img = ax.scatter( x, y, c=z, **kwargs )
    cbar = plt.colorbar(img)
    cbar.ax.tick_params(labelsize=5)
    #cbar = plt.colorbar(cm.ScalarMappable(norm = norm), ax=img)
    #cbar.ax.set_ylabel('Density')

    return ax



def tuple_cmp (a,b):
    if a[0] <= b[0]:
        return -1
    else:
        return 1


# read num file
def read_num_file_new (fname, chr_choice):
    ID_pos = {}
    ID_score1, ID_score2 = {}, {}
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        if cols[-1] == '*':
            #print "here"
            continue
        chr = cols[1]
        if cols[1] != chr_choice:
            continue
        #ID = chr + ':' + cols[0]
        ID = cols[0]
        pos = int(cols[2])
        sup1 = float(cols[3])
        sup2 = float(cols[4])
        control = float(cols[5])
        if sup1 * sup2 * control <= 0:
            continue
        score1 = -np.log(sup1/control)
        score2 = -np.log(sup2/control)
        ID_pos[ID] = pos
        ID_score1[ID] = score1
        ID_score2[ID] = score2
    return ID_pos, ID_score1, ID_score2

def get_corr(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = np.average(x)
    avg_y = np.average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff
    return diffprod / np.sqrt(xdiff2 * ydiff2)

# read raw count file and convert it to RPKM
def read_RPKM (fname, field_name, gtf_fname, chr_list=None):
    gID_field_values = load_file.read_GTF (gtf_fname, chr_list=chr_list, mode="gene")

    # get total exon length
    gID_exonlen = {}
    for gID in gID_field_values:
        length = 0
        for start, end in gID_field_values[gID]['exons']:
            length +=  end - start + 1
        gID_exonlen[gID] = length

    field_gID_counts = load_file.read_tabular_file (fname, mode="col")
    vgID_counts = field_gID_counts[field_name]

    gID_counts = {}
    total_counts = 0.0
    for vgID in vgID_counts:
        gID = vgID.split('.')[0] # strip off version
        counts = vgID_counts[vgID]
        gID_counts[gID] = counts # exclude the case of zero counts
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


def read_tsv (fname, chr_choices=None):
    First = True
    geneID_FPKM = {}
    #for fname in fnames:
    for line in open(fname):
        if First:
            First = False
            continue
        cols = line.strip().split()
        geneID, FPKM = cols[0], float(cols[6])
        geneID = geneID.split('.')[0]
        #assert geneID not in geneID_FPKM
        geneID_FPKM[geneID] = FPKM
    return geneID_FPKM

def read_bivalent (fname, chr_choices=None):
    gname_info = {}
    btype_gnames = {}
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        _, gname, chr, st, ed, strand, btype = cols[:7]
        if chr_choices and chr not in chr_choices:
            continue
        #assert gname not in gname_info
        if gname in gname_info:
            continue
        st, ed = int(st), int(ed)
        gname_info[gname] = {}
        gname_info[gname]['chr'] = chr
        gname_info[gname]['range'] = (st, ed)
        gname_info[gname]['strand'] = strand
        gname_info[gname]['btype'] = btype

        if btype not in btype_gnames:
            btype_gnames[btype] = []
        btype_gnames[btype].append(gname)
    return gname_info, btype_gnames

# set data path
path = "/home/spark159/../../media/spark159/sw/"

# chr list for analysis
#chr_list = ['chr' + str(i) for i in range(1, 23)]

chr_list = ['chr' + str(i) for i in range(1, 20)]
chr_list += ['chrX']
#chr_list += ['chrX', 'chrY']
#chr_list = ['chr9']
#chr_list.remove('chr13')

# mouse CD4 T cell lineage commitment genes
mCD4Tcell_gnames = ['Tbx21', 'Gata3', 'Rorc', 'Ifng', 'Il17']
#mCD4Tcell_gnames = [gname.upper() for gname in mCD4Tcell_gnames]

# stem cell marker genes
ESC_tf_cores =  ['Pou5f1', 'Sox2', 'KLF4', 'Nanog']
ESC_tf_others = ['Zfp42', 'UTF1', 'ZFX', 'TBN', 'FoxD3', 'HMGA2', 'NAC1', 'NR6A1', 'Stat3', 'LEF1', 'TCF3', 'Sall4', 'Fbxo15', 'L1TD1', 'Gdf3', 'Dppa5', 'Dppa4', 'Dppa2', 'Dppa3']
ESC_tf_cores = [gname.upper() for gname in ESC_tf_cores]
ESC_tf_others = [gname.upper() for gname in ESC_tf_others]
ESC_gnames = ESC_tf_cores + ESC_tf_others

# Pancreatic cancer cell marker genes
#PC_gnames = ['BRCA1', 'BRCA2', 'PALB2', 'CDKN2A', 'ATM', 'TP53', 'STK11', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'EPCAM']
PC_gnames = ['COL1A1', 'KRT17', 'ceacaM5', 'S100P', 'COL10A1', 'SerPinB5', 'GJB2', 'COL17A1', 'cXcl5', 'TMPRSS4', 'SDR16C5', 'CTHRC1', 'COL11A1', 'SLC6A14', 'MMP11', 'SULF1', 'Fn1', 'POSTN', 'ccl18', 'Muc4']
PC_gnames = [gname.upper() for gname in PC_gnames]

# read gene expression file
#gID_FPKM1 = read_tsv(path+"ENCFF174OMR.tsv") # H1-hESC
#gID_FPKM2 = read_tsv(path+"ENCFF345SHY.tsv") # GM12878
#gID_FPKM2 = load_file.read_RPKM("GSE63124_all_gene_raw_readcounts.txt", "Homo_sapiens.GRCh37.87.gtf", chr_list = chr_list)

# read gtf file
#gID_field_values = load_file.read_GTF (path+"ENCFF159KBI.gtf", chr_list=chr_list, mode="gene") #human
gID_field_values = load_file.read_GTF (path+"gencodeM21pri-UCSC-tRNAs-ERCC-phiX.gtf", mode="gene", chr_list=chr_list) #mouse
#gID_field_values = load_file.read_GTF (path+"gencodeM21pri-UCSC-tRNAs-ERCC-phiX.gtf", mode="gene", chr_list=chr_list, strip_ver=False) #mouse #no version strip

## read txt RNA_seq file and get geneID_FPKM (mouse CD8 T cell)
#gname_FPKM = read_txtRNA_seq("GSM3721901_W_effector_batch_1_1.txt")
#gID_FPKM = {}
#for gID in gID_field_values:
#    gname = gID_field_values[gID]["geneName"]
#    try:
#        FPKM = gname_FPKM[gname]
#        gID_FPKM[gID] = FPKM
#    except:
#        continue
#gID_FPKM1 = gID_FPKM
#gID_FPKM2 = gID_FPKM

# get old ODC data RPKM
gID_FPKM1 = read_RPKM (path + 'GSE157596_ODC_counts.tsv', 'TH2_WT2', path+"gencodeM21pri-UCSC-tRNAs-ERCC-phiX.gtf", chr_list=chr_list)

gID_FPKM2 = read_RPKM (path + 'GSE157596_ODC_counts.tsv', 'TH2_KO2', path+"gencodeM21pri-UCSC-tRNAs-ERCC-phiX.gtf", chr_list=chr_list)

    


# calculate score for all genes
gID_mscore1 = {}
gID_mscore2 = {}
gID_mscore3 = {}
for chr in chr_list:
    print "processing " + chr
    
    #fname1 = 'mCD8Tcell_%s_gene_score1' % (chr)
    #fname2 = 'mCD8Tcell_%s_gene_score2' % (chr)

    fname1 = 'mCD8Tcell_%s_corrected_gene_score1' % (chr)
    fname2 = 'mCD8Tcell_%s_corrected_gene_score2' % (chr)
    fname3 = 'mCD8Tcell_%s_corrected_gene_score3' % (chr)

    #fname1 = 'mCD8Tcell_%s_corrected_gene_score1_ver' % (chr)
    #fname2 = 'mCD8Tcell_%s_corrected_gene_score2_ver' % (chr)
    #fname3 = 'mCD8Tcell_%s_corrected_gene_score3_ver' % (chr)

    #fname1 = 'mCD8Tcell_%s_corrected_upstream_gene_score1' % (chr)
    #fname2 = 'mCD8Tcell_%s_corrected_upstream_gene_score2' % (chr)


    try:
        chr_gID_mscore1 = pickle.load(open(path+fname1 + ".pickle", "rb"))
        chr_gID_mscore2 = pickle.load(open(path+fname2 + ".pickle", "rb"))
        chr_gID_mscore3 = pickle.load(open(path+fname3 + ".pickle", "rb"))

    except:

        # read gene locations
        chr_gID_field_values = load_file.read_GTF (path+"gencodeM21pri-UCSC-tRNAs-ERCC-phiX.gtf", chr_list=[chr], mode="gene")
        #chr_gID_field_values = load_file.read_GTF (path+"gencodeM21pri-UCSC-tRNAs-ERCC-phiX.gtf", chr_list=[chr], mode="gene", strip_ver=False) # no version strip
        ##chr_gID_field_values2 = load_file.read_GTF_old ("Homo_sapiens.GRCh37.87.gtf", chr_list=[chr], mode="gene")

        gIDs = chr_gID_field_values.keys()

        print len(gIDs)

        gID_ginterval = {}
        for gID in gIDs:
            try:
                TSS = chr_gID_field_values[gID]['TSS']
                interval = (TSS-2500, TSS+2500)
            except:
                continue
            gID_ginterval[gID] = interval

        del chr_gID_field_values

        ## read NCP scores
        #anot_fname1 = "mCD8T_WT-NCP_sp_%s_anot.cn" % (chr)
        #nID_chr1, nID_pos1, name_nID_value1 = load_file.read_anot_file(path+anot_fname1)
        #nID_score1 = name_nID_value1['/home/spark159/scratch4-tha4/sangwoo/2022_09_23_mouseCD8_detail/mCD8T-WT-NCP-sp-8']
        #del nID_chr1, name_nID_value1

        #anot_fname2 = "mCD8T_inht-NCP_sp_%s_anot.cn" % (chr)
        #nID_chr2, nID_pos2, name_nID_value2 = load_file.read_anot_file(path+anot_fname2)
        #nID_score2 = name_nID_value2['/home/spark159/scratch4-tha4/sangwoo/2022_09_23_mouseCD8_detail/mCD8T-inht-NCP-sp-8']
        #del nID_chr2, name_nID_value2


        # read NCP scores (from corrected NCP number files)
        num_fname1 = "mCD8T_WT-NCP_sp_%s_num.cn" % (chr)
        nID_pos1, _, nID_score1 = read_num_file_new (path+num_fname1, chr_choice=chr)

        num_fname2 = "mCD8T_inht-NCP_sp_%s_num.cn" % (chr)
        nID_pos2, _, nID_score2 = read_num_file_new (path+num_fname2, chr_choice=chr)

        num_fname3 = "mCD8T_KO-NCP_sp_%s_num.cn" % (chr)
        nID_pos3, _, nID_score3 = read_num_file_new (path+num_fname3, chr_choice=chr)

        del _


        # get scores near TSS
        data_sets = [(gID_ginterval, nID_pos1, nID_score1), (gID_ginterval, nID_pos2, nID_score2), (gID_ginterval, nID_pos3, nID_score3)]
        output_sets = []
        for i in range(len(data_sets)):
            gID_ginterval, nID_pos, nID_score = data_sets[i]
            ginterval_dict = Interval_dict.double_hash(gID_ginterval, 100000, 250000000)
            gID_scores = {}
            for nID in nID_pos:
                pos = nID_pos[nID]
                score = nID_score[nID]
                for gID in ginterval_dict.find(pos):
                    if gID not in gID_scores:
                        gID_scores[gID] = []
                    gID_scores[gID].append(score)
            output_sets.append(gID_scores)

        gID_scores1, gID_scores2, gID_scores3 = output_sets
        del ginterval_dict, data_sets, output_sets

        gIDs = list(set(gID_scores1.keys()) & set(gID_scores2.keys()) & set(gID_scores3.keys()))

        # get average scores near TSS
        chr_gID_mscore1, chr_gID_mscore2, chr_gID_mscore3 = {}, {}, {}
        for gID in gIDs:
            if len(gID_scores1[gID]) >= 10:
                chr_gID_mscore1[gID] = np.mean(gID_scores1[gID])
            if len(gID_scores2[gID]) >= 10:
                chr_gID_mscore2[gID] = np.mean(gID_scores2[gID])
            if len(gID_scores3[gID]) >= 10:
                chr_gID_mscore3[gID] = np.mean(gID_scores3[gID])


        # save pickle file
        pickle.dump(chr_gID_mscore1, open(path+fname1 + ".pickle", "wb"))
        pickle.dump(chr_gID_mscore2, open(path+fname2 + ".pickle", "wb"))
        pickle.dump(chr_gID_mscore3, open(path+fname3 + ".pickle", "wb"))


    # update total data set
    gID_mscore1.update(chr_gID_mscore1)
    gID_mscore2.update(chr_gID_mscore2)
    gID_mscore3.update(chr_gID_mscore3)

    print len(set(chr_gID_mscore1.keys()) & set(chr_gID_mscore2.keys()) & set(chr_gID_mscore3.keys()))

    print 

# finalize common gene set
gIDs = list(set(gID_mscore1) & set(gID_mscore2) & set(gID_mscore3) & set(gID_FPKM1.keys()) & set(gID_FPKM2.keys()))
#gIDs = list(set(gID_mscore1) & set(gID_mscore2) & set(gID_mscore3))

print 'Total gene count:' + str(len(gIDs))


## standardization of scores
scores1 = [gID_mscore1[gID] for gID in gIDs]
scores2 = [gID_mscore2[gID] for gID in gIDs]
scores3 = [gID_mscore3[gID] for gID in gIDs]
score_mean1, score_std1 = np.mean(scores1), np.std(scores1)
score_mean2, score_std2 = np.mean(scores2), np.std(scores2)
score_mean3, score_std3 = np.mean(scores3), np.std(scores3)

for gID in gIDs:
    gID_mscore1[gID] = float(gID_mscore1[gID] - score_mean1) / score_std1
    gID_mscore2[gID] = float(gID_mscore2[gID] - score_mean2) / score_std2
    gID_mscore3[gID] = float(gID_mscore3[gID] - score_mean3) / score_std3

    #gID_mscore1[gID] = float(gID_mscore1[gID] - score_mean1)
    #gID_mscore2[gID] = float(gID_mscore2[gID] - score_mean2)


## map gname to gIDs (not unique)
gname_gIDs = {}
for gID in gIDs:
    gname = gID_field_values[gID]['geneName']
    if gname not in gname_gIDs:
        gname_gIDs[gname] = []
    gname_gIDs[gname].append(gID)

### averaging values by gname
gnames = gname_gIDs.keys()
gname_score1, gname_score2, gname_score3 = {}, {}, {}
gname_FPKM1, gname_FPKM2 = {}, {}
for gname in gnames:
    score1 = np.mean([gID_mscore1[gID] for gID in gname_gIDs[gname]])
    score2 = np.mean([gID_mscore2[gID] for gID in gname_gIDs[gname]])
    score3 = np.mean([gID_mscore3[gID] for gID in gname_gIDs[gname]])
    gname_score1[gname] = score1
    gname_score2[gname] = score2
    gname_score3[gname] = score3
    FPKM1 = np.mean([gID_FPKM1[gID] for gID in gname_gIDs[gname]])
    FPKM2 = np.mean([gID_FPKM2[gID] for gID in gname_gIDs[gname]])
    gname_FPKM1[gname] = FPKM1
    gname_FPKM2[gname] = FPKM2


# write all gnames as background
f = open("gname_bg.txt", 'w')
for gname in gnames:
    print >> f, gname
f.close()


# by gname
diff_gnames = set([])
up_gnames1, up_gnames2 = [], []
down_gnames1, down_gnames2 =[], []
gname_test_list = [gname_score2, gname_score3]
gname_control = gname_score1
test_names = ['inht', 'KO']

#for gID_test, tname in zip(gID_test_list, test_names):
for gname_test, tname in zip(gname_test_list, test_names):

    X, Y = [], []
    C = []

    for gname in gnames:
        X.append(gname_control[gname])
        Y.append(gname_test[gname])
        C.append(np.log2(1+gname_FPKM2[gname]) - np.log2(1+gname_FPKM1[gname]))
    C = stats.zscore(C)

    #for gID in gIDs:
    #    X.append(gID_control[gID])
    #    Y.append(gID_test[gID])
    #    C.append(np.log2(1+gID_FPKM2[gID]) - np.log2(1+gID_FPKM1[gID]))
    #    #C.append(1)
    ##C = stats.zscore(C)
    #mean_C = np.mean(C)
    #std_C = np.std(C)

    # draw all genes
    fig = plt.figure()
    #plt.scatter(X, Y, c=C, cmap='Spectral', vmin=-3, vmax=3, alpha=0.3, s=3)
    plt.scatter(X, Y, c=C, cmap=pastel_jet_div, vmin=-5, vmax=5, alpha=0.3, s=2)
    #plt.scatter(X, Y, c=C, cmap=pastel_jet_div, vmin=mean_C-std_C, vmax=mean_C+std_C, alpha=0.3, s=2)
    #plt.plot(X, Y, 'k.', alpha=0.2, markersize=1.5)

    #for gID in mCD4Tcell_gID_gname:
    #    gname = mCD4Tcell_gID_gname[gID]
    #    x, y = gID_control[gID], gID_test[gID]
    #    plt.plot(x, y, 'kx', markersize=5, alpha=1, zorder=10, mew=1.5)
    #    plt.annotate(gname, (x, y), color='black', zorder=40, size=8, weight='bold')

    for gname in mCD4Tcell_gnames:
        try:
            x, y = gname_control[gname], gname_test[gname]
        except:
            continue
        plt.plot(x, y, 'kx', markersize=5, alpha=1, zorder=10, mew=1.5)
        plt.annotate(gname, (x, y), color='black', zorder=40, size=8, weight='bold')


    #plt.plot([min(X), max(X)], [min(Y), max(Y)], 'k--', alpha=0.7)
    #plt.plot([-2.5, 2], [-3.5, 3.5], 'k--', alpha=0.5)
    #plt.plot([-4, 3.5], [-4, 3.5], 'r--', alpha=0.7)
    plt.plot([-5, 3.5], [-5, 3.5], 'k--', alpha=0.7)
    plt.title("Condensability near TSS (5kb)")
    plt.xlabel('Mouse CD8 T cell (WT)')
    #plt.ylabel('A38-41 hPanc')
    plt.ylabel('Mouse CD8 T cell ' + '(ODC %s)' % (tname))
    #plt.ylim([-2.5, 2.5])
    #plt.xlim([-3, 2.5])
    #plt.ylim([-3, 2.5])
    #plt.xlim([-4, 3.5])
    #plt.ylim([-4, 3.5])
    plt.xlim([-5, 3.5])
    plt.ylim([-5, 3.5])
    cbar = plt.colorbar()
    #cbar.ax.set_ylabel('Gene expression (A38-41 hPanc - H1 hESC)', rotation=-90, va="bottom")
    #cbar.ax.set_ylabel('Gene expression (GM12878 - H1 hESC)', rotation=-90, va="bottom")
    cbar.ax.set_ylabel('Gene expression (%s - WT)' % (tname), rotation=-90, va="bottom")
    #plt.savefig('A38VSH1hESC_all.png')
    plt.savefig('WTVSODC%s_all.png' % (tname))
    #plt.show()
    plt.close()


    # write up rank file for GSEA
    dscore_gname = []
    for gname in gnames:
        dscore = gname_test[gname] - gname_control[gname]
        dscore_gname.append((dscore, gname))
    dscore_gname = sorted(dscore_gname, reverse=True)
    
    f = open("%s-WT.rnk" % (tname), 'w')
    for dscore, gname in dscore_gname:
        print >> f, "%s\t%f" % (gname, dscore)
    f.close()

    
    # find up/down differential genes by rank (top 15%)
    sample_size = int(len(gnames)*0.15)

    up_gnames, down_gnames = [], []

    for dscore, gname in dscore_gname[:sample_size]:
        up_gnames.append(gname)

    for dscore, gname in dscore_gname[::-1][:sample_size]:
        down_gnames.append(gname)


    f = open("rank_up_" + tname + ".txt", 'w')
    for gname in up_gnames:
        #print >> f, gID_field_values[gID]['geneName']
        print >> f, gname
    f.close()

    f = open("rank_down_" + tname + ".txt", 'w')
    for gname in down_gnames:
        #print >> f, gID_field_values[gID]['geneName']
        print >> f, gname
    f.close()


    diff_gnames |= set(up_gnames)
    diff_gnames |= set(down_gnames)

    if tname == 'inht':
        up_gnames1 += up_gnames
        down_gnames1 += down_gnames
    else:
        up_gnames2 += up_gnames
        down_gnames2 += down_gnames


    # plot scatter with differential genes
    middleX, middleY = [], []
    for gname in list(set(gnames)-set(up_gnames + down_gnames)):
        middleX.append(gname_control[gname])
        middleY.append(gname_test[gname])

    upX, upY = [], []
    for gname in up_gnames:
        upX.append(gname_control[gname])
        upY.append(gname_test[gname])

    downX, downY = [], []
    for gname in down_gnames:
        downX.append(gname_control[gname])
        downY.append(gname_test[gname])

    fig = plt.figure()
    plt.plot(middleX, middleY, color='tab:gray', marker='.', alpha=0.15, markersize=1.5)
    plt.plot(upX, upY, 'b.', alpha=0.15, markersize=1.5, label='Up')
    plt.plot(downX, downY, 'r.', alpha=0.15, markersize=1.5, label='Down')
    plt.plot([-5, 3.5], [-5, 3.5], 'k--', alpha=0.7)
    plt.title("Condensability near TSS (5kb)")
    plt.xlabel('Mouse CD8 T cell (WT)')
    plt.ylabel('Mouse CD8 T cell (ODC %s)' % (tname))
    plt.xlim([-5, 3.5])
    plt.ylim([-5, 3.5])
    leg = plt.legend(loc='best', numpoints=1, prop={'size': 10})
    for lh in leg.legendHandles:
        lh._legmarker.set_markersize(10)
        lh._legmarker.set_alpha(1)
    plt.savefig('WTVSODC%s_diff.png' % (tname), dpi=500, bbox_inches='tight')
    #plt.show()
    plt.close()

    # plot bar graph of GO analysis
    #fname1 = "GO_top_raw.csv"
    #fname2 = "GO_top.csv"
    #fname3 = "GO_bott_raw.csv"
    #fname4 = "GO_bott.csv"

    #fname1 = "GO_%s_up.csv" % (tname)
    #fname2 = "GO_%s_down.csv" % (tname)

    fname1 = "new_GO_%s_up.csv" % (tname)
    fname2 = "new_GO_%s_down.csv" % (tname)

    draw_GO(fname1, 10, 'tab:green', alpha=0.7)
    draw_GO(fname2, 10, 'tab:orange', alpha=0.7)

# plot dscore VS dscore plot
X, Y = [], []
for gname in gnames:
    x = gname_score2[gname] - gname_score1[gname]
    y = gname_score3[gname] - gname_score1[gname]
    X.append(x)
    Y.append(y)

upX1, upY1 = [], []
for gname in up_gnames1:
    x = gname_score2[gname] - gname_score1[gname]
    y = gname_score3[gname] - gname_score1[gname]
    upX1.append(x)
    upY1.append(y)

upX2, upY2 = [], []
for gname in up_gnames2:
    x = gname_score2[gname] - gname_score1[gname]
    y = gname_score3[gname] - gname_score1[gname]
    upX2.append(x)
    upY2.append(y)

downX1, downY1 = [], []
for gname in down_gnames1:
    x = gname_score2[gname] - gname_score1[gname]
    y = gname_score3[gname] - gname_score1[gname]
    downX1.append(x)
    downY1.append(y)

downX2, downY2 = [], []
for gname in down_gnames2:
    x = gname_score2[gname] - gname_score1[gname]
    y = gname_score3[gname] - gname_score1[gname]
    downX2.append(x)
    downY2.append(y)

fig = plt.figure()
density_scatter(X, Y, bins = [20,20], s=3, cmap=pastel_jet, ax=plt.gca(), sort=False)
#plt.plot(X, Y, ',')
plt.plot([-1.5, 1.5], [-1.5, 1.5], 'k--', ms=2, alpha=0.7)
#plt.plot(upX1, upY1, '.', ms=1.5, label='UP1', alpha=0.2)
#plt.plot(upX2, upY2, '.', ms=1.5, label='UP2', alpha=0.2)
#plt.plot(downX1, downY1, '.', ms=1.5, label='DOWN1', alpha=0.2)
#plt.plot(downX2, downY2, '.', ms=1.5, label='DOWN1', alpha=0.2)
#plt.axvline(x=min(upX1), color='r', linestyle='--')
#plt.axvline(x=max(downX1), color='r', linestyle='--')
#plt.axhline(y=min(upY2), color='r', linestyle='--')
#plt.axhline(y=max(downY2), color='r', linestyle='--')
plt.xlim([-2,2])
plt.ylim([-2,2])
#leg = plt.legend(loc='best', numpoints=1, prop={'size': 6})
#for lh in leg.legendHandles:
#    lh._legmarker.set_markersize(15)
#    lh._legmarker.set_alpha(1)
plt.xlabel("$\Delta$zscore (+inht-WT)")
plt.ylabel("$\Delta$zscore (ODC KO-WT)")
plt.title("Condensability near TSS (5kb)")
plt.savefig("zdscoreVSzdscore.png", dpi=500, bbox_inches='tight')
#plt.show()
plt.close()


sys.exit(1)

# PCA projection
data = []
for gID in gIDs:
    state = [gID_mscore2[gID] - gID_mscore1[gID],
             gID_mscore3[gID] - gID_mscore1[gID]]
    data.append(state)
pca = PCA().fit(data)
Xr = pca.transform(data)

gID_PCA = {}
PCA_gID = []
for i in range(len(Xr)):
    gID = gIDs[i]
    gID_PCA[gID] = -Xr[i][0]
    PCA_gID.append((-Xr[i][0], gID))

PCA_gID = sorted(PCA_gID)
sample_size = int(len(gIDs)*0.15)
up_gIDs_pca = []
down_gIDs_pca = []
up_PCAs = []
down_PCAs = []

for PCA, gID in PCA_gID[:sample_size]:
    down_gIDs_pca.append(gID)
    down_PCAs.append(PCA)

for PCA, gID in PCA_gID[::-1][:sample_size]:
    up_gIDs_pca.append(gID)
    up_PCAs.append(PCA)

f = open("rank_up_pca.txt", 'w')
for gID in up_gIDs_pca:
    print >> f, gID_field_values[gID]['geneName']
f.close()

f = open("rank_down_pca.txt", 'w')
for gID in down_gIDs_pca:
    print >> f, gID_field_values[gID]['geneName']
f.close()

fig = plt.figure()
plt.hist(gID_PCA.values(), bins=500)
#plt.hist(up_PCAs, bins=500, label='up')
#plt.hist(down_PCAs, bins=500, label='down')
plt.axvline(x=min(up_PCAs), color='k', linestyle='--')
plt.axvline(x=max(down_PCAs), color='k', linestyle='--')
plt.xlim([-2.5, 2.5])
plt.xlabel("PCA 1st component")
plt.ylabel("Counts")
plt.title("PCA projection")
#plt.legend()
plt.savefig("PCAproject.png")
plt.close()

fname1 = "GO_pca_up.csv"
fname2 = "GO_pca_down.csv"

draw_GO(fname1, 10, 'tab:green', alpha=0.7)
draw_GO(fname2, 10, 'tab:orange', alpha=0.7)


# heat map
sample_size = int(len(gIDs)*0.005)

sample_gIDs = []
for PCA, gID in PCA_gID[:sample_size]:
    sample_gIDs.append(gID)

for PCA, gID in PCA_gID[::-1][:sample_size]:
    sample_gIDs.append(gID)

sample_data = []
for gID in sample_gIDs:
    state = [gID_mscore2[gID] - gID_mscore1[gID],
             gID_mscore3[gID] - gID_mscore1[gID]]
    sample_data.append(state)

fig = plt.figure(figsize=(3,10))
plt.imshow(sample_data, cmap='coolwarm', aspect='auto', vmin=-1.5, vmax=1.5)
#plt.show()
plt.close()

# Hierarchial clustering of differential genes
diff_gIDs = list(diff_gIDs)
sample_data = []
sample_gIDs = []
for gID in random.sample(diff_gIDs, 50):
    sample_gIDs.append(gID)
    state = [gID_mscore2[gID] - gID_mscore1[gID],
             gID_mscore3[gID] - gID_mscore1[gID]]
    sample_data.append(state)

hmap = sns.clustermap(sample_data,
                      method='average',
                      metric='euclidean',
                      figsize=(3,10),
                      cbar_kws=None,
                      row_cluster=True,
                      col_cluster=False,
                      dendrogram_ratio=0.2,
                      colors_ratio=0.03,
                      tree_kws=None,
                      cmap='coolwarm',
                      center=0,
                      xticklabels=['+ODCinht/WT', 'KO/WT'],
                      cbar_pos=None)

new_labels = []
for tick_label in hmap.ax_heatmap.axes.get_yticklabels():
    idx = int(tick_label.get_text())
    gID = sample_gIDs[idx]
    gname = gID_field_values[gID]["geneName"]
    tick_label.set_text(gname)
    new_labels.append(tick_label)

hmap.ax_heatmap.axes.set_yticklabels(new_labels, fontsize=5, rotation=0)
plt.gca().xaxis.tick_top()
plt.xticks(rotation=70)
plt.savefig("diffgID_hmap.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()


sys.exit(1)





diff_gIDs = list(diff_gIDs)
data = []
for gID in random.sample(diff_gIDs, 20):
    state = [gID_mscore1[gID], gID_mscore2[gID], gID_mscore3[gID]]
    data.append(state)

hmap = sns.clustermap(data,
                      method='average',
                      metric='euclidean',
                      figsize=(5,10),
                      cbar_kws=None,
                      row_cluster=True,
                      col_cluster=False,
                      dendrogram_ratio=0.2,
                      colors_ratio=0.03,
                      tree_kws=None,
                      cmap='Spectral_r',
                      center=0,
                      xticklabels=['WT', '+inht', 'KO'],
                      cbar_pos=None)

new_labels = []
for tick_label in hmap.ax_heatmap.axes.get_yticklabels():
    idx = int(tick_label.get_text())
    gID = diff_gIDs[idx]
    gname = gID_field_values[gID]["geneName"]
    tick_label.set_text(gname)
    new_labels.append(tick_label)

hmap.ax_heatmap.axes.set_yticklabels(new_labels, rotation=90)

plt.gca().xaxis.tick_top()
plt.xticks(rotation=70)
plt.savefig("diffgID_hmap.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()


 

sys.exit(1)

# find Ensemble Gene ID for ESC marker genes / PC marker genes
ESC_gname_gIDs = {gname :[] for gname in ESC_gnames}
PC_gname_gIDs = {gname:[] for gname in PC_gnames}
mCD4Tcell_gname_gIDs = {gname:[] for gname in mCD4Tcell_gnames}
HOX_gname_gID = {}

for gID in gIDs:
    gname = gID_field_values[gID]['geneName'].upper()
    if gname.startswith('HOX'):
        HOX_gname_gID[gname] = gID
    try:
        ESC_gname_gIDs[gname].append(gID)
    except:
        pass
    try:
        PC_gname_gIDs[gname].append(gID)
    except:
        pass
    try:
        mCD4Tcell_gname_gIDs[gname].append(gID)
    except:
        pass
    
ESC_gID_gname = {}
for gname in ESC_gname_gIDs:
    if len(ESC_gname_gIDs[gname]) == 1:
        gID = ESC_gname_gIDs[gname][0]
        assert gID not in ESC_gID_gname
        ESC_gID_gname[gID] = gname

PC_gID_gname = {}
for gname in PC_gname_gIDs:
    if len(PC_gname_gIDs[gname]) == 1:
        gID = PC_gname_gIDs[gname][0]
        assert gID not in PC_gID_gname
        PC_gID_gname[gID] = gname

mCD4Tcell_gID_gname = {}
for gname in mCD4Tcell_gname_gIDs:
    if len(mCD4Tcell_gname_gIDs[gname]) == 1:
        gID = mCD4Tcell_gname_gIDs[gname][0]
        assert gID not in mCD4Tcell_gID_gname
        mCD4Tcell_gID_gname[gID] = gname


HOX_gID_gname = {}
for gname in HOX_gname_gID:
    gID = HOX_gname_gID[gname]
    HOX_gID_gname[gID] = gname

# find bivalent genes
gname_binfo, btype_gnames = read_bivalent("Bivalent_info.csv")
btype_gID_gname = {}
for gID in gIDs:
    gname = gID_field_values[gID]['geneName']
    try:
        btype = gname_binfo[gname]['btype']
        if btype not in btype_gID_gname:
            btype_gID_gname[btype] = {}
        assert gID not in btype_gID_gname[btype]
        btype_gID_gname[btype][gID] = gname
    except:
        continue


# set in-IDs and out-IDs
in_gIDs = list(set(gIDs) - set(ESC_gID_gname.keys())) # all genes except ESC marker
out_gIDs = ESC_gID_gname.keys() # ESC marker genes


# H1 FPKM vs GM FPKM
X, Y = [], []
for gID in gIDs:
    X.append(np.log2(1+gID_FPKM1[gID]))
    Y.append(np.log2(1+gID_FPKM2[gID]))
fig = plt.figure()
plt.plot(X, Y, '.')
plt.title('RNA-seq data comparison')
plt.xlabel('H1 hESC logFPKM')
#plt.ylabel('A38-41 hPanc logFPKM')
plt.ylabel('GM12878 logFPKM')
#plt.show()
plt.close()

# H1 FPKM vs H1 score
inX = [np.log2(1+gID_FPKM1[gID]) for gID in in_gIDs]
inY = [gID_mscore1[gID] for gID in in_gIDs]
outX = [np.log2(1+gID_FPKM1[gID]) for gID in out_gIDs]
outY = [gID_mscore1[gID] for gID in out_gIDs]
gnames = [ESC_gID_gname[gID] for gID in out_gIDs]

# get correaltion
X, Y = inX + outX, inY + outY
print 'Spearman corr: ', statis.get_spearman_corr(X, Y)

# polynomial fitting
#feature_list = [[x, x**2, x**3] for x in X]
#test_list = [[y] for y in Y]
#reg = linear_model.Ridge(alpha=0.5)
#reg.fit(feature_list, test_list)
#Xrange = np.linspace(min(X), max(X), num=100)
#Ypred = reg.predict([[x, x**2, x**3] for x in Xrange])
#Ypred = [value[0] for value in Ypred]

group_gIDs = statis.partition({gID:np.log2(1+gID_FPKM1[gID]) for gID in gIDs}, 500)
meanX, meanY = [], []
for i in range(len(group_gIDs)):
    meanX.append(np.median([np.log2(1+gID_FPKM1[gID]) for gID in group_gIDs[i]]))
    meanY.append(np.median([gID_mscore1[gID] for gID in group_gIDs[i]]))

#fig = plt.figure(figsize=(2.8,2.4))
fig = plt.figure()
plt.plot(inX, inY, color='tab:blue', marker=',', linestyle='none', markersize=3, alpha=0.3)
for x, y, gname in zip(outX, outY, gnames):
    #plt.plot(x, y, 'r.', markersize=3, alpha=1, zorder=10)
    #if gname in ESC_tf_cores:
    if gname in ESC_gnames:
        if gname in ESC_tf_cores:
            plt.plot(x, y, 'r.', markersize=3, alpha=1, zorder=10)
            plt.annotate(gname, (x, y), color='red', zorder=40, size=8, weight='bold')
        else:
            plt.plot(x, y, 'k.', markersize=3, alpha=1, zorder=10)
            plt.annotate(gname, (x, y), color='black', zorder=40, size=8)
#sns.kdeplot(data=[[X[i], Y[i]] for i in range(len(X))])
#plt.plot(Xrange, Ypred, 'k--', alpha=0.7)
plt.plot(meanX, meanY, 'k--', alpha=0.7)
plt.title('H1-hESC (Near TSS)', fontsize=12)
plt.xlabel('Gene expression (logFPKM)', fontsize=12)
plt.ylabel('Condensabiltiy (A.U.)', fontsize=12)
plt.gca().tick_params(axis='both', which='major', labelsize=8)
plt.gca().tick_params(axis='both', which='minor', labelsize=8)
plt.ylim([-2.5, 2.5])
#plt.savefig('nearTSScondvsExpreesion_hESC.svg', format='svg', bbox_inches='tight')
plt.savefig('nearTSScondvsExpreesion_hESC.png', bbox_inches='tight')
#plt.show()
plt.close()

# Panc FPKM vs GM score
inX = [np.log2(1+gID_FPKM2[gID]) for gID in in_gIDs]
inY = [gID_mscore2[gID] for gID in in_gIDs]
outX = [np.log2(1+gID_FPKM2[gID]) for gID in out_gIDs]
outY = [gID_mscore2[gID] for gID in out_gIDs]
gnames = [ESC_gID_gname[gID] for gID in out_gIDs]

X, Y = inX + outX, inY + outY
print 'Spearman corr: ', statis.get_spearman_corr(X, Y)

# polynomial fitting
#feature_list = [[x, x**2, x**3] for x in X]
#test_list = [[y] for y in Y]
#reg = linear_model.Ridge(alpha=0.5)
#reg.fit(feature_list, test_list)
#Xrange = np.linspace(min(X), max(X), num=100)
#Ypred = reg.predict([[x, x**2, x**3] for x in Xrange])
#Ypred = [value[0] for value in Ypred]

group_gIDs = statis.partition({gID:np.log2(1+gID_FPKM2[gID]) for gID in gIDs}, 500)
meanX, meanY = [], []
for i in range(len(group_gIDs)):
    meanX.append(np.median([np.log2(1+gID_FPKM2[gID]) for gID in group_gIDs[i]]))
    meanY.append(np.median([gID_mscore2[gID] for gID in group_gIDs[i]]))

fig = plt.figure()
plt.plot(inX, inY, color='tab:blue', marker=',', linestyle='none', markersize=3, alpha=0.3)
for x, y, gname in zip(outX, outY, gnames):
    #plt.plot(x, y, 'r.', markersize=3, alpha=1, zorder=10)
    #if gname in ESC_tf_cores:
    if gname in ESC_gnames:
        if gname in ESC_tf_cores:
            plt.plot(x, y, 'r.', markersize=3, alpha=1, zorder=10)
            plt.annotate(gname, (x, y), color='red', zorder=40, size=8, weight='bold')
        else:
            plt.plot(x, y, 'k.', markersize=3, alpha=1, zorder=10)
            plt.annotate(gname, (x, y), color='black', zorder=40, size=8)
#plt.plot(Xrange, Ypred, 'k--', alpha=0.7)
plt.plot(meanX, meanY, 'k--', alpha=0.7)
#plt.title('A38-41 hPanc')
#plt.title('A38-41 hPanc (Near TSS)', fontsize=12)
plt.title('GM12878 (Near TSS)', fontsize=12)
plt.xlabel('Gene expression (logFPKM)', fontsize=12)
plt.ylabel('Condensabiltiy (A.U.)', fontsize=12)
plt.gca().tick_params(axis='both', which='major', labelsize=8)
plt.gca().tick_params(axis='both', which='minor', labelsize=8)
plt.ylim([-3.5, 3.5])
#plt.savefig('nearTSScondvsExpreesion_hPanc.svg', format='svg', bbox_inches='tight')
#plt.savefig('nearTSScondvsExpreesion_GM.svg', format='svg', bbox_inches='tight')
plt.savefig('nearTSScondvsExpreesion_GM.png', bbox_inches='tight')
#plt.show()
plt.close()


## H1 score VS GM score
X, Y = [], []
C = []
for gID in gIDs:
    X.append(gID_mscore1[gID])
    Y.append(gID_mscore2[gID])
    C.append(np.log2(1+gID_FPKM2[gID]) - np.log2(1+gID_FPKM1[gID]))
C = stats.zscore(C)

outX = [gID_mscore1[gID] for gID in out_gIDs]
outY = [gID_mscore2[gID] for gID in out_gIDs]
gnames = [ESC_gID_gname[gID] for gID in out_gIDs]


# draw all genes

# custom diverging colormap with white background
pastel_jet = LinearSegmentedColormap.from_list('white_viridis',
                                             [(0, 'darkblue'),
                                              (0.1, 'blue'),
                                              (0.2, 'tab:blue'),
                                              (0.4, 'tab:cyan'),
                                              (0.5, 'ivory'),
                                              (0.6, 'tab:orange'),
                                              (0.8, 'tab:red'),
                                              (0.9, 'red'),
                                              (1, 'darkred')
                                             ], N=256)


fig = plt.figure()
#plt.scatter(X, Y, c=C, cmap='Spectral', vmin=-3, vmax=3, alpha=0.3, s=3)
plt.scatter(X, Y, c=C, cmap=pastel_jet, vmin=-5, vmax=5, alpha=0.3, s=2)
#plt.plot(X, Y, 'k.', alpha=0.2, markersize=1.5)
        
for gID in mCD4Tcell_gID_gname:
    gname = mCD4Tcell_gID_gname[gID]
    x, y = gID_mscore1[gID], gID_mscore2[gID]
    plt.plot(x, y, 'kx', markersize=5, alpha=1, zorder=10, mew=1.5)
    plt.annotate(gname, (x, y), color='black', zorder=40, size=8, weight='bold')

#plt.plot([min(X), max(X)], [min(Y), max(Y)], 'k--', alpha=0.7)
#plt.plot([-2.5, 2], [-3.5, 3.5], 'k--', alpha=0.5)
#plt.plot([-4, 3.5], [-4, 3.5], 'r--', alpha=0.7)
plt.plot([-5, 3.5], [-5, 3.5], 'k--', alpha=0.7)
plt.title("Condensability near TSS (5kb)")
plt.xlabel('Mouse CD8 T cell (WT)')
#plt.ylabel('A38-41 hPanc')
plt.ylabel('Mouse CD8 T cell (+ODC inhibitor)')
#plt.ylim([-2.5, 2.5])
#plt.xlim([-3, 2.5])
#plt.ylim([-3, 2.5])
#plt.xlim([-4, 3.5])
#plt.ylim([-4, 3.5])
plt.xlim([-5, 3.5])
plt.ylim([-5, 3.5])
cbar = plt.colorbar()
#cbar.ax.set_ylabel('Gene expression (A38-41 hPanc - H1 hESC)', rotation=-90, va="bottom")
#cbar.ax.set_ylabel('Gene expression (GM12878 - H1 hESC)', rotation=-90, va="bottom")
cbar.ax.set_ylabel('Gene expression (+inht - WT)', rotation=-90, va="bottom")
#plt.savefig('A38VSH1hESC_all.png')
plt.savefig('WTVSODCinht_all.png')
#plt.show()
plt.close()


# plot scatter with differential genes
middleX, middleY = [], []
for gID in list(set(gIDs)-set(up_gIDs + down_gIDs)):
    middleX.append(gID_mscore1[gID])
    middleY.append(gID_mscore2[gID])

upX, upY = [], []
for gID in up_gIDs:
    upX.append(gID_mscore1[gID])
    upY.append(gID_mscore2[gID])

downX, downY = [], []
for gID in down_gIDs:
    downX.append(gID_mscore1[gID])
    downY.append(gID_mscore2[gID])

fig = plt.figure()
plt.plot(middleX, middleY, color='tab:gray', marker='.', alpha=0.15, markersize=1.5)
plt.plot(upX, upY, 'b.', alpha=0.15, markersize=1.5, label='Up')
plt.plot(downX, downY, 'r.', alpha=0.15, markersize=1.5, label='Down')
plt.plot([-5, 3.5], [-5, 3.5], 'k--', alpha=0.7)
plt.title("Condensability near TSS (5kb)")
plt.xlabel('Mouse CD8 T cell (WT)')
plt.ylabel('Mouse CD8 T cell (+ODC inhibitor)')
plt.xlim([-5, 3.5])
plt.ylim([-5, 3.5])
leg = plt.legend(loc='best', numpoints=1, prop={'size': 6})
for lh in leg.legendHandles:
    lh._legmarker.set_markersize(15)
    lh._legmarker.set_alpha(1)
plt.savefig('WTVSODCinht_diff.png')
#plt.show()
plt.close()



    

# score difference vs mean gene expression (MA-plot)
X, Y = [] ,[]
gID_dscore = {}
for gID in gIDs:
    #mean_FPKM = 0.5*(np.log(1+gID_FPKM1[gID]) + np.log(1+gID_FPKM2[gID]))
    #mean_FPKM = 0.5*(gID_FPKM1[gID] + gID_FPKM2[gID])
    mean_score = 0.5*(gID_mscore1[gID] + gID_mscore2[gID])
    score_diff = gID_mscore2[gID] - gID_mscore1[gID]
    X.append(mean_score)
    Y.append(score_diff)
    gID_dscore[gID] = score_diff

fig = plt.figure()
#plt.plot(X, Y, 'k.', markersize=2, alpha=0.2)
plt.scatter(X, Y, c=C, cmap=pastel_jet, vmin=-5, vmax=5, alpha=0.3, s=2)
plt.axhline(y=0, linestyle='--', color='r')
plt.title("Condensability mean vs difference")
plt.xlabel('mean condensability')
#plt.ylabel('A38-41 hPanc - H1 hESC')
plt.ylabel('+ODC inht - WT')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Gene expression (+inht - WT)', rotation=-90, va="bottom")
plt.savefig('MA_plot_all.png')
#plt.show()
plt.close()



# draw only bivalent genes
for btype in btype_gID_gname:
    Bt_gID_gname = btype_gID_gname[btype]
    
    fig = plt.figure()
    plt.scatter(X, Y, c='lightgrey', alpha=0.05, s=1)

    Bt_X, Bt_Y = [], []
    Bt_C = []
    for gID in btype_gID_gname[btype]:
        gname = btype_gID_gname[btype][gID]
        x, y = gID_mscore1[gID], gID_mscore2[gID]
        Bt_X.append(x)
        Bt_Y.append(y)
        Bt_C.append(np.log2(1+gID_FPKM2[gID]) - np.log2(1+gID_FPKM1[gID]))
    plt.scatter(Bt_X, Bt_Y, c=Bt_C, s=3, alpha=0.5, zorder=10, vmin=-3, vmax=3,  cmap='bwr')

    #plt.plot([min(X), max(X)], [min(Y), max(Y)], 'k--', alpha=0.7)
    plt.plot([-2.5, 2], [-3.5, 3.5], 'k--', alpha=0.5)
    plt.title("Condensability near TSS (5kb) %s" % (btype))
    plt.xlabel('H1 hESC')
    plt.ylabel('GM12878')
    plt.xlim([-2.5, 2])
    plt.ylim([-3.5, 3.5])
    #plt.xlim([-2.5, 2.5])
    #plt.ylim([-2.5, 2.5])
    cbar = plt.colorbar()
    #cbar.ax.set_ylabel('Gene expression (A38-41 hPanc - H1 hESC)', rotation=-90, va="bottom")
    cbar.ax.set_ylabel('Gene expression (GM12878 - H1 hESC)', rotation=-90, va="bottom")
    #plt.savefig('A38VSH1hESC_%s.png' % (btype))
    plt.savefig('GMVSH1hESC_%s.png' % (btype))
    #plt.show()
    plt.close()


# draw only HOX genes
fig = plt.figure()
plt.scatter(X, Y, c='lightgrey', alpha=0.05, s=1)

HOX_X, HOX_Y = [], []
HOX_C = []
for gID in HOX_gID_gname:
    gname = HOX_gID_gname[gID]
    x, y = gID_mscore1[gID], gID_mscore2[gID]
    HOX_X.append(x)
    HOX_Y.append(y)
    HOX_C.append(np.log2(1+gID_FPKM2[gID]) - np.log2(1+gID_FPKM1[gID]))
plt.scatter(HOX_X, HOX_Y, c=HOX_C, s=12, alpha=1, zorder=10, vmin=-3, vmax=3, edgecolor='k', linewidth=0.5, cmap='bwr')

#plt.plot([min(X), max(X)], [min(Y), max(Y)], 'k--', alpha=0.7)
plt.plot([-2.5, 2], [-3.5, 3.5], 'k--', alpha=0.5)
plt.title("Condensability near TSS (5kb) only HOX genes")
plt.xlabel('H1 hESC')
#plt.ylabel('A38-41 hPanc')
plt.ylabel('GM12878')
#plt.xlim([-2.5, 2.5])
#plt.ylim([-2.5, 2.5])
plt.xlim([-2.5, 2])
plt.ylim([-3.5, 3.5])
cbar = plt.colorbar()
#cbar.ax.set_ylabel('Gene expression (A38-41 hPanc - H1 hESC)', rotation=-90, va="bottom")
#plt.savefig('A38VSH1hESC_HOX.png')
cbar.ax.set_ylabel('Gene expression (GM12878 - H1 hESC)', rotation=-90, va="bottom")
plt.savefig('GMVSH1hESC_HOX.png')
#plt.show()
plt.close()


# draw only ES markers
fig = plt.figure()
plt.scatter(X, Y, c='lightgrey', alpha=0.05, s=1)

ESC_X, ESC_Y = [], []
ESC_C = []
for gID in ESC_gID_gname:
    gname = ESC_gID_gname[gID]
    x, y = gID_mscore1[gID], gID_mscore2[gID]
    ESC_X.append(x)
    ESC_Y.append(y)
    ESC_C.append(np.log2(1+gID_FPKM2[gID]) - np.log2(1+gID_FPKM1[gID]))
plt.scatter(ESC_X, ESC_Y, c=ESC_C, s=12, alpha=1, zorder=10, vmin=-3, vmax=3, edgecolor='k', linewidth=0.5, cmap='bwr')

#plt.plot([min(X), max(X)], [min(Y), max(Y)], 'k--', alpha=0.7)
plt.plot([-2.5, 2], [-3.5, 3.5], 'k--', alpha=0.5)
plt.title("Condensability near TSS (5kb) only ES markers")
plt.xlabel('H1 hESC')
#plt.ylabel('A38-41 hPanc')
plt.ylabel('GM12878')
#plt.xlim([-2.5, 2.5])
#plt.ylim([-2.5, 2.5])
plt.xlim([-2.5, 2])
plt.ylim([-3.5, 3.5])
cbar = plt.colorbar()
#cbar.ax.set_ylabel('Gene expression (A38-41 hPanc - H1 hESC)', rotation=-90, va="bottom")
#plt.savefig('A38VSH1hESC_ESm.png')
cbar.ax.set_ylabel('Gene expression (GM12878 - H1 hESC)', rotation=-90, va="bottom")
plt.savefig('GMVSH1hESC_ESm.png')
#plt.show()
plt.close()





### H1 score VS Panc score
#inX = [gID_mscore1[gID] for gID in in_gIDs]
#inY = [gID_mscore2[gID] for gID in in_gIDs]
#outX = [gID_mscore1[gID] for gID in out_gIDs]
#outY = [gID_mscore2[gID] for gID in out_gIDs]
#gnames = [ESC_gID_gname[gID] for gID in out_gIDs]
#
#X = inX + outX
#Y = inY + outY
#
#fig = plt.figure()
#plt.plot(inX, inY, color='tab:blue', marker=',', linestyle='none', markersize=3, alpha=0.3)
#for x, y, gname in zip(outX, outY, gnames):
#    plt.plot(x, y, 'r.', markersize=3, alpha=1, zorder=10)
#    if gname in ESC_tf_cores:
#        plt.annotate(gname, (x, y), color='red', zorder=40, size=8, weight='bold')
#plt.plot([min(X), max(X)], [min(Y), max(Y)], 'k--', alpha=0.7)
#plt.title("Condensability near TSS (5kb)")
#plt.xlabel('H1 hESC')
#plt.ylabel('A38-41 hPanc')
#plt.xlim([-2.5, 2.5])
#plt.ylim([-2.5, 2.5])
#plt.show()
#plt.close()

#sys.exit(1)

#data_list = []
#for gID in gIDs:
#    data = [gID_mscore1[gID], gID_mscore2[gID]]
#    data_list.append(data)

#clf = EllipticEnvelope()
#outcheck = clf.fit_predict(data_list)


#fig = plt.figure()
#for i in range(len(data_list)):
#    data = data_list[i]
#    if outcheck[i] < 0:
#        #mark = 'r.'
#        mark = 'k.'
#    else:
#        mark = 'k.'
#    plt.plot([data[0]], [data[1]], mark, markersize=3, alpha=0.5)

#plt.title("Condensability near TSS (5kb)")
#plt.xlabel('H1 hESC')
#plt.ylabel('A38-41 hPanc')
#plt.show()
#plt.close()

#sys.exit(1)
    

#inX, inY = [], []
#outX, outY = [], []
#for i in range(len(X), len(Y)):
#    if outcheck[i] < 0:
#        outX.append(X[i])
#        outY.append(Y[i])
#    else:
#        inX.append(X[i])
#        inY.append(Y[i])

#reg = linear_model.Ridge(alpha=0.5)
#reg.fit (feature_list, test_list)
#Ypred = reg.predict(feature_list)
#Ypred = [ value[0] for value in Ypred]


#fig = plt.figure()
#plt.plot(X, Y, 'k.', markersize=3, alpha=0.5)
#plt.plot(X, Ypred, 'r--')
#plt.title("Condensability near TSS (5kb)")
#plt.xlabel('H1 hESC')
#plt.ylabel('A38-41 hPanc')
#plt.show()
#plt.close()

# gene expression difference VS score differences
inX = [np.log2(1+gID_FPKM2[gID]) - np.log2(1+gID_FPKM1[gID]) for gID in in_gIDs]
inY = [gID_mscore2[gID] - gID_mscore1[gID] for gID in in_gIDs]
outX = [np.log2(1+gID_FPKM2[gID]) - np.log2(1+gID_FPKM1[gID]) for gID in out_gIDs]
outY = [gID_mscore2[gID] - gID_mscore1[gID] for gID in out_gIDs]
gnames = [ESC_gID_gname[gID] for gID in out_gIDs]

fig = plt.figure()
plt.plot(inX, inY, color='tab:blue', marker='.', linestyle='none', markersize=2, alpha=0.2)
for x, y, gname in zip(outX, outY, gnames):
    plt.plot(x, y, 'r.', markersize=2, alpha=1, zorder=10)
    if gname in ESC_tf_cores:
        plt.annotate(gname, (x, y), color='red', zorder=40, size=8, weight='bold')
#plt.xlabel('logFPKM (A38-41 hPanc - H1 hESC)')
#plt.ylabel('Condensability (A38-41 hPanc - H1 hESC)')
plt.xlabel('logFPKM (GM12878 - H1 hESC)')
plt.ylabel('Condensability (GM12878 - H1 hESC)')
plt.axvline(x=0, color='k', linestyle='--', alpha=0.7)
plt.axhline(y=0, color='k', linestyle='--', alpha=0.7)
#plt.xlim([-2.5, 2.5])
#plt.ylim([-2.5, 2.5])
#plt.show()
plt.close()

sys.exit(1)

# score difference vs mean gene expression (MA-plot)
X, Y = [] ,[]
gID_dscore = {}
for gID in gIDs:
    #mean_FPKM = 0.5*(np.log(1+gID_FPKM1[gID]) + np.log(1+gID_FPKM2[gID]))
    mean_FPKM = 0.5*(gID_FPKM1[gID] + gID_FPKM2[gID])
    score_diff = gID_mscore2[gID] - gID_mscore1[gID]
    X.append(np.log2(1+mean_FPKM))
    Y.append(score_diff)
    gID_dscore[gID] = score_diff

fig = plt.figure()
plt.plot(X, Y, 'k.', markersize=3, alpha=0.3)
plt.axhline(y=0, linestyle='--', color='r')
plt.title("Condensability difference vs mean gene expression")
plt.xlabel('mean log2(1+FPKM))')
#plt.ylabel('A38-41 hPanc - H1 hESC')
plt.ylabel('GM12878 - H1 hESC')
plt.show()
plt.close()


# score difference histogram and partitions
# Partition by score
med = np.median(gID_dscore.values())
std = np.std(gID_dscore.values())
lines = [med-0.5*std-i*std for i in range(3)] + [med+0.5*std+i*std for i in range(3)]
lines = sorted(lines)
p_num = len(lines)+1
p_range = []
for i in range(p_num):
    if i == 0:
        st = -np.inf
        ed = lines[i]
    elif i == p_num-1:
        st = ed
        ed = np.inf
    else:
        st = ed
        ed = lines[i]
    p_range.append((st, ed))
fig = plt.figure()
plt.hist(gID_dscore.values(), bins=100)
for line in lines:
    plt.axvline(x=line, color='k', linestyle='--')
num_rom = {1:'I', 2:'II', 3:'III', 4:'IV', 5:'V', 6:'VI', 7:'VII'}
for i in range(p_num):
    st, ed = p_range[i]
    if i == 0:
        x = np.mean([-3, ed])
    elif i == p_num - 1:
        x = np.mean([st, 3])
    else:
        x = np.mean([st, ed])
    plt.text(x, 10000, num_rom[i+1], fontsize=20, va='center', ha='center')
plt.xlim([-3,3])
#plt.title("Chromosome1")
plt.title("All genes")
#plt.xlabel("Condensability difference (A38-41 hPanc - H1 hESC)")
plt.xlabel("Condensability difference (GM12878 - H1 hESC)")
plt.ylabel("Gene Counts")
#plt.savefig("partition_hist.png")
plt.show()
plt.close()

p_IDs = [[] for i in range(p_num)]
for ID in gID_dscore:
    dscore = gID_dscore[ID]
    for i in range(p_num):
        st, ed = p_range[i]
        if dscore >= st and dscore < ed:
            break
    p_IDs[i].append(ID)

for i in range(len(p_IDs)):
    f = open("output_" + str(i) + ".txt", "w")
    for ID in p_IDs[i]:
        print >> f, ID
    f.close()
