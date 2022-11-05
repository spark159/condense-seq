import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
import pickle
from sklearn import linear_model
from sklearn.neighbors import LocalOutlierFactor
from sklearn.covariance import EllipticEnvelope
import seaborn as sns
from scipy import stats
from matplotlib.colors import LinearSegmentedColormap

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
chr_list = ['chr' + str(i) for i in range(1, 23)]
chr_list += ['chrX']
#chr_list += ['chrX', 'chrY']
#chr_list = ['chr9']
#chr_list.remove('chr13')

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
gID_FPKM1 = read_tsv(path+"ENCFF174OMR.tsv") # H1-hESC
gID_FPKM2 = read_tsv(path+"ENCFF345SHY.tsv") # GM12878
#gID_FPKM2 = load_file.read_RPKM("GSE63124_all_gene_raw_readcounts.txt", "Homo_sapiens.GRCh37.87.gtf", chr_list = chr_list)

# read gtf file
gID_field_values = load_file.read_GTF (path+"ENCFF159KBI.gtf", chr_list=chr_list, mode="gene")

# calculate score for all genes
gID_mscore1 = {}
gID_mscore2 = {}
for chr in chr_list:
    print "processing " + chr
    
    #fname1 = '%s_gene_score1' % (chr)
    #fname2 = '%s_gene_score2' % (chr)

    fname1 = '%s_gene_corrected_score1' % (chr)
    fname2 = '%s_gene_corrected_score2' % (chr)


    try:
        chr_gID_mscore1 = pickle.load(open(path+fname1 + ".pickle", "rb"))
        chr_gID_mscore2 = pickle.load(open(path+fname2 + ".pickle", "rb"))

    except:

        # read gene locations
        chr_gID_field_values1 = load_file.read_GTF (path+"ENCFF159KBI.gtf", chr_list=[chr], mode="gene")
        chr_gID_field_values2 = load_file.read_GTF (path+"ENCFF159KBI.gtf", chr_list=[chr], mode="gene")
        #chr_gID_field_values2 = load_file.read_GTF_old ("Homo_sapiens.GRCh37.87.gtf", chr_list=[chr], mode="gene")

        gIDs = list(set(chr_gID_field_values1.keys()) & set(chr_gID_field_values2.keys()))

        gID_ginterval1, gID_ginterval2 = {}, {}
        for gID in gIDs:
            try:
                TSS1 = chr_gID_field_values1[gID]['TSS']
                interval1 = (TSS1-2500, TSS1+2500)
            except:
                continue
            try:
                TSS2 = chr_gID_field_values2[gID]['TSS']
                interval2 = (TSS2-2500, TSS2+2500)
            except:
                continue
            gID_ginterval1[gID] = interval1
            gID_ginterval2[gID] = interval2

        del chr_gID_field_values1, chr_gID_field_values2

        # read NCP scores
        #anot_fname1 = "H1_NCP_sp_%s_anot.cn" % (chr)
        #nID_chr1, nID_pos1, name_nID_value1 = load_file.read_anot_file(path+anot_fname1)
        #nID_score1 = name_nID_value1['work/2021_06_07_H1_sp_detail/H1-NCP-sp-8']
        #del nID_chr1, name_nID_value1

        #anot_fname2 = "GM_NCP_sp_%s_anot.cn" % (chr)
        #nID_chr2, nID_pos2, name_nID_value2 = load_file.read_anot_file(path+anot_fname2)
        #nID_score2 = name_nID_value2['/home/spark159/scratch4-tha4/sangwoo/2022_09_08_GM_sp_H1_HP1a_deep/GM-NCP-sp-8']
        #del nID_chr2, name_nID_value2

        # read NCP scores (from corrected NCP number files)
        num_fname1 = "H1_NCP_sp_%s_num.cn" % (chr)
        #num_fname1 = "mCD8T_WT-NCP_sp_%s_num.cn" % (chr)
        nID_pos1, _, nID_score1 = read_num_file_new (path+num_fname1, chr_choice=chr)

        num_fname2 = "GM_NCP_sp_%s_num.cn" % (chr)
        #num_fname2 = "mCD8T_inht-NCP_sp_%s_num.cn" % (chr)
        nID_pos2, _, nID_score2 = read_num_file_new (path+num_fname2, chr_choice=chr)

        del _


        # get scores near TSS
        data_sets = [(gID_ginterval1, nID_pos1, nID_score1), (gID_ginterval2, nID_pos2, nID_score2)]
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

        gID_scores1, gID_scores2 = output_sets
        del ginterval_dict, data_sets, output_sets

        gIDs = list(set(gID_scores1.keys()) & set(gID_scores2.keys()))

        # get average scores near TSS
        chr_gID_mscore1, chr_gID_mscore2 = {}, {}
        for gID in gIDs:
            if len(gID_scores1[gID]) >= 10:
                chr_gID_mscore1[gID] = np.mean(gID_scores1[gID])
            if len(gID_scores2[gID]) >= 10:
                chr_gID_mscore2[gID] = np.mean(gID_scores2[gID])

        # save pickle file
        pickle.dump(chr_gID_mscore1, open(path+fname1 + ".pickle", "wb"))
        pickle.dump(chr_gID_mscore2, open(path+fname2 + ".pickle", "wb"))

    # update total data set
    gID_mscore1.update(chr_gID_mscore1)
    gID_mscore2.update(chr_gID_mscore2)

    print len(set(chr_gID_mscore1.keys()) & set(chr_gID_mscore2.keys()))

    print 

# finalize common gene set
gIDs = list(set(gID_mscore1) & set(gID_mscore2) & set(gID_FPKM1.keys()) & set(gID_FPKM2.keys()))
print 'Total gene count:' + str(len(gIDs))

# standardization of scores
scores1 = [gID_mscore1[gID] for gID in gIDs]
scores2 = [gID_mscore2[gID] for gID in gIDs]
score_mean1, score_std1 = np.mean(scores1), np.std(scores1)
score_mean2, score_std2 = np.mean(scores2), np.std(scores2)

for gID in gIDs:
    gID_mscore1[gID] = float(gID_mscore1[gID] - score_mean1) / score_std1
    gID_mscore2[gID] = float(gID_mscore2[gID] - score_mean2) / score_std2


# find Ensemble Gene ID for ESC marker genes / PC marker genes
ESC_gname_gIDs = {gname :[] for gname in ESC_gnames}
PC_gname_gIDs = {gname:[] for gname in PC_gnames}
HOX_gname_gID = {}
for gID in gIDs:
    gname = gID_field_values[gID]['geneName'].upper()
    if gname.startswith('HOX'):
        HOX_gname_gID[gname] = gID
    try:
        ESC_gname_gIDs[gname].append(gID)
        print gname
    except:
        pass
    try:
        PC_gname_gIDs[gname].append(gID)
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

fig = plt.figure(figsize=(2.8,2.4))
#fig = plt.figure()
#plt.plot(inX, inY, color='tab:blue', marker=',', linestyle='none', markersize=3, alpha=0.3)
plt.plot(inX, inY, color='tab:blue', marker=',', linestyle='none', markersize=35, alpha=0.5)
for x, y, gname in zip(outX, outY, gnames):
    #plt.plot(x, y, 'r.', markersize=3, alpha=1, zorder=10)
    #if gname in ESC_tf_cores:
    if gname in ESC_gnames:
        if gname in ESC_tf_cores:
            plt.plot(x, y, 'r.', markersize=3, alpha=1, zorder=10)
            plt.annotate(gname, (x, y), color='red', zorder=40, size=8, weight='bold')
        else:
            #plt.plot(x, y, 'k.', markersize=3, alpha=1, zorder=10)
            #plt.annotate(gname, (x, y), color='black', zorder=40, size=8)
            plt.plot(x, y, 'r.', markersize=3, alpha=1, zorder=10)
            
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
plt.savefig('nearTSScondvsExpreesion_hESC.png', bbox_inches='tight', dpi=300)
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
plt.scatter(X, Y, c=C, cmap=pastel_jet, vmin=-5, vmax=5, alpha=0.2, s=2)
        
for gID in ESC_gID_gname:
    gname = ESC_gID_gname[gID]
    x, y = gID_mscore1[gID], gID_mscore2[gID]
    plt.plot(x, y, 'kx', markersize=5, alpha=1, zorder=10, mew=1.5)
    if gname in ESC_tf_cores:
        plt.annotate(gname, (x, y), color='black', zorder=40, size=8, weight='bold')

#plt.plot([-2.5, 2], [-3.5, 3.5], 'k--', alpha=0.5)
#plt.plot([min(X), max(X)], [min(Y), max(Y)], 'k--', alpha=0.7)
#plt.plot([-3, 3], [-3, 3], 'k--', alpha=0.7)
plt.plot([-3, 3], [-3, 3], 'k--', alpha=0.7)
#plt.plot([0, 7], [0, 7], 'k--', alpha=0.7)
plt.title("Condensability near TSS (5kb)")
plt.xlabel('H1 hESC')
#plt.ylabel('A38-41 hPanc')
plt.ylabel('GM12878')
#plt.ylim([-2.5, 2.5])
#plt.xlim([-2.5, 2])
#plt.ylim([-3.5, 3.5])
plt.xlim([-3, 3])
plt.ylim([-3, 3])
#plt.xlim([0, 7])
#plt.ylim([0, 7])
cbar = plt.colorbar()
#cbar.ax.set_ylabel('Gene expression (A38-41 hPanc - H1 hESC)', rotation=-90, va="bottom")
cbar.ax.set_ylabel('Gene expression (GM12878 - H1 hESC)', rotation=-90, va="bottom")
#plt.savefig('A38VSH1hESC_all.png')
plt.savefig('GMVSH1hESC_all.png', bbox_inches='tight', dpi=300)
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
