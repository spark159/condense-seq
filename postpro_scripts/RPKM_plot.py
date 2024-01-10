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


# RPKM analysis
gID_field_values, field_gID_values = load_file.read_GTF ("ENCFF159KBI.gtf", chr_list=['chr1'], mode="both")
"""
gID_exons = field_gID_values['exons']
gID_exonlen = {}
for gID in gID_exons:
    length = 0
    for start, end in gID_exons[gID]:
        length +=  end - start + 1
    gID_exonlen[gID] = length

gID_exp_counts, exp_gID_counts = load_file.read_tabular_file ("/home/spark159/../../media/spark159/sw/dataforcondense/GSE63124_all_gene_raw_readcounts.txt", mode="both")
gID_counts1 = exp_gID_counts['38-Per_rep1']
gID_counts2 = exp_gID_counts['38-Per_rep2']

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
"""
gID_RPKM = read_tsv("ENCFF174OMR.tsv")
#gID_RPKM = read_tsv("ENCFF910OBU.tsv")
#gID_RPKM = read_tsv("ENCFF395XDK.tsv")

#temp_dict = {}
#for gID in gID_RPKM:
#    if gID in gID_field_values:
#        RPKM = gID_RPKM[gID]
#        if RPKM <= 0:
#            continue
#        temp_dict[gID] = gID_RPKM[gID]
#gID_RPKM = temp_dict

gID_ginterval = {}
for gID in gID_RPKM.keys():
    try:
        TSS = gID_field_values[gID]['TSS']
    except:
        continue
    #TTS = gID_field_values[gID]['TTS']
    strand = gID_field_values[gID]['strand']
    interval = (TSS-2500, TSS+2500)
    #interval = (TSS-500, TSS+500)
    #if strand == '+':
        #interval = (TSS, TTS)
        #interval = (TSS-250, TSS)
        #interval = (TSS, TSS+2500)
    #else:
        #interval = (TTS, TSS)
        #interval = (TSS, TSS+250)
        #interval = (TSS-2500, TSS)
    gID_ginterval[gID] = interval

ginterval_dict = Interval_dict.double_hash(gID_ginterval, 100000, 250000000)

nID_chr, nID_pos, name_nID_value = load_file.read_anot_file("H1_NCP_sp_chr1_anot.cn")
#nID_chr, nID_pos, name_nID_value = load_file.read_anot_file("/home/spark159/../../media/spark159/sw/sp_spd_tests_detail/hg19_chr1_167win25step_anot.cn")

nID_CG = name_nID_value['CNumber(CpG)']
nID_me = name_nID_value['meCNumber(CpG)']
nID_mefrac = {}
for nID in nID_CG:
    CG = nID_CG[nID]
    if CG <= 0:
        continue
    me = nID_me[nID]
    mefrac = float(me) / CG
    nID_mefrac[nID] = mefrac
name_nID_value['meCpG density'] = nID_mefrac

name_gID_values = {}
for nID in nID_pos:
    pos = nID_pos[nID]
    gIDs = ginterval_dict.find(pos)
    if len(gIDs) <= 0:
        continue
    for name in name_nID_value:
        if name == 'Sequence':
            continue
        if name not in name_gID_values:
            name_gID_values[name] = {}
        try:
            value = name_nID_value[name][nID]
        except:
            continue
        if np.isnan(value):
            continue
        for gID in gIDs:
            if gID not in name_gID_values[name]:
                name_gID_values[name][gID] = []
            name_gID_values[name][gID].append(value)

name_gID_mean = {}
for name in name_gID_values:
    if name not in name_gID_mean:
        name_gID_mean[name] = {}
    for gID in name_gID_values[name]:
        assert gID not in name_gID_mean[name]
        name_gID_mean[name][gID] = np.mean(name_gID_values[name][gID])


names = name_gID_mean.keys()
X_list, Y_list = [], []
for i in range(len(names)):
    name = names[i]
    X, Y = [], []
    for gID in gID_RPKM:
        #counts = gID_counts[gID]
        #if counts <= 0:
            #continue
        RPKM = gID_RPKM[gID]
        if RPKM <= 0:
            continue
        #logRPKM = RPKM
        logRPKM = math.log(RPKM)
        #logRPKM = math.log(RPKM+1)
        try:
            mean = name_gID_mean[name][gID]
        except:
            continue
        #X.append(RPKM)
        X.append(logRPKM)
        #if name.startswith('k'):
        #    Y.append(math.log(mean+1))
        #else:
        #    Y.append(mean+1)
        Y.append(mean)
    X_list.append(X)
    Y_list.append(Y)


corr_list = []
name_list = []
for i in range(len(names)):
    name = names[i].split('/')[-1]
    #if name == 'H1-NCP-sp4':
    #    continue
    #if name == 'H1-NCP-sp8':
    #    name = 'Condensability'
    name_list.append(name)
    X, Y = X_list[i], Y_list[i]
    corr = get_corr(X, Y)
    corr_list.append(corr)
    print name, corr
    feature_list = [[x] for x in X]
    test_list = [[y] for y in Y]
    reg = linear_model.Ridge(alpha=0.5)
    reg.fit (feature_list, test_list)
    Ypred = reg.predict(feature_list)
    Ypred = [ value[0] for value in Ypred]
    fig = plt.figure()
    plt.plot(X, Y, '.', alpha=0.5)
    plt.plot(X, Ypred, 'r--')
    xloc, yloc = np.mean([np.median(X), max(X)]), np.mean([np.median(Y), max(Y)])
    plt.text(xloc, yloc, str(round(corr,3)), fontsize=20, va='center', ha='center')
    #plt.xlabel("logRPKM")
    plt.xlabel("logFPKM")
    #plt.xscale('log', basex=2)
    plt.ylabel(name)
    plt.title("Near TSS (5kb)")
    plt.savefig("RPKM_nearTSS_scatter_" + name + ".png",bbox_inches='tight')
    #plt.show()
    plt.close()


fig = plt.figure()
plt.bar(range(len(corr_list)), corr_list, width=0.5, color='g')
plt.xticks(range(len(corr_list)), name_list, rotation=90)
plt.ylabel("Pearson Correlation")
plt.axhline(y=0, color='k', linestyle='--')
plt.title("Correlation with FPKM")
plt.savefig("bar_corr.png",bbox_inches='tight')
#plt.show()
plt.close()
