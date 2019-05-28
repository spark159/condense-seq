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

# RPKM analysis
gID_field_values, field_gID_values = load_file.read_GTF ("data/Homo_sapiens.GRCh37.87.gtf", "chr1", mode="both")

gID_exons = field_gID_values['exons']
gID_exonlen = {}
for gID in gID_exons:
    length = 0
    for start, end in gID_exons[gID]:
        length +=  end - start + 1
    gID_exonlen[gID] = length

gID_exp_counts, exp_gID_counts = load_file.read_tabular_file ("data/GSE63124_all_gene_raw_readcounts.txt", mode="both")
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

gID_ginterval = {}
for gID in gID_RPKM.keys():
    TSS = gID_field_values[gID]['TSS']
    TTS = gID_field_values[gID]['TTS']
    strand = gID_field_values[gID]['strand']
    interval = (TSS-500, TSS+500)
    #if strand == '+':
        #interval = (TSS, TTS)
        #interval = (TSS-250, TSS)
        #interval = (TSS, TSS+2500)
    #else:
        #interval = (TTS, TSS)
        #interval = (TSS, TSS+250)
        #interval = (TSS-2500, TSS)
    gID_ginterval[gID] = interval

ginterval_dict = Interval_dict.double_hash(gID_ginterval, 10000, 250000000)

#nID_chr, nID_pos, name_nID_value = load_file.read_anot_file("data/hg19_chr1_167win25step_anot.cn")
nID_chr, nID_pos, name_nID_value = load_file.read_anot_file("data/hg19_chr1_171_everything_anot.cn")
#nID_score1 = name_nID_value['work/condense_seq/sp9_hg19_chr1']
#nID_score2 = name_nID_value['work/condense_seq/sp10_hg19_chr1']

name_gID_values = {}
for nID in nID_pos:
    pos = nID_pos[nID]
    gIDs = ginterval_dict.find(pos)
    if len(gIDs) <= 0:
        continue
    for name in name_nID_value:
        if name not in name_gID_values:
            name_gID_values[name] = {}
        value = name_nID_value[name][nID]
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
        counts = gID_counts[gID]
        #if counts <= 0:
            #continue
        RPKM = gID_RPKM[gID]
        logRPKM = math.log(RPKM)
        try:
            mean = name_gID_mean[name][gID]
        except:
            continue
        X.append(logRPKM)
        Y.append(mean)
    X_list.append(X)
    Y_list.append(Y)


corr_list = []
for i in range(len(names)):
    name = names[i].split('/')[-1]
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
    plt.xlabel("logRPKM")
    plt.ylabel(name)
    plt.title("Near TSS (1kb)")
    plt.savefig("RPKM_nearTSS_scatter_" + name + ".png",bbox_inches='tight')
    #plt.show()
    plt.close()


#fig = plt.figure()
#plt.bar(range(len(corr_list)), corr_list)
#plt.xticks(range(len(corr_list)), names, rotation=45)
#plt.savefig("bar_corr.png",bbox_inches='tight')
#plt.show()
#plt.close()
