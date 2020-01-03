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

def read_bincountfile (fname, chr_list=None):
    First = True
    for line in open(fname):
        if First:
            cols = line.strip().split()
            names = [name.rsplit('.')[-2] for name in cols[4:-1]]
            chr_binID_counts = [{} for i in range(len(names))]
            chr_binID_range = {}
            chr_binID_GC = {}
            First = False
            continue
        line = line.strip()
        if not line:
            continue
        cols = line.strip().split()
        ID, chr, st, ed = cols[:4]
        if chr_list != None and chr not in chr_list:
            continue
        ID = int(ID)
        st, end = int(st), int(ed)
        GC = float(cols[-1])
        if chr not in chr_binID_range:
            chr_binID_range[chr] = []
        chr_binID_range[chr].append((st, ed))
        if chr not in chr_binID_GC:
            chr_binID_GC[chr] = []
        chr_binID_GC[chr].append(GC)
        datas = cols[4:-1]
        for i in range(len(datas)):
            data = float(datas[i])
            chr_binID_count = chr_binID_counts[i]
            if chr not in chr_binID_count:
                chr_binID_count[chr] = []
            chr_binID_count[chr].append(data)
    return names, chr_binID_counts, chr_binID_range, chr_binID_GC

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
gID_field_values, field_gID_values = load_file.read_GTF ("/home/spark159/../../media/spark159/sw/dataforcondense/Homo_sapiens.GRCh37.87.gtf", mode="both")

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
    #if RPKM <= 0:
    #    continue
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

fname = "NCP_Spermidine(3+)_1kb"
bin_size = 1000
names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile("/home/spark159/Downloads/" + fname + "_bin.cn")
chr_binID_control = chr_binID_counts[-1]

gID_rcounts = []
gID_GC = {}
for i in range(len(names)-1):
    chr_binID_count = chr_binID_counts[i]
    gID_rcount = {}
    for gID in gID_RPKM:
        chr = gID_field_values[gID]['chr']
        if chr not in chr_binID_count:
            continue
        st, ed = gID_ginterval[gID]
        st_binID, ed_binID = st / bin_size, ed / bin_size
        if st_binID == ed_binID:
            control = chr_binID_control[chr][st_binID]
            if control <= 0:
                continue
            value = chr_binID_count[chr][st_binID] / control
            GC = chr_binID_GC[chr][st_binID]
        else:
            total = 0.0
            control = 0.0
            GC = 0.0
            for k in range(st_binID, ed_binID+1):
                if k == st_binID:
                    control += ((k+1)*bin_size - st) * chr_binID_control[chr][k]
                    total += ((k+1)*bin_size - st) * chr_binID_count[chr][k]
                    GC += ((k+1)*bin_size - st) * chr_binID_GC[chr][k]
                elif k == ed_binID:
                    control += (ed - k*bin_size) * chr_binID_control[chr][k]
                    total += (ed - k*bin_size) * chr_binID_count[chr][k]
                    GC += (ed - k*bin_size) * chr_binID_GC[chr][k] 
                else:
                    control += bin_size * chr_binID_control[chr][k]
                    total += bin_size * chr_binID_count[chr][k]
                    GC += bin_size * chr_binID_GC[chr][k]
            if control <= 0:
                continue
            value = float(total) / control
            GC = float(GC) / (ed - st)
        gID_rcount[gID] = value
        gID_GC[gID] = GC
    gID_rcounts.append(gID_rcount)

for i in range(len(names)-1):
    gID_rcount = gID_rcounts[i]
    X, Y = [], []
    Z = []
    for gID in gID_rcount:
        if gID_rcount[gID] <= 0:
            continue
        X.append(np.log(gID_rcount[gID]))
        Y.append(np.log(gID_RPKM[gID]))
        Z.append(gID_GC[gID])

    feature_list = [[x] for x in X]
    test_list = [[y] for y in Y]
    reg = linear_model.Ridge(alpha=0.5)
    reg.fit (feature_list, test_list)
    Ypred = reg.predict(feature_list)
    Ypred = [ value[0] for value in Ypred]
    corr = get_corr(X, Y)    
    fig = plt.figure()
    plt.text(0,0,str(round(corr,3)), fontsize=25)
    plt.plot(X, Y, 'g,', alpha=0.5)
    plt.plot(X, Ypred, 'r--')
    plt.xlabel("log Normalized Counts")
    plt.ylabel("log RPKM")
    plt.title("Titration " + str(i+1))
    plt.savefig(fname+"_"+str(i+1)+".png")
    #plt.show()
    plt.close()

                
"""
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
"""
