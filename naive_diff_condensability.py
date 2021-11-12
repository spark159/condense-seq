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


# chr list for analysis
chr_list = ['chr' + str(i) for i in range(1, 23)]
chr_list += ['chrX', 'chrY']
#chr_list = ['chr9']
#chr_list.remove('chr13')

# stem cell marker genes
ESC_tf_cores =  ['Pou5f1', 'Sox2', 'KLF4', 'Nanog']
ESC_tf_others = ['Zfp42', 'UTF1', 'ZFX', 'TBN', 'FoxD3', 'HMGA2', 'NAC1', 'NR6A1', 'Stat3', 'LEF1', 'TCF3', 'Sall4', 'Fbxo15', 'L1TD1', 'Gdf3', 'Dppa5', 'Dppa4', 'Dppa2', 'Dppa3']
ESC_tf_cores = [gname.upper() for gname in ESC_tf_cores]
ESC_tf_others = [gname.upper() for gname in ESC_tf_others]
ESC_gnames = ESC_tf_cores + ESC_tf_others

# read gene expression file
gID_FPKM1 = read_tsv("ENCFF174OMR.tsv")
gID_FPKM2 = load_file.read_RPKM("GSE63124_all_gene_raw_readcounts.txt", "Homo_sapiens.GRCh37.87.gtf", chr_list = chr_list)

# read gtf file
gID_field_values = load_file.read_GTF ("ENCFF159KBI.gtf", chr_list=chr_list, mode="gene")

# calculate score for all genes
gID_mscore1 = {}
gID_mscore2 = {}
for chr in chr_list:
    print "processing " + chr
    
    fname1 = '%s_gene_score1' % (chr)
    fname2 = '%s_gene_score2' % (chr)

    try:
        chr_gID_mscore1 = pickle.load(open(fname1 + ".pickle", "rb"))
        chr_gID_mscore2 = pickle.load(open(fname2 + ".pickle", "rb"))

    except:

        # read gene locations
        chr_gID_field_values1 = load_file.read_GTF ("ENCFF159KBI.gtf", chr_list=[chr], mode="gene")
        chr_gID_field_values2 = load_file.read_GTF_old ("Homo_sapiens.GRCh37.87.gtf", chr_list=[chr], mode="gene")

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
        anot_fname1 = "H1_NCP_sp_%s_anot.cn" % (chr)
        nID_chr1, nID_pos1, name_nID_value1 = load_file.read_anot_file(anot_fname1)
        nID_score1 = name_nID_value1['work/2021_06_07_H1_sp_detail/H1-NCP-sp-8']
        del nID_chr1, name_nID_value1

        anot_fname2 = "hg19_%s_NCP_ics_anot.cn" % (chr)
        nID_chr2, nID_pos2, name_nID_value2 = load_file.read_anot_file(anot_fname2)
        nID_score2 = name_nID_value2['data/sp_spd_tests_detail/sp8']
        del nID_chr2, name_nID_value2

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
        pickle.dump(chr_gID_mscore1, open(fname1 + ".pickle", "wb"))
        pickle.dump(chr_gID_mscore2, open(fname2 + ".pickle", "wb"))

    # update total data set
    gID_mscore1.update(chr_gID_mscore1)
    gID_mscore2.update(chr_gID_mscore2)

    print len(set(chr_gID_mscore1.keys()) & set(chr_gID_mscore2.keys()))

    print 

# finalize common gene set
gIDs = list(set(gID_mscore1) & set(gID_mscore2) & set(gID_FPKM1.keys()) & set(gID_FPKM2.keys()))
print 'Total gene count:' + str(len(gIDs))

# find Ensemble Gene ID for ESC marker genes
ESC_gname_gIDs = {gname :[] for gname in ESC_gnames}
for gID in gIDs:
    gname = gID_field_values[gID]['geneName']
    try:
        ESC_gname_gIDs[gname].append(gID)
    except:
        pass
    
ESC_gID_gname = {}
for gname in ESC_gname_gIDs:
    if len(ESC_gname_gIDs[gname]) == 1:
        gID = ESC_gname_gIDs[gname][0]
        assert gID not in ESC_gID_gname
        ESC_gID_gname[gID] = gname

# set in-IDs and out-IDs
in_gIDs = list(set(gIDs) - set(ESC_gID_gname.keys())) # all genes except ESC marker
out_gIDs = ESC_gID_gname.keys() # ESC marker genes


# H1 FPKM vs Panc FPKM
X, Y = [], []
for gID in gIDs:
    X.append(np.log2(1+gID_FPKM1[gID]))
    Y.append(np.log2(1+gID_FPKM2[gID]))
fig = plt.figure()
plt.plot(X, Y, '.')
plt.title('RNA-seq data comparison')
plt.xlabel('H1 hESC logFPKM')
plt.ylabel('A38-41 hPanc logFPKM')
#plt.show()
plt.close()

# H1 FPKM vs H1 score
inX = [np.log2(1+gID_FPKM1[gID]) for gID in in_gIDs]
inY = [gID_mscore1[gID] for gID in in_gIDs]
outX = [np.log2(1+gID_FPKM1[gID]) for gID in out_gIDs]
outY = [gID_mscore1[gID] for gID in out_gIDs]
gnames = [ESC_gID_gname[gID] for gID in out_gIDs]

fig = plt.figure()
plt.plot(inX, inY, 'k,', markersize=3, alpha=0.3)
for x, y, gname in zip(outX, outY, gnames):
    plt.plot(x, y, 'r.', markersize=3, alpha=1, zorder=10)
    if gname in ESC_tf_cores:
        plt.annotate(gname, (x, y), color='red', zorder=40, size=5)
plt.title('H1 hESC')
plt.xlabel('Gene expression (logFPKM)')
plt.ylabel('Condensabiltiy (A.U.)')
plt.show()
plt.close()

# Panc FPKM vs Panc score
inX = [np.log2(1+gID_FPKM2[gID]) for gID in in_gIDs]
inY = [gID_mscore2[gID] for gID in in_gIDs]
outX = [np.log2(1+gID_FPKM2[gID]) for gID in out_gIDs]
outY = [gID_mscore2[gID] for gID in out_gIDs]
gnames = [ESC_gID_gname[gID] for gID in out_gIDs]

fig = plt.figure()
plt.plot(inX, inY, 'k,', markersize=3, alpha=0.3)
for x, y, gname in zip(outX, outY, gnames):
    plt.plot(x, y, 'r.', markersize=3, alpha=1, zorder=10)
    if gname in ESC_tf_cores:
        plt.annotate(gname, (x, y), color='red', zorder=40, size=5)
plt.title('A38-41 hPanc')
plt.xlabel('Gene expression (logFPKM)')
plt.ylabel('Condensabiltiy (A.U.)')
plt.show()
plt.close()


## H1 score VS Panc score
inX = [gID_mscore1[gID] for gID in in_gIDs]
inY = [gID_mscore2[gID] for gID in in_gIDs]
outX = [gID_mscore1[gID] for gID in out_gIDs]
outY = [gID_mscore2[gID] for gID in out_gIDs]
gnames = [ESC_gID_gname[gID] for gID in out_gIDs]

fig = plt.figure()
plt.plot(inX, inY, 'k.', markersize=3, alpha=0.3)
for x, y, gname in zip(outX, outY, outgnames):
    plt.plot(x, y, 'r.', markersize=3, alpha=1, zorder=10)
    if gname in ESC_tf_cores:
        plt.annotate(gname, (x, y), color='red', zorder=40, size=5)
plt.title("Condensability near TSS (5kb)")
plt.xlabel('H1 hESC')
plt.ylabel('A38-41 hPanc')
plt.show()
plt.close()

sys.exit(1)

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
plt.ylabel('A38-41 hPanc - H1 hESC')
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
plt.xlabel("Condensability difference (A38-41 hPanc - H1 hESC)")
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
