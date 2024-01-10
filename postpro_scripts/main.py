import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math

"""
# gene type analysis
genome_size = load_file.read_genome_size("data/hg19.fa")
ID_feature_intervals, feature_ID_intervals = load_file.read_GTF ("data/Homo_sapiens.GRCh37.87.gtf", "chr1", mode="both")
TSS_range = 500
TTS_range = 500
Prom_range = 5000
dID_interval = {}
for ID in ID_feature_intervals:
    feature_intervals = ID_feature_intervals[ID]
    if feature_intervals['geneType'] != "protein_coding":
        continue
    TSS = feature_intervals['TSS']
    TTS = feature_intervals['TTS']
    TSS_st = max(TSS - TSS_range, 0)
    TSS_ed = min(TSS + TSS_range, genome_size['chr1'])
    TTS_st = max(TTS - TTS_range, 0) 
    TTS_ed = min(TTS + TTS_range, genome_size['chr1'])
    strand = feature_intervals['strand']
    if strand == '+':
       Body_st = TSS_ed
       Body_ed = TTS_st
       Prom_st = max(TSS - Prom_range, 0)
       Prom_ed = TSS_st
    if strand == '-':
       Body_st = TTS_ed
       Body_ed = TSS_st
       Prom_st = TSS_ed
       Prom_ed = min(TSS + Prom_range, genome_size['chr1'])
    if TSS_ed - TSS_st > 0:
        dID_interval["TSS:" + ID] = [TSS_st, TSS_ed]
    if TTS_ed - TTS_st > 0:
        dID_interval["TTS:" + ID] = [TTS_st, TTS_ed]
    if Body_ed - Body_st > 0:
        dID_interval["Body:" + ID] = [Body_st, Body_ed]
    if Prom_ed - Prom_st > 0:
        dID_interval["Prom:" + ID] = [Prom_st, Prom_ed]
dinterval_dict = Interval_dict.double_hash(dID_interval, 100000, genome_size['chr1'])

nID_chr, nID_pos, name_nID_value = load_file.read_anot_file("data/hg19_chr1_167win25step_anot.cn")
nID_score1 = name_nID_value['work/condense_seq/sp9_hg19_chr1']
nID_score2 = name_nID_value['work/condense_seq/sp10_hg19_chr1']
nID_AT = name_nID_value['ATcontent']
nID_meGC = name_nID_value['meGCNumber']
nID_CpG = name_nID_value['CpGNumber']
nID_k27ac = name_nID_value['k27ac']
nID_k27me3 = name_nID_value['k27me3']
nID_k36me3 = name_nID_value['k36me3']
nID_k4me3 = name_nID_value['k4me3']
nID_k9ac = name_nID_value['k9ac']
nID_k9me2 = name_nID_value['k9me2']
nID_k9me3 = name_nID_value['k9me3']


#frac = [(4**i) for i in range(1,11)]
#frac = frac[::-1]
#group1 = statis.quantile(nID_score1, len(frac), frac=frac)
#group2 = statis.quantile(nID_score2, len(frac), frac=frac)

mark_domain_values = {}

domain_list1 = {}
domain_list2 = {}
domain_AT = {}
domain_meGC = {}
domain_CpG = {}
domain_k27ac = {}
domain_k27me3 = {}
domain_k36me3 = {}
domain_k4me3 = {}
domain_k9ac = {}
domain_k9me2 = {}
domain_k9me3 = {}

for nID in nID_score1:
    score1 = nID_score1[nID]
    score2 = nID_score2[nID]
    pos = nID_pos[nID]
    AT = nID_AT[nID]
    IDs = dinterval_dict.find(pos)
    if len(IDs) <= 0:
        if "Intergenic" not in domain_list1:
            domain_list1["Intergenic"] = []
            domain_list2["Intergenic"] = []
            domain_AT["Intergenic"] = []
            domain_meGC["Intergenic"] = []
            domain_CpG["Intergenic"] = []
            domain_k27ac["Intergenic"] = []
            domain_k27me3["Intergenic"] = []
            domain_k36me3["Intergenic"] = []
            domain_k4me3["Intergenic"] = []
            domain_k9ac["Intergenic"] = []
            domain_k9me2["Intergenic"] = []
            domain_k9me3["Intergenic"] = []
        domain_list1["Intergenic"].append(score1)
        domain_list2["Intergenic"].append(score2)
        domain_AT["Intergenic"].append(AT)
        domain_meGC["Intergenic"].append(nID_meGC[nID])
        domain_CpG["Intergenic"].append(nID_CpG[nID])
        domain_k27ac["Intergenic"].append(nID_k27ac[nID])
        domain_k27me3["Intergenic"].append(nID_k27me3[nID])
        domain_k36me3["Intergenic"].append(nID_k36me3[nID])
        domain_k4me3["Intergenic"].append(nID_k4me3[nID])
        domain_k9ac["Intergenic"].append(nID_k9ac[nID])
        domain_k9me2["Intergenic"].append(nID_k9me2[nID])
        domain_k9me3["Intergenic"].append(nID_k9me3[nID])
        continue
    for ID in IDs:
        domain = ID.split(':')[0]
        if domain not in domain_list1:
            domain_list1[domain] = []
            domain_list2[domain] = []
            domain_AT[domain] = []
            domain_meGC[domain] = []
            domain_CpG[domain] = []
            domain_k27ac[domain] = []
            domain_k27me3[domain] = []
            domain_k36me3[domain] = []
            domain_k4me3[domain] = []
            domain_k9ac[domain] = []
            domain_k9me2[domain] = []
            domain_k9me3[domain] = []
        domain_list1[domain].append(score1)
        domain_list2[domain].append(score2)
        domain_AT[domain].append(AT)
        domain_meGC["Intergenic"].append(nID_meGC[nID])
        domain_CpG["Intergenic"].append(nID_CpG[nID])
        domain_k27ac["Intergenic"].append(nID_k27ac[nID])
        domain_k27me3["Intergenic"].append(nID_k27me3[nID])
        domain_k36me3["Intergenic"].append(nID_k36me3[nID])
        domain_k4me3["Intergenic"].append(nID_k4me3[nID])
        domain_k9ac["Intergenic"].append(nID_k9ac[nID])
        domain_k9me2["Intergenic"].append(nID_k9me2[nID])
        domain_k9me3["Intergenic"].append(nID_k9me3[nID])

fig = plt.figure()


#fig = plt.figure()
#domains = domain_list1.keys()
#pos_list1 = [i - 0.2 for i in range(len(domains))]
#pos_list2 = [i + 0.2 for i in range(len(domains))]
#plt.boxplot([domain_AT[domain] for domain in domains], 0, "", positions=pos_list1, widths=0.3)
#plt.boxplot([domain_list1[domain] for domain in domains], 0, "", positions=pos_list1, widths=0.3)
#plt.boxplot([domain_list2[domain] for domain in domains], 0, "", positions=pos_list2, widths=0.3)
#plt.xticks(range(len(domains)), domains)
#plt.xlim([-0.5, len(domains)-0.5])
#plt.show()
#plt.close()


domain_meanscore1 = {}
domain_meanscore2 = {}
domain_error1 = {}
domain_error2 = {}

for domain in domain_list1:
    domain_meanscore1[domain] = np.mean(domain_list1[domain])
    domain_meanscore2[domain] = np.mean(domain_list2[domain])
    domain_error1[domain] = np.std(domain_list1[domain])/np.sqrt(len(domain_list1[domain]))
    domain_error2[domain] = np.std(domain_list2[domain])/np.sqrt(len(domain_list2[domain]))

X = []
Y1, Y2 = [], []
for ID in ID_RPKM1:
    counts1 = ID_counts1[ID]
    if counts1 <= 0:
        continue
    RPKM1 = ID_RPKM1[ID]
    logRPKM1 = math.log(RPKM1)
    try:
        score1 = ID_meanscore1[ID]
        score2 = ID_meanscore2[ID]
    except:
        continue
    #X.append(RPKM1)
    X.append(logRPKM1)
    #X.append(np.log2(counts1))
    Y1.append(score1)
    Y2.append(score2)

fig = plt.figure()
plt.plot(X, Y1, '.', label='sp9')
plt.plot(X, Y2, '.', label='sp10')
plt.xlabel("logRPKM")
plt.ylabel("condensability")
plt.legend()
plt.show()
plt.close()



# RPKM analysis
ID_feature_intervals, feature_ID_intervals = load_file.read_GTF ("data/Homo_sapiens.GRCh37.87.gtf", "chr1", mode="both")

ID_exons = feature_ID_intervals['exons']
ID_exonlen = {}
for ID in ID_exons:
    length = 0
    for start, end in ID_exons[ID]:
        length +=  end - start + 1
    ID_exonlen[ID] = length

ID_exp_counts, exp_ID_counts = load_file.read_tabular_file ("data/GSE63124_all_gene_raw_readcounts.txt", mode="both")
ID_counts1 = exp_ID_counts['38-Per_rep1']
ID_counts2 = exp_ID_counts['38-Per_rep2']

ID_RPKM1 = {}
ID_RPKM2 = {}
for ID in ID_exonlen:
    try:
        RPKM1 = float(ID_counts1[ID])/ID_exonlen[ID]
        RPKM2 = float(ID_counts2[ID])/ID_exonlen[ID]
    except:
        continue
    ID_RPKM1[ID] = RPKM1
    ID_RPKM2[ID] = RPKM2

ID_ginterval = {}
for ID in ID_RPKM1.keys():
    TSS = ID_feature_intervals[ID]['TSS']
    TTS = ID_feature_intervals[ID]['TTS']
    strand = ID_feature_intervals[ID]['strand']
    interval = (TSS-500, TSS+500)
    #if strand == '+':
        #interval = (TSS, TTS)
        #interval = (TSS-250, TSS)
    #else:
        #interval = (TTS, TSS)
        #interval = (TSS, TSS+250)
    ID_ginterval[ID] = interval

ginterval_dict = Interval_dict.double_hash(ID_ginterval, 10000, 250000000)

nID_chr, nID_pos, name_nID_value = load_file.read_anot_file("data/hg19_chr1_167win25step_anot.cn")
nID_score1 = name_nID_value['work/condense_seq/sp9_hg19_chr1']
nID_score2 = name_nID_value['work/condense_seq/sp10_hg19_chr1']
#field_nID_score = load_file.read_tabular_file ("data/score.cn", mode="col")
#nID_score1 = field_nID_value['work/condense_seq/sp9_hg19_chr1']
#nID_score2 = field_nID_value['work/condense_seq/sp10_hg19_chr1']

ID_list1 = {}
ID_list2 = {}

for nID in nID_score1:
    score1 = nID_score1[nID]
    score2 = nID_score2[nID]
    pos = nID_pos[nID]
    IDs = ginterval_dict.find(pos)
    for ID in IDs:
        if ID not in ID_list1:
            ID_list1[ID] = []
        ID_list1[ID].append(score1)
        if ID not in ID_list2:
            ID_list2[ID] = []
        ID_list2[ID].append(score2)

ID_meanscore1 = {}
ID_meanscore2 = {}

for ID in ID_list1:
    ID_meanscore1[ID] = np.mean(ID_list1[ID])
    ID_meanscore2[ID] = np.mean(ID_list2[ID])

X = []
Y1, Y2 = [], []
for ID in ID_RPKM1:
    counts1 = ID_counts1[ID]
    if counts1 <= 0:
        continue
    RPKM1 = ID_RPKM1[ID]
    logRPKM1 = math.log(RPKM1)
    try:
        score1 = ID_meanscore1[ID]
        score2 = ID_meanscore2[ID]
    except:
        continue
    #X.append(RPKM1)
    X.append(logRPKM1)
    #X.append(np.log2(counts1))
    Y1.append(score1)
    Y2.append(score2)

fig = plt.figure()
plt.plot(X, Y1, '.', label='sp9')
plt.plot(X, Y2, '.', label='sp10')
plt.xlabel("logRPKM")
plt.ylabel("Condensability")
plt.legend()
plt.show()
plt.close()
"""

"""
#GC/me/PTM analysis
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("data/hg19_chr1_171_everything_anot.cn")
#ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("hg19_chr1_171_anot.cn")
ID_score1 = name_ID_value['work/condense_seq/sp9_hg19_chr1']
ID_score2 = name_ID_value['work/condense_seq/sp10_hg19_chr1']
ID_AT = name_ID_value['ATcontent']
ID_CpG = name_ID_value['CpGNumber']
ID_me = name_ID_value['meGCNumber']

ID_chip_list, chip_names = [], []
for name in name_ID_value:
    if name.startswith('k'):
        ID_value = name_ID_value[name]
        chip_names.append(name)
        ID_chip_list.append(ID_value)

for ID in ID_AT:
    ID_AT[ID] = ID_AT[ID]*100

ID_meCpGfrac = {}
ID_meGCfrac = {}
for ID in ID_CpG:
    CpG = ID_CpG[ID]
    me = ID_me[ID]
    GC = 100-ID_AT[ID]
    if CpG <= 0:
        value = 0.0
    else:
        float(me)/(171*CpG)
    ID_meCpGfrac[ID] = value
    ID_meGCfrac[ID] = float(me)/GC


new_ID_score2 = statis.neutralize_score_by_target(ID_score2, ID_AT)

#graphics.Scatter_plot(ID_AT, ID_score1, note='sp9')
graphics.Scatter_plot(ID_AT, ID_score2, note='sp10')

frac = [(4**i) for i in range(1,11)]
frac = frac[::-1]

#graphics.PartitionMeanplot(ID_AT, ID_score1, ID_me, frac, note="sp9_me")
graphics.PartitionMeanplot(ID_AT, ID_score2, ID_me, frac, note="sp10_me")
graphics.PartitionScatterplot(ID_AT, ID_score2, ID_me, frac, note="sp10_me")
graphics.PartitionBoxplot(new_ID_score2, ID_me, frac, note="sp10_me")


for i in range(len(ID_chip_list)):
    ID_chip = ID_chip_list[i]
    name = chip_names[i]
    #graphics.PartitionMeanplot(ID_AT, ID_score1, ID_chip, frac, note="sp9_" + name)
    graphics.PartitionMeanplot(ID_AT, ID_score2, ID_chip, frac, note="sp10_" + name)
    graphics.PartitionScatterplot(ID_AT, ID_score2, ID_chip, frac, note="sp10_" + name)
    graphics.PartitionBoxplot(new_ID_score2, ID_chip, frac, note="sp10_" + name)
"""
