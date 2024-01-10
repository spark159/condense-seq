import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math

def pair_boxplot (key_values1, key_values2, ylabel1='', ylabel2='Condensability (A.U.)', title=None, keys=None, note=""):
    if keys:
        assert set(keys) <= set(key_values1.keys())
        assert set(keys) <= set(key_values2.keys())
    else:
        assert set(key_values1,keys()) == set(key_values2.keys())
        keys = key_values1.keys()
    fig, ax1 = plt.subplots()
    pos_list1 = [i - 0.2 for i in range(len(keys))]
    bp1 = ax1.boxplot([key_values1[key] for key in keys], 0, "", positions=pos_list1, widths=0.3, patch_artist=True, boxprops=dict(facecolor="pink"))
    ax1.set_ylabel(ylabel1, color='r')
    ax1.tick_params('y', colors='r')
    ax2 = ax1.twinx()
    pos_list2 = [i + 0.2 for i in range(len(keys))]
    bp2 = ax2.boxplot([key_values2[key] for key in keys], 0, "", positions=pos_list2, widths=0.3, patch_artist=True, boxprops=dict(facecolor="lightblue"))
    ax2.set_ylabel(ylabel2, color='b')
    ax2.tick_params('y', colors='b')
    plt.legend([bp1["boxes"][0], bp2["boxes"][0]], [ylabel1, ylabel2], loc='upper right')
    if title:
        plt.title(title)
    plt.xticks(range(len(keys)), keys)
    plt.xlim([-0.5, len(keys)-0.5])
    fig.tight_layout()
    plt.savefig("pairbox_" + note + ".png",bbox_inches='tight')
    #plt.show()
    plt.close()

# genomic domain analysis
genome_size = load_file.read_genome_size("data/hg19.fa")
ID_field_values, field_ID_values = load_file.read_GTF ("data/Homo_sapiens.GRCh37.87.gtf", "chr1", mode="both")

geneType_num = {}
ID_geneType = field_ID_values['geneType']
for ID in ID_geneType:
    geneType = ID_geneType[ID]
    if geneType not in geneType_num:
        geneType_num[geneType] = 0
    geneType_num[geneType] += 1

num_geneType = []
for geneType, num in geneType_num.items():
    num_geneType.append([num, geneType])
num_geneType = sorted(num_geneType, reverse=True)

for num, geneType in num_geneType:
    print geneType, num


"""
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

#nID_chr, nID_pos, name_nID_value = load_file.read_anot_file("data/hg19_chr1_167win25step_anot.cn")
nID_chr, nID_pos, name_nID_value = load_file.read_anot_file("data/hg19_chr1_171_everything_anot.cn")

name_domain_values = {}
for name in name_nID_value:
    domain_values = {}
    nID_value = name_nID_value[name]
    for nID in nID_value:
        pos = nID_pos[nID]
        value = nID_value[nID]
        IDs = dinterval_dict.find(pos)
        if len(IDs) <= 0:
            if "Intergenic" not in domain_values:
                domain_values["Intergenic"] = []
            domain_values["Intergenic"].append(value)
        for ID in IDs:
            domain = ID.split(':')[0]
            if domain not in domain_values:
                domain_values[domain] = []
            domain_values[domain].append(value)
    assert name not in name_domain_values
    name_domain_values[name] = domain_values

domain_score2 = name_domain_values["work/condense_seq/sp10_hg19_chr1"]
for name in name_domain_values:
    if name == "work/condense_seq/sp10_hg19_chr1":
        continue
    domain_values = name_domain_values[name]
    pair_boxplot(domain_values, domain_score2, ylabel1=name, keys=['Prom', 'TSS', 'Body', 'TTS', 'Intergenic'], title='Protein coding genes', note=name.split('/')[-1])
"""
