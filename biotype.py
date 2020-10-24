import sys
import Interval_dict
import numpy as np
import matplotlib.pyplot as plt
import load_file
import statis

path = './data/'

# read GTF file
gID_field_values, field_gID_values = load_file.read_GTF(path+"gencode.v19.annotation.gtf_withproteinids", chr_list=["chr1"], mode="both")

# sort by biotype
gID_biotype = field_gID_values['geneType']

biotype_gIDs = {}
for gID, biotype in gID_biotype.items():
    if biotype not in biotype_gIDs:
        biotype_gIDs[biotype] = []
    biotype_gIDs[biotype].append(gID)

selected_biotypes = []
for biotype in biotype_gIDs:
    if len(biotype_gIDs[biotype]) < 50:
        continue
    print biotype, len(biotype_gIDs[biotype])
    selected_biotypes.append(biotype)

gID_ginterval = {}

gID_ginterval = {}
for gID in gID_field_values:
    TSS = gID_field_values[gID]['TSS']
    TTS = gID_field_values[gID]['TTS']
    strand = gID_field_values[gID]['strand']
    #interval = (TSS-500, TSS+500)
    if strand == '+':
        #interval = (TSS, TTS)
        interval = (TSS-500, TSS+500)
    else:
        #interval = (TTS, TSS)
        interval = (TSS-500, TSS+500)
    biotype = gID_biotype[gID]
    if biotype not in selected_biotypes:
        continue
    gID_ginterval[gID + ':' + biotype] = interval

ginterval_dict = Interval_dict.double_hash(gID_ginterval, 10000, 250000000)

ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+ "hg19_chr1_NCP_anot.cn")
ID_score1 = name_ID_value["data/sp_spd_tests_detail/sp7"]
ID_score2 = name_ID_value["data/sp_spd_tests_detail/sp8"]

biotype_IDs = {}

for ID in ID_pos:
    pos = ID_pos[ID]
    find_gIDs = ginterval_dict.find(pos)
    if len(find_gIDs) <= 0:
        if 'Intergenic' not in biotype_IDs:
            biotype_IDs['Intergenic'] = []
        biotype_IDs['Intergenic'].append(ID)
        continue
    for gID in find_gIDs:
        biotype = gID.split(':')[1]
        if biotype not in biotype_IDs:
            biotype_IDs[biotype] = []
        biotype_IDs[biotype].append(ID)

biotype_list = ['protein_coding', 'pseudogene', 'antisense', 'sense_intronic', 'rRNA', 'lincRNA', 'snRNA', 'snoRNA', 'miRNA', 'misc_RNA', 'Intergenic']
scores_list = []

for biotype in biotype_list:
    IDs = biotype_IDs[biotype]
    temp = [ ID_score1[ID] for ID in IDs]
    scores_list.append(temp)

fig = plt.figure()
plt.boxplot(scores_list, 0, "")
plt.title("Condensability by gene type")
#plt.xlabel('non-coding DNA type')
plt.ylabel('Condensability (A.U.)')
plt.xticks(range(1, len(biotype_list)+1), biotype_list, rotation=45)
plt.tight_layout()
plt.savefig("biotype_pbox.png")
#plt.show()
plt.close()

del ginterval_dict

# sort protein-coding genes by gene anatomy (promoter, TSS, 5'UTR, CSS, gene_body(exon, intron), CTS, 3'UTR, TTS)
UTR5len_list = []
UTR3len_list = []
gID_A_interval = {}
for gID in biotype_gIDs['protein_coding']:
    strand = gID_field_values[gID]['strand']
    TSS = gID_field_values[gID]['TSS']
    TTS = gID_field_values[gID]['TTS']
    try:
        CSS = gID_field_values[gID]['CSS']
    except:
        CSS = None
    try:
        CTS = gID_field_values[gID]['CTS']
    except:
        CSS = None
    exons = gID_field_values[gID]['exons']
    introns = [] 

    TSS_interval = (TSS-50, TSS+50)
    gID_A_interval[gID + ":" + "TSS"] = TSS_interval
    if CSS != None:
        CSS_interval = (CSS-50, CSS+50)
        gID_A_interval[gID + ":" + "CSS"] = CSS_interval
    if CTS != None:
        CTS_interval = (CTS-50, CTS+50)
        gID_A_interval[gID + ":" + "CTS"] = CTS_interval
    TTS_interval = (TTS-50, TTS+50)
    gID_A_interval[gID + ":" + "TTS"] = TTS_interval

    if strand == '+':
        Prom_interval = (TSS-500, TSS)
        if CSS != None:
            UTR5_interval = (TSS, CSS)
        else:
            UTR5_interval = None
        if CTS != None:
            UTR3_interval = (CTS, TTS)
        else:
            UTR3_interval = None
        for i in range(len(exons)+1):
            if i == 0:
                intron = [TSS, exons[i][0]]
            elif i == len(exons):
                intron = [exons[i-1][1], TTS]
            else:
                intron = [exons[i-1][1], exons[i][0]]
            introns.append(intron)

    elif strand == '-':
        Prom_interval = (TSS, TSS+500)
        if CSS != None:
            UTR5_interval = (CSS, TSS)
        else:
            UTR5_interval = None
        if CTS != None:
            UTR3_interval = (TTS, CTS)
        else:
            UTR3_interval = None
        for i in range(len(exons)+1):
            if i == 0:
                intron = [TTS, exons[i][0]]
            elif i == len(exons):
                intron = [exons[i-1][1], TSS]
            else:
                intron = [exons[i-1][1], exons[i][0]]
            introns.append(intron)

    gID_A_interval[gID + ":" + "Prom"] = Prom_interval
    if UTR5_interval != None:
        gID_A_interval[gID + ":" + "5UTR"] = UTR5_interval
        UTR5len_list.append(UTR5_interval[1] - UTR5_interval[0])
    if UTR3_interval != None:
        gID_A_interval[gID + ":" + "3UTR"] = UTR3_interval
        UTR3len_list.append(UTR3_interval[1] - UTR3_interval[0])
        
    for i in range(len(exons)):
        gID_A_interval[gID + ":" + "exon" + '_' + str(i)] = tuple(exons[i])
    for i in range(len(introns)):
        gID_A_interval[gID + ":" + "intron" + '_' + str(i)] = tuple(introns[i])


GAinterval_dict = Interval_dict.double_hash(gID_A_interval, 10000, 250000000)

anatomy_IDs = {}

for ID in ID_pos:
    pos = ID_pos[ID]
    find_gIDs = GAinterval_dict.find(pos)
    for gID in find_gIDs:
        anatomy = gID.split(':')[1].split('_')[0]
        if anatomy not in anatomy_IDs:
            anatomy_IDs[anatomy] = []
        anatomy_IDs[anatomy].append(ID)


anatomy_list = ['Prom', 'TSS', '5UTR', 'CSS', 'exon', 'intron', 'CTS', '3UTR', 'TTS']
scores_list = []

for anatomy in anatomy_list:
    IDs = anatomy_IDs[anatomy]
    temp = [ID_score1[ID] for ID in IDs]
    scores_list.append(temp)

fig = plt.figure()
plt.boxplot(scores_list, 0, "")
plt.title("Condensability by anatomy of protein-coding genes")
#plt.xlabel('Anatomy')
plt.ylabel('Condensability (A.U.)')
plt.xticks(range(1, len(anatomy_list)+1), anatomy_list, rotation=45)
plt.tight_layout()
plt.savefig("anatomy_pbox.png")
#plt.show()
plt.close()
