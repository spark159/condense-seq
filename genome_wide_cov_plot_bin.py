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

def read_bincountfile (fname, chr_list=None):
    First = True
    for line in open(fname):
        if First:
            cols = line.strip().split()
            names = [name.rsplit('.')[-2] for name in cols[4:-1]]
            #names = [name.rsplit('.')[-2] for name in cols[4:-2]]
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
        st, ed = int(st), int(ed)
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

def read_bintlenfile (fname, chr_list=None):
    First = True
    for line in open(fname):
        if First:
            cols = line.strip().split()
            names = [name.rsplit('.')[-2] for name in cols[4:-2]]
            chr_binID_counts = [{} for i in range(len(names))]
            chr_binID_range = {}
            chr_binID_GC = {}
            chr_binID_tlen = {}
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
        st, ed = int(st), int(ed)
        GC = float(cols[-2])
        tlen = float(cols[-1])
        if chr not in chr_binID_range:
            chr_binID_range[chr] = []
        chr_binID_range[chr].append((st, ed))
        if chr not in chr_binID_GC:
            chr_binID_GC[chr] = []
        chr_binID_GC[chr].append(GC)
        if chr not in chr_binID_tlen:
            chr_binID_tlen[chr] = []
        chr_binID_tlen[chr].append(tlen)
        datas = cols[4:-2]
        for i in range(len(datas)):
            data = float(datas[i])
            chr_binID_count = chr_binID_counts[i]
            if chr not in chr_binID_count:
                chr_binID_count[chr] = []
            chr_binID_count[chr].append(data)
    return names, chr_binID_counts, chr_binID_range, chr_binID_GC, chr_binID_tlen


path = ""
#fname = "H1_NCP-new_spd_10kb_bin.cn"
#fname = "NCP_Spermine(4+)_10kb_bin.cn"

#fname = "H1_DNA_sp_10kb_bin.cn"
#fname = "H1_NCP_sp_10kb_bin.cn"
#fname = "H1_DNA_spd_10kb_bin.cn"
#fname = "H1_NCP_spd_10kb_bin.cn"
#fname = "H1_DNA_CoH_10kb_bin.cn"
#fname = "H1_NCP_CoH_10kb_bin.cn
#fname = "H1_DNA_PEG_10kb_bin.cn"
#fname = "H1_NCP_PEG_10kb_bin.cn"
#fname = "H1_NCP_Mg_10kb_bin.cn"
#fname = "H1_NCP_Ca_10kb_bin.cn"

#fname = "H1_NCP_sp_1kb_bin.cn"
#fname = "H1_NCP_sp_1kb_tlen_bin.cn"

# proteins
#fname = "H1_DNA_HP1a_10kb_bin.cn"
#fname = "H1_NCP_HP1a_10kb_bin.cn"
#fname = "H1_DNA_HP1bSUV_10kb_bin.cn"
#fname = "H1_NCP_HP1bSUV_10kb_bin.cn"
#fname = "H1_DNA_LKH_10kb_bin.cn"
#fname = "H1_NCP_LKH_10kb_bin.cn"
#fname = "H1_DNA_Ki67_10kb_bin.cn"
#fname = "H1_NCP_Ki67_10kb_bin.cn"
#fname = "H1_DNA_FUS_10kb_bin.cn"
#fname = "H1_NCP_FUS_10kb_bin.cn"

#path = "../H1_protein_qc_again/"

# protein filling
#path = ""
#fname = "H1_new-NCP_HP1a_10kb_bin.cn"
#fname = "H1_new-NCP_LKH_10kb_bin.cn"
#fname = "H1_new-NCP_Ki67_10kb_bin.cn"
#fname = "H1_NCP_FUS_10kb_bin.cn"

# GM NCP
#path=""
#fname = "GM_NCP_sp_10kb_bin.cn"
#fname = "GM_NCP_spd_10kb_bin.cn"
#fname = "GM_NCP_CoH_10kb_bin.cn"
#fname = "GM_NCP_PEG_10kb_bin.cn"

# Progeria NCP
path = '/media/spark159/sw/'
fname = 'HGPS_NCP_sp_bin.cn'

# mouse CD8 T cell
#path= ""
#path = '/media/spark159/sw/'
#fname = "mCD8T_WT-NCP_sp_10kb_bin.cn"
#fname = "mCD8T_inht-NCP_sp_10kb_bin.cn"
#fname = 'mCD8T_KO-NCP_sp_bin.cn'

# replicates QC
path = "/home/spark159/../../storage/replicates/"
#path = "/home/spark159/../../storage/"

#fname = 'H1_NCP_sp_10kb_bin.cn'
#fname = 'GM_NCP_sp_10kb_bin.cn'
#fname = 'mCD8T_WT-NCP_sp_10kb_bin.cn'
#fname = 'mCD8T_inht-NCP_sp_10kb_bin.cn'
#fname = 'mCD8T_KO-NCP_sp_10kb_bin.cn'
#fname = 'H1_NCP_HP1a_10kb_bin.cn'
#fname = 'H1_NCP_LKH_10kb_bin.cn'
#fname = 'H1_NCP_Ki67_10kb_bin.cn'
#fname = 'H1_NCP_spd_10kb_bin.cn'
#fname = 'H1_NCP_CoH_10kb_bin.cn'
#fname = 'H1_NCP_PEG_10kb_bin.cn'
#fname = 'H1_NCP_Ca_10kb_bin.cn'
#fname = 'H1_DNA_HP1a_10kb_bin.cn'
#fname = 'H1_DNA_LKH_10kb_bin.cn'
#fname = 'H1_DNA_Ki67_10kb_bin.cn'
#fname = 'H1_NCP_HP1bSUV_10kb_2_bin.cn'
#fname = 'H1_NCP_HP1bSUV_10kb_3_bin.cn'
#fname = 'H1_NCP_HP1bTRIM_10kb_1_bin.cn'
#fname = 'H1_NCP_HP1bTRIM_10kb_2_bin.cn'
#fname = 'H1_DNA_HP1bSUV_10kb_2_bin.cn'
#fname = 'H1_DNA_HP1bTRIM_10kb_1_bin.cn'
#fname = 'H1_DNA_HP1bTRIM_10kb_2_bin.cn'
fname = 'H1_DNA_HP1a_10kb_2_bin.cn'
fname = 'H1_DNA_HP1a_10kb_3_bin.cn'
fname = 'H1_NCP_PEG_10kb_2_bin.cn'
fname = 'H1_NCP_PEG_10kb_3_bin.cn'
fname = 'mCD8T_KO-NCP_sp_10kb_bin.cn'




#names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(path+fname, chr_list=['chr1'])
names, chr_binID_counts, chr_binID_range, chr_binID_GC, chr_binID_tlen = read_bintlenfile(path+fname, chr_list=['chr1'])

#chr_binID_control = chr_binID_counts[-1]

for k in range(len(chr_binID_counts)):
    chr_binID_count = chr_binID_counts[k]
    ID_pos = {}
    ID_count = {}
    for i in range(len(chr_binID_count['chr1'])):
        start, end = chr_binID_range['chr1'][i]
        pos = (start+end)/2
        ID_pos[i] = pos
        ID_count[i] = chr_binID_count['chr1'][i]
    mean = np.median(ID_count.values())
    std = np.std(ID_count.values())
    if k == len(chr_binID_counts) - 1:
        title='Chromosome1 t#0'
        note = "_" + fname.rsplit('.',1)[0]+ '_' + str(0)
    else:
        title='Chromosome1 t#' + str(k+1)
        note = "_" + fname.rsplit('.',1)[0]+ '_' + str(k+1)
        
    #graphics.draw_along_genome (ID_pos, [ID_count], win=1000, labels=['Coverage'], ylabel='Read counts per 10kb', ylim=[max(mean-3*std, -5), mean+2*std], title=title, scatt=True, note=note)
    graphics.draw_along_genome (ID_pos, [ID_count], win=1000, labels=['Coverage'], ylabel='Read counts per 10kb', ylim=[max(mean-5*std, -5), mean+4*std], title=title, scatt=True, note=note)




#for k in range(12):
#    fname = "%d_100kb_bin.cn" % (k+1)
#    names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(path+fname, chr_list=['chr1'])
#    chr_binID_control = chr_binID_counts[0]
#    ID_pos = {}
#    ID_count = {}
#    for i in range(len(chr_binID_control['chr1'])):
#        start, end = chr_binID_range['chr1'][i]
#        pos = (start+end)/2
#        ID_pos[i] = pos
#        ID_count[i] = chr_binID_control['chr1'][i]
#    mean = np.median(ID_count.values())
#    std = np.std(ID_count.values())
#    graphics.draw_along_genome (ID_pos, [ID_count], win=100, labels=['Input'], ylabel='Read counts per 100kb', ylim=[max(mean-3*std, -5), mean+2*std], title='Chromosome 1', scatt=True, note="_" + fname.rsplit('.',1)[0])



#names, chr_binID_counts, chr_binID_range, chr_binID_GC, chr_binID_tlen = read_bintlenfile(fname, chr_list=['chr1'])
#chr_binID_control = chr_binID_counts[5]

#for k in range(6):
#    chr_binID_control = chr_binID_counts[k]
#    ID_pos = {}
#    ID_count = {}
#    for i in range(len(chr_binID_control['chr1'])):
#        start, end = chr_binID_range['chr1'][i]
#        pos = (start+end)/2
#        ID_pos[i] = pos
#        ID_count[i] = chr_binID_control['chr1'][i]
#    graphics.draw_along_genome (ID_pos, [ID_count], win=1000, labels=['Input'], ylabel='Read counts per 10kb', ylim=[-5,160], title='Chromosome 1', scatt=True, note="_" + fname.rsplit('.',1)[0] + '_' + str(k))

    #ID_counts = chr_binID_control['chr1']
    #graphics.draw_along_genome (ID_pos, [ID_count], win=1000, labels=['Input'], ylabel='Read counts per 10kb', ylim=[-5,180], title='Chromosome 1', scatt=True, note="_" + fname.rsplit('.',1)[0])
    #graphics.draw_along_genome (ID_pos, [ID_count], win=1000, labels=['Input'], ylabel='Read counts per 1kb',ylim=[-5, 200],  title='Chromosome 1', scatt=True, note="_" + fname.rsplit('.',1)[0])

    #ID_tlen = chr_binID_tlen['chr1']
    #graphics.draw_along_genome (ID_pos, [ID_tlen], win=1000, labels=['Input'], ylabel='Read length', ylim=[100,250], title='Chromosome 1', scatt=True, note="_" + fname.rsplit('.',1)[0])
