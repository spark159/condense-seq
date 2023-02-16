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
import matplotlib as mpl


def read_rlen(fname, chr_list=None):
    rlen_count = {}
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        line = line.strip()
        if not line:
            continue
        cols = line.strip().split()
        chr, rlen, count = cols
        if chr_list and chr not in chr_list:
            continue
        rlen, count = int(rlen), int(count)
        if rlen not in rlen_count:
            rlen_count[rlen] = 0
        rlen_count[rlen] += count
    return rlen_count

#path = "/home/spark159/../../media/spark159/sw/sp_spd_tests_detail/"
#path = ""
#path = "/home/spark159/../../media/spark159/sw/"

#fnames = ["H1_DNA_spd_H1-DNA-spd-0_rlen.txt", "H1_DNA_sp_H1-DNA-sp-0_rlen.txt", "H1_NCP-new_spd_H1-NCP-new-spd-0_rlen.txt", "H1_NCP-new_sp_H1-NCP-new-sp-0_rlen.txt"]
#names = ["H1-DNA-spd-0", "H1-DNA-sp-0", "H1-NCP-spd-0", "H1-NCP-sp-0"]

#fnames = ["H1_NCP-new_sp_H1-NCP-new-sp-%s_rlen.txt" % (i) for i in range(6)]
#names = ["H1-NCP-sp-%s" % (i) for i in range(6)]

#fnames = ["H1_NCP-new_spd_H1-NCP-new-spd-%s_rlen.txt" % (i) for i in range(6)]
#names = ["H1-NCP-spd-%s" % (i) for i in range(6)] 

#fnames = ["H1_DNA_sp_H1-DNA-sp-%s_rlen.txt" % (i) for i in range(6)]
#names = ["H1-DNA-sp-%s" % (i) for i in range(6)] 

#fnames = ["H1_DNA_spd_H1-DNA-spd-%s_rlen.txt" % (i) for i in range(6)]
#names = ["H1-DNA-sp-%s" % (i) for i in range(6)]

#fnames = ["H1_DNA_sp_H1-DNA-sp-%s_rlen.txt" % (i) for i in range(6)]
#names = ["H1-DNA-sp-%s" % (i) for i in range(6)] 

#fnames = ["H1_NCP_sp_H1-NCP-sp-%s_rlen.txt" % (i) for i in range(10)]
#names = ["H1-NCP-sp-%s" % (i) for i in range(10)] 

#fnames = ["H1_DNA_spd_H1-DNA-spd-%s_rlen.txt" % (i) for i in range(6)]
#names = ["H1-DNA-spd-%s" % (i) for i in range(6)] 

#fnames = ["H1_NCP_spd_H1-NCP-spd-%s_rlen.txt" % (i) for i in range(10)]
#names = ["H1-NCP-spd-%s" % (i) for i in range(10)] 

#fnames = ["H1_DNA_CoH_H1-DNA-CoH-%s_rlen.txt" % (i) for i in range(6)]
#names = ["H1-DNA-CoH-%s" % (i) for i in range(6)] 

#fnames = ["H1_NCP_CoH_H1-NCP-CoH-%s_rlen.txt" % (i) for i in range(8)]
#names = ["H1-NCP-CoH-%s" % (i) for i in range(8)]

#fnames = ["H1_DNA_PEG_H1-DNA-PEG-%s_rlen.txt" % (i) for i in range(7)]
#names = ["H1-DNA-PEG-%s" % (i) for i in range(7)] 

#fnames = ["H1_NCP_PEG_H1-NCP-PEG-%s_rlen.txt" % (i) for i in range(8)]
#names = ["H1-NCP-PEG-%s" % (i) for i in range(8)]

#fnames = ["H1_NCP_Mg_H1-NCP-Mg-%s_rlen.txt" % (i) for i in range(8)]
#names = ["H1-NCP-Mg-%s" % (i) for i in range(8)]

#fnames = ["H1_NCP_Ca_H1-NCP-Ca-%s_rlen.txt" % (i) for i in range(8)]
#names = ["H1-NCP-Ca-%s" % (i) for i in range(8)]

#fnames = ["H1_NCP_Ca_H1-NCP-Ca-%s_rlen.txt" % (i) for i in range(8)]
#names = ["H1-NCP-Ca-%s" % (i) for i in range(8)]

#fnames = ["H1_NCP_sp_H1-NCP-sp-%s_rlen.txt" % (i) for i in [0, 4, 8]]
#names = ["H1-NCP-sp-%s" % (i) for i in [0, 4, 8]]


# proteins
"""
fnames = ["H1_DNA_HP1a_H1-DNA-HP1a-%s_rlen.txt" % (i) for i in range(6)]
names = ["H1-DNA-HP1a-%s" % (i) for i in range(6)]

fnames = ["H1_NCP_HP1a_H1-NCP-HP1a-%s_rlen.txt" % (i) for i in range(6)]
names = ["H1-NCP-HP1a-%s" % (i) for i in range(6)]

fnames = ["H1_DNA_HP1bSUV_H1-DNA-HP1bSUV-%s_rlen.txt" % (i) for i in range(6)]
names = ["H1-DNA-HP1bSUV-%s" % (i) for i in range(6)]

fnames = ["H1_NCP_HP1bSUV_H1-NCP-HP1bSUV-%s_rlen.txt" % (i) for i in range(6)]
names = ["H1-NCP-HP1bSUV-%s" % (i) for i in range(6)]

fnames = ["H1_DNA_LKH_H1-DNA-LKH-%s_rlen.txt" % (i) for i in range(6)]
names = ["H1-DNA-LKH-%s" % (i) for i in range(6)]

fnames = ["H1_NCP_LKH_H1-NCP-LKH-%s_rlen.txt" % (i) for i in range(6)]
names = ["H1-NCP-LKH-%s" % (i) for i in range(6)]

fnames = ["H1_DNA_Ki67_H1-DNA-Ki67-%s_rlen.txt" % (i) for i in range(6)]
names = ["H1-DNA-Ki67-%s" % (i) for i in range(6)]

fnames = ["H1_NCP_Ki67_H1-NCP-Ki67-%s_rlen.txt" % (i) for i in range(6)]
names = ["H1-NCP-Ki67-%s" % (i) for i in range(6)]

fnames = ["H1_DNA_FUS_H1-DNA-FUS-%s_rlen.txt" % (i) for i in range(6)]
names = ["H1-DNA-FUS-%s" % (i) for i in range(6)]

fnames = ["H1_NCP_FUS_H1-NCP-FUS-%s_rlen.txt" % (i) for i in range(6)]
names = ["H1-NCP-FUS-%s" % (i) for i in range(6)]
"""

#path = "protein_qc/"

#fnames = ["%d_%d_rlen.txt" % (i, i) for i in range(1, 5)]
#names = ["H1-NCP-HP1a-qc-%s" % (i) for i in range(4)]
#note = 'HP1a_qc'

#fnames = ["%d_%d_rlen.txt" % (i, i) for i in range(5, 11)]
#names = ["H1-NCP-LKH-qc-%s" % (i) for i in range(6)]
#note = 'LKH_qc'

#fnames = ["%d_%d_rlen.txt" % (i, i) for i in range(11, 17)]
#names = ["H1-NCP-Ki67-qc-%s" % (i) for i in range(6)]
#note = 'Ki67_qc'

#fnames = ["%d_%d_rlen.txt" % (i, i) for i in range(17, 23)]
#names = ["H1-NCP-FUS-qc-%s" % (i) for i in range(6)]
#note = 'FUS_qc'

#fnames = ["%d_%d_rlen.txt" % (i, i) for i in range(23, 29)]
#names = ["H1-NCP-old-LKH-qc-%s" % (i) for i in range(6)]
#note = 'oldKLH_qc'

#path = "../H1_protein_qc_again/"

#fnames = ["%d_%d_rlen.txt" % (i, i) for i in range(1, 7)]
#names = ["H1-NCP-HP1a-qc-%s" % (i) for i in range(6)]
#note = 'HP1a_qc'

#fnames = ["%d_%d_rlen.txt" % (i, i) for i in range(7, 13)]
#names = ["H1-NCP-Ki67-qc-%s" % (i) for i in range(6)]
#note = 'Ki67_qc'


# protein fillings
#fnames = ["H1_new-NCP_HP1a_H1-new-NCP-HP1a-%s_rlen.txt" % (i) for i in range(6)]
#names = ["H1-NCP-HP1a-%s" % (i) for i in range(6)]

#fnames = ["H1_new-NCP_LKH_H1-new-NCP-LKH-%s_rlen.txt" % (i) for i in range(6)]
#names = ["H1-NCP-LKH-%s" % (i) for i in range(6)]

#fnames = ["H1_new-NCP_Ki67_H1-new-NCP-Ki67-%s_rlen.txt" % (i) for i in range(6)]
#names = ["H1-NCP-Ki67-%s" % (i) for i in range(6)]

#fnames = ["H1_NCP_FUS_H1-NCP-FUS-%s_rlen.txt" % (i) for i in range(6)]
#names = ["H1-NCP-FUS-%s" % (i) for i in range(6)]

#note='filling'

# GM cells 
#fnames = ["GM_NCP_sp_GM-NCP-sp-%s_rlen.txt" % (i) for i in range(10)]
#names = ["GM-NCP-sp-%s" % (i) for i in range(10)]
#note="GM_NCP_sp"

#fnames = ["GM_NCP_spd_GM-NCP-spd-%s_rlen.txt" % (i) for i in range(10)]
#names = ["GM-NCP-spd-%s" % (i) for i in range(10)]
#note="GM_NCP_spd"

#fnames = ["GM_NCP_CoH_GM-NCP-CoH-%s_rlen.txt" % (i) for i in range(8)]
#names = ["GM-NCP-CoH-%s" % (i) for i in range(8)]
#note="GM_NCP_CoH"

#fnames = ["GM_NCP_PEG_GM-NCP-PEG-%s_rlen.txt" % (i) for i in range(8)]
#names = ["GM-NCP-PEG-%s" % (i) for i in range(8)]
#note="GM_NCP_PEG"

# Progeria cells

#fnames = ["HGPS_NCP_sp_HGPS-NCP-sp-%s_rlen.txt" % (i) for i in range(10)]
#names = ["HGPS-NCP-sp-%s" % (i) for i in range(10)]
#note="HGPS_NCP_sp"

# mouse CD8 T cells
#fnames = ["mCD8T_WT-NCP_sp_mCD8T-WT-NCP-sp-%s_rlen.txt" % (i) for i in range(10)]
#names = ["mCD8T:WT-NCP-sp-%s" % (i) for i in range(10)]
#note="mCD8T:WT_NCP_sp"

#fnames = ["mCD8T_inht-NCP_sp_mCD8T-inht-NCP-sp-%s_rlen.txt" % (i) for i in range(10)]
#names = ["mCD8T:inht-NCP-sp-%s" % (i) for i in range(10)]
#note="mCD8T:inht_NCP_sp"

#fnames = ["mCD8T_KO-NCP_sp_mCD8T-KO-NCP-sp-%s_rlen.txt" % (i) for i in range(10)]
#names = ["mCD8T:KO-NCP-sp-%s" % (i) for i in range(10)]
#note="mCD8T:KO_NCP_sp"


# GM NCP spermine in detail
#fnames = ["GM_NCP_sp_GM-NCP-sp-%s_rlen.txt" % (i) for i in [0, 4, 8]]
#names = ["GM-NCP-sp-%s" % (i) for i in [0, 4, 8]]
#note="GM_NCP_sp_detail"

# H1 NCP HP1a in detail
#fnames = ["H1_new-NCP_HP1a_H1-new-NCP-HP1a-%s_rlen.txt" % (i) for i in [0, 3]]
#names = ["H1-NCP-HP1a-%s" % (i) for i in [0, 3]]
#note="H1_NCP_HP1a_detail"

# H1 DNA HP1a in detail
#fnames = ["H1_DNA_HP1a_H1-DNA-HP1a-%s_rlen.txt" % (i) for i in [0, 3]]
#names = ["H1-DNA-HP1a-%s" % (i) for i in [0, 3]]
#note="H1_DNA_HP1a_detail"


# Mouse CD8 T cell in detail
#fnames = ["mCD8T_WT-NCP_sp_mCD8T-WT-NCP-sp-%s_rlen.txt" % (i) for i in [0, 4, 8]]
#names = ["mCD8T:WT-NCP-sp-%s" % (i) for i in [0, 4, 8]]
#note="mCD8T:WT_NCP_sp"

#fnames = ["mCD8T_inht-NCP_sp_mCD8T-inht-NCP-sp-%s_rlen.txt" % (i) for i in [0, 4, 8]]
#names = ["mCD8T:inht-NCP-sp-%s" % (i) for i in [0, 4, 8]]
#note="mCD8T:inht_NCP_sp"

#fnames = ["mCD8T_KO-NCP_sp_mCD8T-KO-NCP-sp-%s_rlen.txt" % (i) for i in [0, 4, 8]]
#names = ["mCD8T:KO-NCP-sp-%s" % (i) for i in [0, 4, 8]]
#note="mCD8T:KO_NCP_sp"


# replicate QC
path = "/home/spark159/../../storage/replicates/"

fnames = ["H1_NCP_sp_H1-NCP-sp-%s_rlen.txt" % (i) for i in range(10)]
names = ["H1-NCP-sp-%s" % (i) for i in range(10)]
note='H1_NCP_sp_qc'

fnames = ["GM_NCP_sp_GM-NCP-sp-%s_rlen.txt" % (i) for i in range(10)]
names = ["GM-NCP-sp-%s" % (i) for i in range(10)]
note='GM_NCP_sp_qc'

fnames = ["mCD8T_WT-NCP_sp_mCD8T-WT-NCP-sp-%s_rlen.txt" % (i) for i in range(10)]
names = ["mCD8T-WT-NCP-sp-%s" % (i) for i in range(10)]
note='mCD8T_WT-NCP_sp_qc'

fnames = ["mCD8T_inht-NCP_sp_mCD8T-inht-NCP-sp-%s_rlen.txt" % (i) for i in range(10)]
names = ["mCD8T-inht-NCP-sp-%s" % (i) for i in range(10)]
note='mCD8T_inht-NCP_sp_qc'

fnames = ["mCD8T_KO-NCP_sp_mCD8T-KO-NCP-sp-%s_rlen.txt" % (i) for i in range(10)]
names = ["mCD8T-KO-NCP-sp-%s" % (i) for i in range(10)]
note='mCD8T_KO-NCP_sp_qc'

fnames = ["H1_NCP_HP1a_H1-NCP-HP1a-%s_rlen.txt" % (i) for i in range(6)]
names = ["H1-NCP-HP1a-%s" % (i) for i in range(6)]
note='H1_NCP_HP1a_qc'

fnames = ["H1_NCP_LKH_H1-NCP-LKH-%s_rlen.txt" % (i) for i in range(6)]
names = ["H1-NCP-LKH-%s" % (i) for i in range(6)]
note='H1_NCP_LKH_qc'

fnames = ["H1_NCP_Ki67_H1-NCP-Ki67-%s_rlen.txt" % (i) for i in range(6)]
names = ["H1-NCP-Ki67-%s" % (i) for i in range(6)]
note='H1_NCP_Ki67_qc'

fnames = ["H1_NCP_spd_H1-NCP-spd-%s_rlen.txt" % (i) for i in [0, 6]]
names = ["H1-NCP-spd-%s" % (i) for i in [0, 6]]
note='H1_NCP_spd_qc'

fnames = ["H1_NCP_CoH_H1-NCP-CoH-%s_rlen.txt" % (i) for i in [0, 5]]
names = ["H1-NCP-CoH-%s" % (i) for i in [0, 5]]
note='H1_NCP_CoH_qc'

fnames = ["H1_NCP_PEG_H1-NCP-PEG-%s_rlen.txt" % (i) for i in [0, 6]]
names = ["H1-NCP-PEG-%s" % (i) for i in [0, 6]]
note='H1_NCP_PEG_qc'

fnames = ["H1_NCP_Ca_H1-NCP-Ca-%s_rlen.txt" % (i) for i in [0, 5]]
names = ["H1-NCP-Ca-%s" % (i) for i in [0, 5]]
note='H1_NCP_Ca_qc'

fnames = ["H1_DNA_HP1a_H1-DNA-HP1a-%s_rlen.txt" % (i) for i in [0, 3]]
names = ["H1-DNA-HP1a-%s" % (i) for i in [0, 3]]
note='H1_DNA_HP1a_qc'

#fnames = ["H1_DNA_LKH_H1-DNA-LKH-%s_rlen.txt" % (i) for i in [0, 4]]
#names = ["H1-DNA-LKH-%s" % (i) for i in [0, 4]]
#note='H1_DNA_LKH_qc'

fnames = ["H1_DNA_Ki67_H1-DNA-Ki67-%s_rlen.txt" % (i) for i in [0, 4]]
names = ["H1-DNA-Ki67-%s" % (i) for i in [0, 4]]
note='H1_DNA_Ki67_qc'


# replicates fill-in
path = "/home/spark159/../../storage/replicates/"

fnames = ["_H1_NCP_HP1bSUV_%s_2_rlen.txt" % (i) for i in range(6)]
names = ["H1_NCP_HP1bSUV_%s_2" % (i) for i in range(6)]
note='H1_NCP_HP1bSUV_2rep_qc'

fnames = ["_H1_NCP_HP1bSUV_%s_3_rlen.txt" % (i) for i in range(6)]
names = ["H1_NCP_HP1bSUV_%s_3" % (i) for i in range(6)]
note='H1_NCP_HP1bSUV_3rep_qc'

fnames = ["_H1_NCP_HP1bTRIM_%s_1_rlen.txt" % (i) for i in range(6)]
names = ["H1_NCP_HP1bTRIM_%s_1" % (i) for i in range(6)]
note='H1_NCP_HP1bTRIM_1rep_qc'

#fnames = ["_H1_NCP_HP1bTRIM_%s_2_rlen.txt" % (i) for i in range(6)]
#names = ["H1_NCP_HP1bTRIM_%s_2" % (i) for i in range(6)]
#note='H1_NCP_HP1bTRIM_2rep_qc'

fnames = ["_H1_DNA_HP1bSUV_%s_2_rlen.txt" % (i) for i in range(6)]
names = ["H1_DNA_HP1bSUV_%s_2" % (i) for i in range(6)]
note='H1_DNA_HP1bSUV_2rep_qc'

fnames = ["_H1_DNA_HP1bTRIM_%s_1_rlen.txt" % (i) for i in range(6)]
names = ["H1_DNA_HP1bTRIM_%s_1" % (i) for i in range(6)]
note='H1_DNA_HP1bTRIM_1rep_qc'

fnames = ["_H1_DNA_HP1bTRIM_%s_2_rlen.txt" % (i) for i in range(6)]
names = ["H1_DNA_HP1bTRIM_%s_2" % (i) for i in range(6)]
note='H1_DNA_HP1bTRIM_2rep_qc'

fnames = ["_H1_DNA_HP1a_%s_2_rlen.txt" % (i) for i in [0, 1, 2, 4]]
names = ["H1_DNA_HP1a_%s_2" % (i) for i in [0, 1, 2, 4]]
note='H1_DNA_HP1a_2rep_qc'

fnames = ["_H1_DNA_HP1a_%s_3_rlen.txt" % (i) for i in [0, 2, 3, 4]]
names = ["H1_DNA_HP1a_%s_3" % (i) for i in [0, 2, 3, 4]]
note='H1_DNA_HP1a_3rep_qc'

fnames = ["_H1_NCP_PEG_%s_2_rlen.txt" % (i) for i in [0, 4, 5, 7]]
names = ["H1_NCP_PEG_%s_2" % (i) for i in [0, 4, 5, 7]]
note='H1_NCP_PEG_2rep_qc'

fnames = ["_H1_NCP_PEG_%s_3_rlen.txt" % (i) for i in [0, 5, 6]]
names = ["H1_NCP_PEG_%s_3" % (i) for i in [0, 5, 6]]
note='H1_NCP_PEG_3rep_qc'

fnames = ["mCD8T_KO-NCP_sp_mCD8T-KO-NCP-sp-%s_rlen.txt" % (i) for i in range(10)]
names = ["mCD8T-KO-NCP-sp-%s" % (i) for i in range(10)]
note='mCD8T_KO-NCP_sp_qc'













rlen_count_list = []
X_list, Y_list = [], []
for fname in fnames:
    rlen_count = read_rlen(path+fname)
    X = sorted(rlen_count.keys())
    Y = [ rlen_count[x] for x in X ]
    rlen_count_list.append(rlen_count)
    X_list.append(X)
    Y_list.append(Y)

color_list = np.linspace(0.01, 0.7, num=len(names))
alpha_list = np.linspace(0.3, 1, num=len(names))[::-1]
cmap = mpl.cm.get_cmap("copper")

fig = plt.figure()
for i in range(len(names)):
    name = names[i]
    X, Y = X_list[i], Y_list[i]
    #plt.plot(X, Y, label=name, alpha=1)
    plt.plot(X, Y, color=cmap(color_list[i]), label=name, alpha=1)
    #plt.plot(X, Y, color='b', label=name, alpha=alpha_list[i])
plt.xticks(range(0,500,50))
#plt.xticks(range(0,500,10))
plt.grid(True)
plt.xlabel("Read length (bp)")
plt.ylabel("Read counts")
plt.title("Read length distribution")
plt.legend()
plt.savefig("rlen_" + note + ".png")
#plt.savefig(path+"rlen.png")
#plt.savefig(fname.split('.')[0]+".png")
#plt.show()
plt.close()

"""
rlen_count1 = read_rlen(path+"hg19_sp1_rlen.txt")
rlen_count2 = read_rlen(path+"hg19_sp7_rlen.txt")
rlen_count3 = read_rlen(path+"hg19_sp8_rlen.txt")
    
X1 = sorted(rlen_count1.keys())
Y1 = [ rlen_count1[x] for x in X1 ]
X2 = sorted(rlen_count2.keys())
Y2 = [ rlen_count2[x] for x in X2 ]
X3 = sorted(rlen_count3.keys())
Y3 = [ rlen_count3[x] for x in X3 ]

fig = plt.figure()
plt.plot(X1, Y1, label='Input')
plt.plot(X2, Y2, label='NCP Sp6')
plt.plot(X3, Y3, label='NCP Sp7')
plt.xticks(range(0,500,50))
plt.grid(True)
plt.xlabel("Read length (bp)")
plt.ylabel("Read counts")
plt.title("Read length distribution")
plt.legend()
#plt.show()
plt.savefig(path+"rlen.png")
plt.close()
"""
