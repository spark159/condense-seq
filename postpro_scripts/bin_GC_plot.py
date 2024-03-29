import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
#import seaborn as sns
import copy
import math
import scipy
import seaborn as sns
from scipy import signal

def total (chr_binID_count):
    total = 0
    for chr in chr_binID_count:
        total += sum(chr_binID_count[chr])
    return total

def GC_content(seq):
    seq = seq.upper()
    output = 0.0
    for nt in seq:
        if nt in "GC":
            output += 1.0
        elif nt in 'N':
            output += 0.5
    return output/len(seq)

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


#fname = "DNA_Spermine(4+)_1kb"
#names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile("/home/spark159/Downloads/" + fname + "_bin.cn")

#fname = "H1_NCP_sp_10kb_bin.cn"
#names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

#fname = "H1_NCP-new_spd_10kb_bin.cn"
#names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

#fname = "H1_NCP_Ca_10kb_bin.cn"
#names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

#fname = "H1_NCP_sp_1kb_tlen_bin.cn"
#names, chr_binID_counts, chr_binID_range, chr_binID_GC, chr_binID_tlen = read_bintlenfile(fname)

#fname = "H1_NCP-new_spd_10kb_bin.cn"
#names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

# proteins
"""
fname = "H1_DNA_HP1a_10kb_bin.cn"
names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

fname = "H1_NCP_HP1a_10kb_bin.cn"
names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

fname = "H1_DNA_HP1bSUV_10kb_bin.cn"
names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

fname = "H1_NCP_HP1bSUV_10kb_bin.cn"
names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

fname = "H1_DNA_LKH_10kb_bin.cn"
names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

fname = "H1_NCP_LKH_10kb_bin.cn"
names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

fname = "H1_DNA_Ki67_10kb_bin.cn"
names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

fname = "H1_NCP_Ki67_10kb_bin.cn"
names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

fname = "H1_DNA_FUS_10kb_bin.cn"
names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

fname = "H1_NCP_FUS_10kb_bin.cn"
names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)
"""

#path = "H1_protein_qc_again/"
#fname = path + "HP1_qc_again_100kb_bin.cn"
#names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

#path = "H1_protein_qc_again/"
#fname = path + "Ki67_qc_again_100kb_bin.cn"
#names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

#path = "protein_qc/"
#fname = path + "FUS_qc_100kb_bin.cn"
#names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

#path = "protein_qc/"
#fname = path + "oldLKH_qc_100kb_bin.cn"
#names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

# protein filling
#path = ""
#fname = "H1_new-NCP_HP1a_10kb_bin.cn"
#fname = "H1_new-NCP_LKH_10kb_bin.cn"
#fname = "H1_new-NCP_Ki67_10kb_bin.cn"
#fname = "H1_NCP_FUS_10kb_bin.cn"
#names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)

# GM NCP
#path=""
#fname = "GM_NCP_sp_10kb_bin.cn"
#fname = "GM_NCP_spd_10kb_bin.cn"
#fname = "GM_NCP_CoH_10kb_bin.cn"
#fname = "GM_NCP_PEG_10kb_bin.cn"

# Progeria NCP
#path = '/media/spark159/sw/'
#fname = 'HGPS_NCP_sp_bin.cn'

# mouse CD8 T cell
#path= ""
#path = '/media/spark159/sw/'
#fname = "mCD8T_WT-NCP_sp_10kb_bin.cn"
#fname = "mCD8T_inht-NCP_sp_10kb_bin.cn"
#fname = 'mCD8T_KO-NCP_sp_bin.cn'

# some replicates QC
path = "/home/spark159/../../storage/replicates/"
tnums = range(1, 10)
fname = 'H1_NCP_sp_10kb_bin.cn'
fname = 'GM_NCP_sp_10kb_bin.cn'
fname = 'mCD8T_WT-NCP_sp_10kb_bin.cn'
fname = 'mCD8T_inht-NCP_sp_10kb_bin.cn'
fname = 'mCD8T_KO-NCP_sp_10kb_bin.cn'
fname = 'H1_NCP_spd_10kb_bin.cn'
fname = 'H1_NCP_CoH_10kb_bin.cn'
fname = 'H1_NCP_PEG_10kb_bin.cn'
fname = 'H1_NCP_Ca_10kb_bin.cn'
tnums = range(1, 6)
#fname = 'H1_NCP_HP1a_10kb_bin.cn'
#fname = 'H1_NCP_LKH_10kb_bin.cn'
#fname = 'H1_NCP_Ki67_10kb_bin.cn'
#fname = 'H1_DNA_HP1a_10kb_bin.cn'
#fname = 'H1_DNA_LKH_10kb_bin.cn'
#fname = 'H1_DNA_Ki67_10kb_bin.cn'
fname = 'H1_NCP_HP1bSUV_10kb_2_bin.cn'
fname = 'H1_NCP_HP1bSUV_10kb_3_bin.cn'
fname = 'H1_NCP_HP1bTRIM_10kb_1_bin.cn'
fname = 'H1_NCP_HP1bTRIM_10kb_2_bin.cn'
fname = 'H1_DNA_HP1bSUV_10kb_2_bin.cn'
fname = 'H1_DNA_HP1bTRIM_10kb_1_bin.cn'
fname = 'H1_DNA_HP1bTRIM_10kb_2_bin.cn'
fname = 'H1_DNA_HP1a_10kb_2_bin.cn'
fname = 'H1_DNA_HP1a_10kb_3_bin.cn'
fname = 'H1_NCP_PEG_10kb_2_bin.cn'
fname = 'H1_NCP_PEG_10kb_3_bin.cn'
fname = 'mCD8T_KO-NCP_sp_10kb_bin.cn'






#names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(fname)
names, chr_binID_counts, chr_binID_range, chr_binID_GC, chr_binID_tlen = read_bintlenfile(path+fname)



X = []
Ys = [[] for i in range(len(names)-1)]
GC_rcounts = [ {} for i in range(len(names)-1) ]
control_total = float(total(chr_binID_counts[-1]))
totals = [ float(total(chr_binID_counts[i])) for i in range(len(names)-1) ]

chr_binID_control = chr_binID_counts[-1]
for chr in chr_binID_control:
    for binID in range(len(chr_binID_control[chr])):
        control_count = chr_binID_control[chr][binID] / control_total
        #control_count = chr_binID_control[chr][binID]
        if control_count <= 0:
            continue
        #control_count += 1
        GC = chr_binID_GC[chr][binID]
        #GC = 1.0 - chr_binID_GC[chr][binID]
        X.append(GC*100)
        
        for i in range(len(names)-1):
            #test_count = chr_binID_counts[i][chr][binID]
            test_count = chr_binID_counts[i][chr][binID] / totals[i]
            #test_count = chr_binID_counts[i][chr][binID] + 1
            rcount = float(test_count)/control_count
            #rcount = float(test_count) / 10.0
            #rcount = -np.log(float(test_count)/control_count)
            Ys[i].append(rcount)

            GC_rcount = GC_rcounts[i]
            if GC not in GC_rcount:
                GC_rcount[GC] = []
            GC_rcount[GC].append(rcount)


# GC content vs normalized counts
for i in range(len(Ys)):
    try:
        name = 'Titration point ' + str(tnums[i])
    except:
        name = 'Titration point ' + str(i+1)
    GC_rcount = GC_rcounts[i]
    Y = Ys[i]
    fig = plt.figure()
    #plt.plot(X, Y, '.', markersize=1, alpha=0.02, color='blue')
    plt.plot(X, Y, '.', markersize=1, alpha=0.05, color='blue')
    corr = scipy.stats.pearsonr(X, Y)[0]
    print 'GC content VS Norm counts', corr
    #print 'AT content VS Condensability', corr
    GCs = []
    means, stds = [], []
    for GC in GC_rcount:
        rcount = GC_rcount[GC]
        #if len(rcount) <= 1:
        #    continue
        #GCs.append(GC*100)
        GCs.append(GC*100)
        means.append(np.mean(rcount))
        stds.append(np.std(rcount))
    #plt.errorbar(GCs, means, yerr=stds, fmt='.', color='blue')
    plt.plot(GCs, means, 'r.', markersize=2, alpha=0.5, zorder=10)
    plt.title(name)
    plt.xlabel("GC content (%)")
    #plt.xlabel("AT content (%)")
    plt.ylabel("Normalized Counts")
    #plt.ylabel("Condensability")
    plt.ylim([-0.5, 4.5])
    #plt.ylim([-5, 15])
    plt.savefig(fname + "_" + str(i+1) +".png")
    #plt.show()
    plt.close()




"""
X = []
Ys = [[] for i in range(len(names)-1)]
GC_rcounts = [ {} for i in range(len(names)-1) ]
GC_tlen = {}
Z = []
tlen_rcounts = [ {} for i in range(len(names)-1) ]
totals = [ total(chr_binID_counts[i]) for i in range(len(names)-1) ]

chr_binID_control = chr_binID_counts[-1]
for chr in chr_binID_control:
    for binID in range(len(chr_binID_control[chr])):
        control_count = chr_binID_control[chr][binID]
        #control_count = chr_binID_control[chr][binID]
        if control_count <= 0:
            continue
        #control_count += 1
        GC = chr_binID_GC[chr][binID]
        #GC = 1.0 - chr_binID_GC[chr][binID]
        X.append(GC*100)
        tlen = chr_binID_tlen[chr][binID]
        Z.append(tlen)

        if GC not in GC_tlen:
            GC_tlen[GC] = []
        GC_tlen[GC].append(tlen)
        
        for i in range(len(names)-1):
            test_count = chr_binID_counts[i][chr][binID]
            #test_count = chr_binID_counts[i][chr][binID] + 1
            rcount = float(test_count)/control_count
            #rcount = float(test_count) / 10.0
            #rcount = -np.log(float(test_count)/control_count)
            Ys[i].append(rcount)

            GC_rcount = GC_rcounts[i]
            if GC not in GC_rcount:
                GC_rcount[GC] = []
            GC_rcount[GC].append(rcount)

            tlen_rcount = tlen_rcounts[i]
            if tlen not in tlen_rcount:
                tlen_rcount[tlen] = []
            tlen_rcount[tlen].append(rcount)


# GC content vs normalized counts
for i in range(len(Ys)):
    name = 'Titration point ' + str(i+1)
    GC_rcount = GC_rcounts[i]
    Y = Ys[i]
    fig = plt.figure()
    plt.plot(X, Y, '.', markersize=1, alpha=0.02, color='blue')
    corr = scipy.stats.pearsonr(X, Y)[0]
    print 'GC content VS Norm counts', corr
    #print 'AT content VS Condensability', corr
    GCs = []
    means, stds = [], []
    for GC in GC_rcount:
        rcount = GC_rcount[GC]
        #if len(rcount) <= 1:
        #    continue
        #GCs.append(GC*100)
        GCs.append(GC*100)
        means.append(np.mean(rcount))
        stds.append(np.std(rcount))
    #plt.errorbar(GCs, means, yerr=stds, fmt='.', color='blue')
    plt.plot(GCs, means, 'r.', markersize=2, alpha=0.5, zorder=10)
    plt.title(name)
    plt.xlabel("GC content (%)")
    #plt.xlabel("AT content (%)")
    plt.ylabel("Normalized Counts")
    #plt.ylabel("Condensability")
    plt.ylim([-0.5, 4])
    #plt.ylim([-5, 15])
    #plt.savefig(fname + "_" + str(i+1) +".png")
    plt.show()
    plt.close()


# tlen vs normalized counts
for i in range(len(Ys)):
    name = 'Titration point ' + str(i+1)
    tlen_rcount = tlen_rcounts[i]
    Y = Ys[i]
    fig = plt.figure()
    plt.plot(Z, Y, '.', markersize=1, alpha=0.01, color='blue')
    corr = scipy.stats.pearsonr(Z, Y)[0]
    print 'read length VS Norm counts', corr
    #print 'Read length VS Condensability', corr
    #corr = scipy.stats.pearsonr(X, Y)[0]
    #print corr
    #data = [[Z[k], Y[k]] for k in range(len(Z))]
    #sns.kdeplot(data)
    #plt.show()
    #plt.close()
    tlens = []
    means, stds = [], []
    for tlen in tlen_rcount:
        rcount = tlen_rcount[tlen]
        #if len(rcount) <= 1:
        #    continue
        tlens.append(tlen)
        means.append(np.mean(rcount))
        stds.append(np.std(rcount))
    #plt.errorbar(tlens, means, yerr=stds, fmt='.', color='blue')
    plt.plot(tlens, means, 'r.', markersize=1, alpha=0.03, zorder=10)
    plt.title(name)
    plt.xlabel("Read length(bp)")
    plt.ylabel("Normalized Counts")
    #plt.ylabel("Condensability")
    #plt.ylim([-0.5, 4])
    #plt.ylim([-5, 15])
    #plt.savefig(fname + "_" + str(i+1) +".png")
    plt.show()
    plt.close()

# GC content vs tlen
fig = plt.figure()
plt.plot(X, Z, '.', markersize=1, alpha=0.02, color='blue')
corr = scipy.stats.pearsonr(X, Z)[0]
print 'GC content VS read length', corr
GCs = []
means, stds = [], []
for GC in GC_tlen:
    GCs.append(100*GC)
    means.append(np.mean(GC_tlen[GC]))
    stds.append(np.std(GC_tlen[GC]))
#plt.errorbar(tlens, means, yerr=stds, fmt='.', color='blue')
plt.plot(GCs, means, 'r.', markersize=2, alpha=0.5, zorder=10)
#print 'AT content VS Read length', corr
plt.xlabel("GC content")
#plt.xlabel("AT content(%)")
plt.ylabel("Read length(bp)")
plt.show()
plt.close()
"""




"""
X = []
Ys = [ [] for i in range(len(names)-1)]
chr_binID_control = chr_binID_counts[-1]
for chr in chr_binID_control:
    for binID in range(len(chr_binID_control[chr])):
        control_count = chr_binID_control[chr][binID]
        X.append(control_count)
        for i in range(len(names)-1):
            test_count = chr_binID_counts[i][chr][binID] 
            Ys[i].append(test_count)

for i in range(len(Ys)):
    name = 'Titration point ' + str(i+1)
    Y = Ys[i]
    fig = plt.figure()
    plt.plot(X, Y, '.', alpha=0.5)
    plt.xlabel("Control read count")
    plt.ylabel(name + " read count")
    plt.title(name)
    plt.xlim([-500, 1000])
    plt.ylim([-500, 1000])
    plt.savefig(fname + "_" + str(i+1) +".png")
    #plt.show()
    plt.close()
"""
