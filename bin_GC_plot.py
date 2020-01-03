import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import seaborn as sns
import copy
import math
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

fname = "DNA_Spermine(4+)_1kb"
names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile("/home/spark159/Downloads/" + fname + "_bin.cn")

X = []
Ys = [[] for i in range(len(names)-1)]
GC_rcounts = [ {} for i in range(len(names)-1)]
totals = [ total(chr_binID_counts[i]) for i in range(len(names)-1) ]

chr_binID_control = chr_binID_counts[-1]
for chr in chr_binID_control:
    for binID in range(len(chr_binID_control[chr])):
        control_count = chr_binID_control[chr][binID]
        if control_count <= 0:
            continue
        GC = chr_binID_GC[chr][binID]
        X.append(GC*100)
        for i in range(len(names)-1):
            test_count = chr_binID_counts[i][chr][binID] 
            rcount = float(test_count)/control_count
            #rcount = float(test_count) / 10.0
            Ys[i].append(rcount)
            GC_rcount = GC_rcounts[i]
            if GC not in GC_rcount:
                GC_rcount[GC] = []
            GC_rcount[GC].append(rcount)

for i in range(len(Ys)):
    name = 'Titration point ' + str(i+1)
    GC_rcount = GC_rcounts[i]
    Y = Ys[i]
    fig = plt.figure()
    plt.plot(X, Y, '.', markersize=2, alpha=0.1, color='blue')
    GCs = []
    means, stds = [], []
    for GC in GC_rcount:
        rcount = GC_rcount[GC]
        #if len(rcount) <= 1:
        #    continue
        GCs.append(GC*100)
        means.append(np.mean(rcount))
        stds.append(np.std(rcount))
    #plt.errorbar(GCs, means, yerr=stds, fmt='.', color='blue')
    plt.plot(GCs, means, 'r.', markersize=3, zorder=10)
    plt.title(name)
    plt.xlabel("GC content (%)")
    plt.ylabel("Normalized Counts")
    #plt.ylim([-0.5, 4])
    #plt.ylim([-5, 15])
    #plt.savefig(fname + "_" + str(i+1) +".png")
    plt.show()
    plt.close()




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
