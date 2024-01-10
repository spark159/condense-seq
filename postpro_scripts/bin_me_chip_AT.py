import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import copy
import math
from pyliftover import LiftOver
from Bio import SeqIO
import Interval_dict
import graphics
import pickle

def read_bincountfile (fname, bin_size, chr_list=None):
    First = True
    for line in open(fname):
        if First:
            cols = line.strip().split()
            names = [name.rsplit('.')[-2] for name in cols[3:-1]]
            chr_binID_counts = [{} for i in range(len(names))]
            chr_binID_range = {}
            chr_binID_GC = {}
            First = False
            continue
        line = line.strip()
        if not line:
            continue
        cols = line.strip().split()
        ID, chr, pos = cols[:3]
        if chr_list != None and chr not in chr_list:
            continue
        ID = int(ID)
        st, end = int(pos)-bin_size/2, int(pos)+bin_size/2
        GC = float(cols[-1])
        if chr not in chr_binID_range:
            chr_binID_range[chr] = []
        chr_binID_range[chr].append((st, end))
        if chr not in chr_binID_GC:
            chr_binID_GC[chr] = []
        chr_binID_GC[chr].append(GC)
        datas = cols[3:-1]
        for i in range(len(datas)):
            data = float(datas[i])
            chr_binID_count = chr_binID_counts[i]
            if chr not in chr_binID_count:
                chr_binID_count[chr] = []
            chr_binID_count[chr].append(data)
    return names, chr_binID_counts, chr_binID_range, chr_binID_GC
path = "/home/spark159/../../media/spark159/sw/sp_spd_tests_detail/"
names, chr_binID_counts, chr_binID_range, chr_binID_GC = read_bincountfile(path + "hg19_1kb_bin.cn", 1001)
chr_binID_control = chr_binID_counts[-1]

def read_BS_file (fname, bin_size, chr_list=None):
    ID_me = {}
    First = True
    line_count = -1
    for line in open(fname):
        line_count +=1
        #print line_count
        cols = line.strip().split()
        if First:
            First = False
            continue
        chr, pos, strand, count, total = cols
        if chr_list != None and chr not in chr_list:
            continue
        pos = int(pos) - 1
        count, total = int(count), int(total)
        if total <= 0:
            sig = 0.5
            #continue
        else:
            sig = float(count) / total
        binID = pos / bin_size
        ID = chr + ':' + str(binID)
        if ID not in ID_me:
            ID_me[ID] = []
        ID_me[ID].append(sig)
    ID_mean = {}
    for ID in ID_me:
        ID_mean[ID] = np.mean(ID_me[ID])
    return ID_mean

try:
    f = open("temp.pickle", "rb")
    ID_me = pickle.load(f)
except:
    ID_me = read_BS_file("/home/spark159/../../media/spark159/sw/dataforcondense/GSM1541790_38-Per_rep1.cg.txt", bin_size=10000)
    f = open("temp.pickle", "wb")
    pickle.dump(ID_me, f)
    f.close()


def read_chip_file (fname, bin_size, chr_list=None):
    bin_score = {}
    for line in open(fname):
        cols = line.strip().split()
        chr, st, ed = cols[0], int(cols[1]), int(cols[2])
        if chr_list != None and chr not in chr_list:
            continue
        score = float(cols[4])
        st_bin = int(st) / int(bin_size)
        ed_bin = int(ed) / int(bin_size)
        for k in range(st_bin, ed_bin+1):
            if k == st_bin:
                width = ((st_bin + 1) * bin_size) - st
            elif k == ed_bin:
                width = ed - ed_bin*bin_size
            else:
                width = bin_size
            ID = chr + ':' + str(k)
            if ID not in bin_score:
                bin_score[ID] = 0.0
            bin_score[ID] += score*width
    return bin_score
bin_k27ac = read_chip_file("/home/spark159/../../media/spark159/sw/dataforcondense/GSM1541800_38-Per_1_K27ac_macs2_peaks.bed", bin_size=10000)
bin_k9ac = read_chip_file("/home/spark159/../../media/spark159/sw/dataforcondense/GSM1541804_38-Per_1_K9ac_macs2_peaks.bed", bin_size=10000)
bin_k4me3 = read_chip_file("/home/spark159/../../media/spark159/sw/dataforcondense/GSM1541808_38-Per_1_K4me3_macs2_peaks.bed", bin_size=10000)
bin_k36me3 = read_chip_file("/home/spark159/../../media/spark159/sw/dataforcondense/GSM1541812_38-Per_1_K36me3_solid_rseg.bed", bin_size=10000)
bin_k27me3 = read_chip_file("/home/spark159/../../media/spark159/sw/dataforcondense/GSM1541816_38-Per_1_K27me3_rseg.bed", bin_size=10000)
bin_k9me2 = read_chip_file("/home/spark159/../../media/spark159/sw/dataforcondense/GSM1541820_38-Per_1_K9me2_rseg.bed", bin_size=10000)
bin_k9me3 = read_chip_file("/home/spark159/../../media/spark159/sw/dataforcondense/GSM1541824_38-Per_1_K9me3_rseg.bed", bin_size=10000)
data_list = [ID_me, bin_k27ac, bin_k9ac, bin_k4me3, bin_k36me3, bin_k27me3, bin_k9me2, bin_k9me3]
#data_list = [ID_me]

binID_GC = {}
for chr in chr_binID_GC:
    for binID in range(len(chr_binID_GC[chr])):
        GC = chr_binID_GC[chr][binID] * 100
        ID = chr + ":" + str(binID)
        if ID not in binID_GC:
            binID_GC[ID] = 0.0
        binID_GC[ID] += GC

binID_rcount_list = []
for i in range(len(chr_binID_counts)-1):
    binID_rcount = {}
    chr_binID_count = chr_binID_counts[i]
    for chr in chr_binID_count:
        for binID in range(len(chr_binID_count[chr])):
            control = chr_binID_control[chr][binID]
            #if control <= 0:
                #continue
                #rcount = 1.0
            #else:
                #test = chr_binID_count[chr][binID]
                #rcount = float(test) / control
            test = chr_binID_count[chr][binID]
            rcount = float(test)
            ID = chr + ":" + str(binID)
            if ID not in binID_rcount:
                binID_rcount[ID] = 0
            binID_rcount[ID] += rcount
    binID_rcount_list.append(binID_rcount)

frac=[1 for i in range(1,11)]
data_name = ["CpGme", "k27ac", "k9ac", "k4me3", "k36me3", "k27me3", "k9me2", "k9me3"]
#data_name = ["CpGme"]
for i in range(len(data_list)):
    data = data_list[i]
    for j in range(len(binID_rcount_list)):
        binID_rcount = binID_rcount_list[j]
        #graphics.PartitionScatterplot (binID_GC, binID_rcount, data, frac, xlim=[0, 100], ylim=[-0,3], xlabel="GC content (%)", ylabel = 'Normalized Counts', note=fname + "_" + data_name[i] + "_" + str(j+1), title=data_name[i] + " Titration " + str(j+1))
        graphics.PartitionBoxplot (binID_rcount, data, frac, xlabel="", ylabel='Normalized Counts', title = data_name[i] + " Titration " + str(j+1), note="_" + data_name[i] + "_" + str(j+1))



"""
IDsig_list = [[] for i in range(len(chip_list))]
for ID in ID_rcount1:
    for i in range(len(chip_list)):
        ID_chip = chip_list[i]
        if ID not in ID_chip:
            sig = 0.0
        else:
            sig = ID_chip[ID]
        IDsig_list[i].append([ID, sig])

GC_rcount1, GC_rcount2 = {}, {}
for ID in ID_GC:
    GC = ID_GC[ID]*100
    GC = 100-GC
    if ID in ID_rcount1:
        rcount1 = ID_rcount1[ID]
        if GC not in GC_rcount1:
            GC_rcount1[GC] = []
        GC_rcount1[GC].append(rcount1)
    if ID in ID_rcount2:
        rcount2 = ID_rcount2[ID]
        if GC not in GC_rcount2:
            GC_rcount2[GC] = []
        GC_rcount2[GC].append(rcount2)

GC_mean1 = {}
for GC in GC_rcount1:
    GC_mean1[GC] = np.mean(GC_rcount1[GC])
GC_mean2 = {}
for GC in GC_rcount2:
    GC_mean2[GC] = np.mean(GC_rcount2[GC])
GCmean = [GC_mean1, GC_mean2]

def norm(L):
    total = sum(L)
    return [L[i]/float(total) for i in range(len(L))]

def quantile (bin_sig, num, frac=None):
    def value_cmp(a, b):
        if a[1] <= b[1]:
            return -1
        else:
            return 1
    bin_sig = sorted(bin_sig, cmp=value_cmp)
    if frac == None:
        size_list = [int(math.ceil(len(bin_sig) / float(num)))]*num
    else:
        if sum(frac) != 1:
            frac = norm(frac)
        num = len(frac)
        size_list = [int(round(f*len(bin_sig))) for f in frac]
    size_list[-1] += len(bin_sig) - sum(size_list)
    if size_list[-1] == 0:
        size_list[-2] -= 1
        size_list[-1] += 1
    output = []
    ed = 0
    for i in range(num):
        #st = i*size
        #ed = min((i+1)*size,len(bin_sig))
        size = size_list[i]
        st = ed             
        ed = st + size
        output.append(bin_sig[st:ed])
    return output

frac=[(4**i) for i in range(1,11)]
frac = norm(frac)[::-1]
print frac
print sum(frac)
num = len(frac)
#binBS_list = quantile(IDme, num=num, frac=frac)
alpha_list = np.linspace(0.02, 1, num=num)
color_list = np.linspace(0.01, 1, num=num)
cmap = mpl.cm.get_cmap("jet")
cmarks = ['b.', 'gd', 'y*', 'm^', 'rx']
#unit = "AU"
unit = 'AU'

chip_name = ["k27ac", "k9ac", "k4me3", "k36me3", "k27me3", "k9me2", "k9me3"]
#chip_name = ["k27ac", "k9ac", "k4me3"]
names = [chr_num + "_sp9_" + str(win_size), chr_num + "_sp10_" + str(win_size)]
for u in range(len(chip_list)):
    IDsig = IDsig_list[u]
    binBS_list = quantile(IDsig, num=num, frac=frac)
    for k in range(2):
        ID_rcount = ID_rcount_list[k]
        name = names[k]
        boxdata = []
        for i in range(len(binBS_list)):
            binBS = binBS_list[i]
            print len(binBS)
            if len(binBS) <= 0:
                continue
            minBS, maxBS = round(binBS[0][1],1), round(binBS[-1][1],1)
            if i == 0:
                label = "<" + str(maxBS) + unit
            elif i == len(binBS_list)-1:
                label = ">" + str(minBS) + unit
            else:
                label = str(minBS) + unit + "~" + str(maxBS) + unit
            X = [100*(1-ID_GC[bin]) for bin, BS in binBS]
            Y = [ID_rcount[bin] for bin, BS in binBS]
            plt.figure(1)
            plt.plot(X, Y, '.', color=cmap(color_list[i]), alpha=alpha_list[i], label=label)

            GC_rcount = {}
            x,y,z = [], [], []
            values = []
            for BinID, BS in binBS:
                GC = ID_GC[BinID]
                rcount = ID_rcount[BinID]
                if GC not in GC_rcount:
                    GC_rcount[GC] = []
                GC_rcount[GC].append(rcount)
                values.append(rcount - GCmean[k][(100-GC*100)])
            boxdata.append(values)
            for GC in GC_rcount:
                if len(GC_rcount[GC]) <=1:
                    continue
                x.append((1-GC)*100)
                y.append(np.nanmean(GC_rcount[GC]))
                z.append(np.nanstd(GC_rcount[GC])/np.sqrt(len(GC_rcount[GC])))
            plt.figure(2)
            plt.errorbar(x,y, yerr=z, fmt='.', color=cmap(color_list[i]))
            plt.plot(x,y,'.', color=cmap(color_list[i]), label=label)

        plt.figure(1)
        plt.xlabel('AT content (%)')
        plt.ylabel('Condensability (A.U.)')
        plt.ylim([-2,2.5])
        #plt.ylim([0,4])
        plt.xlim([0, 100])
        plt.title(chr_name)
        #plt.ylim([0,5.5])

        #plt.ylim([0.6,2.2])
        #plt.ylim([0.6,3.5])
        leg = plt.legend(loc='best', numpoints=1, prop={'size': 10})
        for lh in leg.legendHandles:
            lh._legmarker.set_markersize(15)
            lh._legmarker.set_alpha(1)
        plt.savefig(name + "_" + chip_name[u] + ".png")
        #plt.show()
        plt.close()    

        plt.figure(2)
        plt.xlabel('AT content (%)')
        plt.ylabel('Condensability (A.U.)')
        #plt.ylim([0,4])

        #plt.ylim([0,5.5])

        #plt.ylim([0.6,2.2])
        plt.ylim([-0.6, 0.6])
        plt.xlim([0, 100])
        plt.title(chr_name)
        #plt.ylim([0.6,3.5])
        leg = plt.legend(loc='best', numpoints=1, prop={'size': 10})
        for lh in leg.legendHandles:
            lh._legmarker.set_markersize(15)
            lh._legmarker.set_alpha(1)
        plt.savefig(name + "_" + chip_name[u] + "_mean.png")
        #plt.show()
        plt.close()

        plt.figure(3)
        plt.xlabel('Partitions by ' + chip_name[u])
        plt.ylabel('Condensability - GC dependence (A.U.)')
        plt.boxplot(boxdata, 0, "")
        plt.savefig(name + "_"+ chip_name[u] + "_box.png")
        #fig3.show()
        plt.close()
        #plt.ylim([0, 3])

"""
