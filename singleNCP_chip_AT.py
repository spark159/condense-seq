import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import copy
import math
from pyliftover import LiftOver
from Bio import SeqIO
import Interval_dict

def read_NCPcov_file (fname):
    ID_dyad = {}
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            names = cols[3:]
            ID_rcount_list = [{} for i in range(len(names)-1)]
            #bin_GC = {}
            First = False
            continue
        ID, chr, dyad = int(cols[0]), cols[1], int(cols[2])
        ID_dyad[ID] = dyad
        counts = cols[3:]
        if float(counts[-1]) <= 0:
            counts[-1] = 1.0
        for i in range(len(ID_rcount_list)):
            ID_rcount = ID_rcount_list[i]
            #bin_rcount[BinID] = np.log(float(counts[i])/float(counts[-1]))
            if float(counts[i]) <= 0:
                rcount = (float(counts[-1])/float(1.0))
            else:
                rcount = (float(counts[-1])/float(counts[i]))
            ID_rcount[ID] = np.log(rcount)
        #bin_GC[BinID] = float(GC)
    return ID_rcount_list, ID_dyad, names
#ID_rcount_list, ID_dyad, names = read_NCPcov_file("hg19_chr13_171_Ncov.cn")
#ID_rcount_list, ID_dyad, names = read_NCPcov_file("simple_Ncov.cn")
ID_rcount_list, ID_dyad, names = read_NCPcov_file("iNPS_Ncov.cn")
ID_rcount1, ID_rcount2 = ID_rcount_list
win_size = 171
chr_num = 'chr1'
chr_name = "chromosome 1"

st_ID = {}
ID_interval = {}
for ID in ID_dyad:
    dyad = ID_dyad[ID]
    st = dyad-win_size/2
    ed = dyad+win_size/2+1
    st_ID[st] = ID
    ID_interval[ID] = [st, ed]

def get_GC(ref_fname, chr, st_ID, win):
    def GC_content(seq):
        seq = seq.upper()
        output = 0.0
        for nt in seq:
            if nt in "GC":
                output += 1.0
        return output/len(seq)
    seq = ""
    pt = -1
    k = 0
    left = []
    Find = False
    stID = [[st, st_ID[st]] for st in sorted(st_ID.keys())]
    pos, ID = stID[k]
    ID_GC = {}
    for line in open(ref_fname):
        line = line.strip()
        if line.startswith(">"):
            if Find:
                break
            if line[1:] == chr:
                Find = True
            continue
        if Find:
            if len(left) == 0 and pt + len(line) < pos:
                pt += len(line)
                continue
            for i in range(len(left)):
                leftID, seq = left.pop(0)
                ed = min(len(line), win-len(seq))
                seq += line[:ed]
                if len(seq) == win:
                    GC = GC_content(seq)
                    ID_GC[leftID] = GC
                else:
                    left.append([leftID, seq])
            while pt + len(line) >= pos and k < len(stID):
                loc = pos - pt - 1
                seq = line[loc:min(loc+win,len(line))]
                if len(seq) == win:
                    GC = GC_content(seq)
                    ID_GC[ID] = GC
                else:
                    left.append([ID, seq])
                k += 1
                try:
                    pos, ID = stID[k]
                except:
                    None
            if len(left) == 0 and len(ID_GC) == len(stID):
                break
            pt += len(line)
    return ID_GC
ID_GC = get_GC("hg19.fa", chr_num, st_ID, win=win_size)

Int_dict1 = Interval_dict.double_hash(ID_interval, 10000, 1000000000)
Int_dict2 = copy.deepcopy(Int_dict1)
Int_dict3 = copy.deepcopy(Int_dict1)
Int_dict4 = copy.deepcopy(Int_dict1)
Int_dict5 = copy.deepcopy(Int_dict1)
Int_dict6 = copy.deepcopy(Int_dict1)
Int_dict7 = copy.deepcopy(Int_dict1)

def read_chip_file (fname, Int_dict, chr):
    count = -1
    for line in open(fname):
        count += 1
        print count
        cols = line.strip().split()
        chrname, st, ed = cols[0], int(cols[1]), int(cols[2])
        if chrname != chr:
            continue
        score = float(cols[4])
        Int_dict.insert_range(st, ed+1, score)
    print "Chip reading is done"
    return Int_dict.get()

"""
def read_chip_file (fname, st_ID, chr, win):
    ID_sig = {}
    stID = [[st, st_ID[st]] for st in sorted(st_ID.keys())]
    print len(stID)
    pos, ID = stID.pop(0)
    left = []
    count = -1
    for line in open(fname):
        count += 1
        print count
        cols = line.strip().split()
        chrname, st, ed = cols[0], int(cols[1]), int(cols[2])
        if chrname != chr:
            continue
        score = float(cols[4])
        for i in range(len(left)):
            leftpos, leftID = left.pop(0)
            last = min(ed, leftpos+win-1)
            if st > last:
                continue
            length = last - st + 1
            ID_sig[leftID] += score*length
            if (last - leftpos + 1) < win:
                left.append([leftpos, leftID])
        while pos < st:
            last = min(ed, pos+win-1)
            if st <= last:
                length = last - st + 1
                if ID not in ID_sig:
                    ID_sig[ID] = 0.0
                ID_sig[ID] += score*length
                if (last - pos + 1) < win:
                    left.append([pos, ID])
            pos, ID = stID.pop(0)
        while pos >=st and pos <= ed:
            last = min(ed, pos+win-1)
            length = last - pos + 1
            if ID not in ID_sig:
                ID_sig[ID] = 0.0
            ID_sig[ID] += score*length
            if length < win:
                left.append([pos, ID])
            pos, ID = stID.pop(0)
    return ID_sig
"""
bin_k27ac = read_chip_file("GSM1541800_38-Per_1_K27ac_macs2_peaks.bed", Int_dict1, chr=chr_num)
bin_k9ac = read_chip_file("GSM1541804_38-Per_1_K9ac_macs2_peaks.bed", Int_dict2, chr=chr_num)
bin_k4me3 = read_chip_file("GSM1541808_38-Per_1_K4me3_macs2_peaks.bed", Int_dict3, chr=chr_num)
bin_k36me3 = read_chip_file("GSM1541812_38-Per_1_K36me3_solid_rseg.bed", Int_dict4, chr=chr_num)
bin_k27me3 = read_chip_file("GSM1541816_38-Per_1_K27me3_rseg.bed", Int_dict5, chr=chr_num)
bin_k9me2 = read_chip_file("GSM1541820_38-Per_1_K9me2_rseg.bed", Int_dict6, chr=chr_num)
bin_k9me3 = read_chip_file("GSM1541824_38-Per_1_K9me3_rseg.bed", Int_dict7, chr=chr_num)

chip_list = [bin_k27ac, bin_k9ac, bin_k4me3, bin_k36me3, bin_k27me3, bin_k9me2, bin_k9me3]
#chip_list = [bin_k27ac, bin_k9ac, bin_k4me3]
#print bin_k27ac

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

