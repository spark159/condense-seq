import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import copy
import math
from pyliftover import LiftOver
from Bio import SeqIO
import Interval_dict

"""
def read_anot_file(fname):
    ID_dyad = {}
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            names = cols[3:]
            ID_values = [{} for i in range(len(names))]
            chip_names = cols[7:]
            First = False
            continue
        ID, chr, dyad = int(cols[0]), cols[1], int(cols[2])
        ID_dyad[ID] = dyad
        cols = cols[3:]
        cols = [ float(value) for value in cols]
        if sum(cols) == 0:
            continue
        for i in range(len(cols)):
            ID_value = ID_values[i]
            value = float(cols[i])
            if ID not in ID_value:
                ID_value[ID] = value
    ID_metric_list = ID_values[0:2]
    ID_AT = ID_values[2]
    ID_me = ID_values[3]
    ID_chip_list = ID_values[4:]
    return ID_dyad, ID_metric_list, ID_AT, ID_me, ID_chip_list, chip_names
"""
def read_anot_file(fname):
    ID_pos = {}
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            names = cols[3:]
            chip_names = cols[8:]
            ID_values = [{} for i in range(len(names)-1)]
            First = False
            continue
        ID, chr, pos = int(cols[0]), cols[1], int(cols[2])
        ID_pos[ID] = pos
        lastcols = cols[3:6] + cols[7:]
        for i in range(len(lastcols)):
            ID_value = ID_values[i]
            if i == 3:
                try:
                    value = float(lastcols[i])#/float(cols[6])
                except:
                    value = float(lastcols[i])
                    #value = 0
            else:
                value = float(lastcols[i])
            if ID not in ID_value:
                ID_value[ID] = value
    ID_metric_list = ID_values[0:2]
    ID_AT = ID_values[2]
    ID_me = ID_values[3]
    ID_chip_list = ID_values[4:]
    return ID_pos, ID_metric_list, ID_AT, ID_me, ID_chip_list, chip_names

ID_dyad, ID_rcount_list, ID_AT, ID_me, chip_list, chip_name = read_anot_file("data/hg19_chr1_171_new_anot.cn")
#ID_dyad, ID_rcount_list, ID_AT, ID_me, chip_list, chip_name = read_anot_file("data/hg19_chr1_1001win_anot.cn")
ID_rcount1, ID_rcount2 = ID_rcount_list
win_size = 171
chr_num = 'chr1'
chr_name = "chromosome 1"


X1, Y1 = [], []
X2, Y2 = [], []
X3, Y3 = [], []
AT_rcount1, AT_rcount2 = {}, {}
IDme = []

IDs = list(set(ID_rcount1.keys()) | set(ID_rcount2.keys()))

for ID in IDs:
    rcount1 = ID_rcount1[ID]
    rcount2 = ID_rcount2[ID]
    AT = ID_AT[ID]*100
    X1.append(AT)
    Y1.append(rcount1)
    X2.append(AT)
    Y2.append(rcount2)
    if AT not in AT_rcount1:
        AT_rcount1[AT] = []
    if AT not in AT_rcount2:
        AT_rcount2[AT] = []
    AT_rcount1[AT].append(rcount1)
    AT_rcount2[AT].append(rcount2)
    AT = (AT/100) * win_size
    if ID not in ID_me:
        me = 0.0
    else:
        me = ID_me[ID]/AT
    me = me*100
    IDme.append([ID,me])

x1, y1, z1 = [], [], []
x2, y2, z2 = [], [], []
for AT in AT_rcount1:
    x1.append(AT)
    y1.append(np.nanmean(AT_rcount1[AT]))
    z1.append(np.nanstd(AT_rcount1[AT])/np.sqrt(len(AT_rcount1[AT])))
for AT in AT_rcount2:
    x2.append(AT)
    y2.append(np.nanmean(AT_rcount2[AT]))
    z2.append(np.nanstd(AT_rcount2[AT])/np.sqrt(len(AT_rcount2[AT])))

fig = plt.figure()
plt.plot(X1,Y1,'k.', alpha=0.01, label='sp9')
plt.errorbar(x1,y1, yerr=z1, fmt='.')
plt.plot(x1,y1,'.')
plt.xlabel('AT content (%)')
plt.ylabel('Condensability (A.U.)')
plt.title(chr_name)
#plt.ylim([0,4])
plt.ylim([-2,2.5])
plt.xlim([0, 100])
leg = plt.legend(loc='best', numpoints=1, prop={'size': 6})
for lh in leg.legendHandles:
    lh._legmarker.set_markersize(15)
    lh._legmarker.set_alpha(1)
plt.savefig(chr_num+'_sp9_' + str(win_size) + '.png')
plt.close()

fig = plt.figure()
plt.plot(X2,Y2,'k.', alpha=0.01, label='sp10')
plt.errorbar(x2,y2, yerr=z2, fmt='.')
plt.plot(x2,y2,'.')
plt.xlabel('AT content (%)')
plt.ylabel('Condensability (A.U.)')
plt.title(chr_name)
#plt.ylim([0,4])
plt.ylim([-2,2.5])
plt.xlim([0, 100])
leg = plt.legend(loc='best', numpoints=1, prop={'size': 6})
for lh in leg.legendHandles:
    lh._legmarker.set_markersize(15)
    lh._legmarker.set_alpha(1)
plt.savefig(chr_num+'_sp10_' + str(win_size) + '.png')
plt.close()

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
binBS_list = quantile(IDme, num=num, frac=frac)
alpha_list = np.linspace(0.02, 1, num=num)
color_list = np.linspace(0.01, 1, num=num)
cmap = mpl.cm.get_cmap("jet")
cmarks = ['b.', 'gd', 'y*', 'm^', 'rx']
#unit = "AU"
unit = '%'


names = [chr_num + "_sp9_" + str(win_size), chr_num + "_sp10_" + str(win_size)]
for k in range(2):
    ID_rcount = ID_rcount_list[k]
    name = names[k]
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
        X = [100*(ID_AT[bin]) for bin, BS in binBS]
        Y = [ID_rcount[bin] for bin, BS in binBS]
        plt.figure(1)
        plt.plot(X, Y, '.', color=cmap(color_list[i]), alpha=alpha_list[i], label=label)
    
        AT_rcount = {}
        x,y,z = [], [], []
        for BinID, BS in binBS:
            AT = ID_AT[BinID]
            rcount = ID_rcount[BinID]
            if AT not in AT_rcount:
                AT_rcount[AT] = []
            AT_rcount[AT].append(rcount)
        for AT in AT_rcount:
            if len(AT_rcount[AT]) <=1:
                continue
            x.append(AT*100)
            y.append(np.nanmean(AT_rcount[AT]))
            z.append(np.nanstd(AT_rcount[AT])/np.sqrt(len(AT_rcount[AT])))
        plt.figure(2)
        plt.errorbar(x,y, yerr=z, fmt='.', color=cmap(color_list[i]))
        plt.plot(x,y,'.', color=cmap(color_list[i]), label=label)

        
    plt.figure(1)
    plt.xlabel('AT content (%)')
    plt.ylabel('Condensability (A.U.)')
    #plt.ylim([0,4])
    plt.ylim([-2,2.5])
    plt.xlim([0, 100])
    plt.title(chr_name)
    #plt.ylim([0,5.5])

    #plt.ylim([0.6,2.2])
    #plt.ylim([0.6,3.5])
    leg = plt.legend(loc='best', numpoints=1, prop={'size': 6})
    for lh in leg.legendHandles:
        lh._legmarker.set_markersize(15)
        lh._legmarker.set_alpha(1)
    plt.savefig(name + "_me.png")
    #plt.show()
    plt.close()    

    plt.figure(2)
    plt.xlabel('AT content (%)')
    plt.ylabel('Condensability (A.U.)')
    #plt.ylim([0,4])

    #plt.ylim([0,5.5])

    #plt.ylim([0.6,2.2])
    #plt.ylim([-2,2])
    plt.ylim([-0.8, 0.6])
    plt.xlim([0, 100])
    plt.title(chr_name)
    #plt.ylim([0.6,3.5])
    leg = plt.legend(loc='best', numpoints=1, prop={'size': 10})
    for lh in leg.legendHandles:
        lh._legmarker.set_markersize(15)
        lh._legmarker.set_alpha(1)
    plt.savefig(name + "_me_mean.png")
    #plt.show()
    plt.close()


IDsig_list = [[] for i in range(len(chip_list))]
for ID in ID_rcount1:
    for i in range(len(chip_list)):
        ID_chip = chip_list[i]
        if ID not in ID_chip:
            sig = 0.0
        else:
            sig = ID_chip[ID]
        IDsig_list[i].append([ID, sig])


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

#chip_name = ["k27ac", "k9ac", "k4me3"]
names = [chr_num + "_sp9_" + str(win_size), chr_num + "_sp10_" + str(win_size)]
for u in range(len(chip_list)):
    IDsig = IDsig_list[u]
    binBS_list = quantile(IDsig, num=num, frac=frac)
    for k in range(2):
        ID_rcount = ID_rcount_list[k]
        name = names[k]
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
            X = [100*ID_AT[bin] for bin, BS in binBS]
            Y = [ID_rcount[bin] for bin, BS in binBS]
            plt.figure(1)
            plt.plot(X, Y, '.', color=cmap(color_list[i]), alpha=alpha_list[i], label=label)

            AT_rcount = {}
            x,y,z = [], [], []
            for BinID, BS in binBS:
                AT = ID_AT[BinID]
                rcount = ID_rcount[BinID]
                if AT not in AT_rcount:
                    AT_rcount[AT] = []
                AT_rcount[AT].append(rcount)
            for AT in AT_rcount:
                if len(AT_rcount[AT]) <=1:
                    continue
                x.append(AT*100)
                y.append(np.mean(AT_rcount[AT]))
                z.append(np.std(AT_rcount[AT])/np.sqrt(len(AT_rcount[AT])))
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
        leg = plt.legend(loc='best', numpoints=1, prop={'size': 6})
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
        leg = plt.legend(loc='best', numpoints=1, prop={'size': 6})
        for lh in leg.legendHandles:
            lh._legmarker.set_markersize(15)
            lh._legmarker.set_alpha(1)
        plt.savefig(name + "_" + chip_name[u] + "_mean.png")
        #plt.show()
        plt.close()
