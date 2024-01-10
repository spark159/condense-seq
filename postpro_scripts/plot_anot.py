import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import seaborn as sns
import copy
import math
import graphics
import sys

win_size = 171
chr_target = 'chr1'
chr_name = "chromosome 1"

def norm(L):
    total = sum(L)
    return [L[i]/float(total) for i in range(len(L))]

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

#ID_pos, ID_score_list, ID_AT, ID_me, ID_chip_list, chip_names = read_anot_file("hg19_" + chr_target + "_" + str(win_size) + "_anot.cn")
ID_pos, ID_score_list, ID_AT, ID_me, ID_chip_list, chip_names = read_anot_file("data/hg19_chr1_171_new_anot.cn")
#ID_pos, ID_score_list, ID_AT, ID_me, ID_chip_list, chip_names = read_anot_file("data/hg19_chr1_1001win_anot.cn")
ID_score1, ID_score2 = ID_score_list

#X, Y = [], []
#Z = []
#for ID in ID_AT:
#    X.append(ID_AT[ID])
#    Y.append(ID_me[ID])
#    Z.append(ID_me[ID]/(ID_AT[ID]))
#fig = plt.figure()
#plt.plot(X, Y, '.', alpha=0.2)
#plt.plot(X, Z, '.', alpha=0.2)
#plt.show()

#graphics.draw_along_genome(ID_pos, ID_score2, ID_AT, ID_me)

#sys.exit()

def quantile (ID_score, num, frac=None):
    def value_cmp(a, b):
        if a[1] <= b[1]:
            return -1
        else:
            return 1
    IDscore = [[ID, score] for ID, score in ID_score.items()]
    IDscore = sorted(IDscore, cmp=value_cmp)
    if frac == None:
        size_list = [int(math.ceil(len(IDscore) / float(num)))]*num
    else:
        if sum(frac) != 1:
            frac = norm(frac)
        num = len(frac)
        size_list = [int(round(f*len(IDscore))) for f in frac]
    size_list[-1] += len(IDscore) - sum(size_list)
    if size_list[-1] == 0:
        size_list[-2] -= 1
        size_list[-1] += 1
    assert sum(size_list) == len(IDscore)
    output = []
    ed = 0
    for i in range(num):
        #st = i*size
        #ed = min((i+1)*size,len(IDscore))
        size = size_list[i]
        st = ed             
        ed = st + size
        temp = [IDscore[j][0] for j in range(st,ed)]
        output.append(temp)
    return output
frac=[(4**i) for i in range(1,11)]
#frac=[(3**i) for i in range(1,11)]
#frac=[(2**i) for i in range(1,11)]
#frac=[(5*i) for i in range(1,11)]
#frac=[ 1.0 for i in range(100)]
frac = norm(frac)[::-1]
group1 = quantile(ID_score1, len(frac), frac=frac)
group2 = quantile(ID_score2, len(frac), frac=frac)
group = quantile(ID_me, len(frac), frac=frac)


def neutralize_score_by_target (ID_score, ID_target):
    def standardization (data_list):
        if len(data_list) <= 1:
            return data_list
        mean = np.mean(data_list)
        std = np.std(data_list)
        return [ float(data - mean)/std for data in data_list]
    target_scores = {}
    target_IDs = {}
    for ID in ID_score:
        target = ID_target[ID]
        score = ID_score[ID]
        if target not in target_scores:
            target_scores[target] = []
            target_IDs[target] = []
        target_scores[target].append(score)
        target_IDs[target].append(ID)
    new_ID_score = {}
    for target in target_scores:
        new_scores = standardization(target_scores[target])
        IDs = target_IDs[target]
        for i in range(len(IDs)):
            ID = IDs[i]
            new_score = new_scores[i]
            new_ID_score[ID] = new_score
    return new_ID_score

new_ID_score2 = neutralize_score_by_target(ID_score2, ID_AT)

def PartitionBoxplot (ID_score, ID_target, frac, xlabel=""):
    ylabel = 'Condensability - GC dependence (A.U.)'
    frac = sorted(norm(frac), reverse=True)
    group = quantile(ID_target, len(frac), frac=frac)
    boxdata = []
    for IDs in group:
        print len(IDs)
        temp = []
        for ID in IDs:
            temp.append(ID_score[ID])
        boxdata.append(temp)
    fig = plt.figure()
    plt.xlabel('Partitions by ' + xlabel)
    plt.ylabel(ylabel)
    plt.boxplot(boxdata, 0, "")
    plt.savefig(xlabel + '_vs_' + ylabel + "_box.png")
    #plt.show()
    plt.close()

PartitionBoxplot (new_ID_score2, ID_me, frac, xlabel='meGC')
for i in range(len(ID_chip_list)):
    ID_chip = ID_chip_list[i]
    PartitionBoxplot(new_ID_score2, ID_chip, frac, xlabel=chip_names[i])



"""
#frac=[ 1.0 for i in range(10)]
#frac = norm(frac)[::-1]
#group = quantile(ID_AT, len(frac), frac=frac)
fig = plt.figure()
X = []; Y = []
Yerror = []
for i in range(len(group)):
    X.append(i)
    temp = []
    for ID in group[i]:
        temp.append(new_ID_score2[ID])
    Y.append(np.mean(temp))
    Yerror.append(np.std(temp)/np.sqrt(len(temp)))
plt.plot(X, Y, 'o')
plt.errorbar(X,Y, yerr=Yerror, fmt='o')
plt.show()


X, Y, Z = [], [], []
for ID in ID_score2:
    X.append(ID_AT[ID])
    Y.append(ID_score2[ID])
    Z.append(ID_me[ID])

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(X, Y, Z, c = 'b', marker='o')
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
plt.show()

X, Y = [], []
Xerror, Yerror = [], []
for IDs in group:
    tempX, tempY = [], []
    for ID in IDs:
        tempX.append(ID_AT[ID]*100)
        tempY.append(ID_score2[ID])
        #tempX.append(ID_score2[ID]*100)
        #tempY.append(ID_chip_list[2][ID])
        #if ID_AT[ID] >= 1:
        #    assert ID_me[ID] <= 0
        #    tempY.append(0.0)
        #else:
        #    tempY.append((float(ID_me[ID])/(win_size*(1-ID_AT[ID])))*100)
    X.append(np.mean(tempX))
    Xerror.append(np.std(tempX)/np.sqrt(len(tempX)))
    #Xerror.append(np.std(tempX))
    Y.append(np.mean(tempY))
    Yerror.append(np.std(tempY)/np.sqrt(len(tempY)))
    #Yerror.append(np.std(tempY))

fig = plt.figure()
cmap = mpl.cm.get_cmap('jet')
for i in range(len(frac)):
    rank = sum(frac[0:i+1])
    plt.plot(X[i],Y[i],'o', color=cmap(rank),  label="Partition" + str(i+1))
    plt.errorbar(X[i],Y[i], xerr=Xerror[i], color=cmap(rank), fmt='o')
    plt.errorbar(X[i],Y[i], yerr=Yerror[i], color=cmap(rank), fmt='o')
plt.xlabel('AT content (%)')
plt.ylabel('CpG methylation')
plt.title(chr_target)
#plt.legend()
#plt.ylim([0,4])
#plt.ylim([-2,2.5])
#plt.xlim([0, 100])
#leg = plt.legend(loc='best', numpoints=1, prop={'size': 10})
#for lh in leg.legendHandles:
#    lh._legmarker.set_markersize(15)
#    lh._legmarker.set_alpha(1)
#plt.savefig(chr_num+'_sp9_' + str(win_size) + '.png')
#plt.close()
plt.show()

"""
