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
from pyliftover import LiftOver

def get_corr(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = np.average(x)
    avg_y = np.average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff
    return diffprod / np.sqrt(xdiff2 * ydiff2)

#lo = LiftOver("hg18", "hg19")

def read_data (fname, chr_target):
    ID_pos, ID_score = {}, {}
    count = -1
    for line in open(fname):
        count +=1
        if count < 3:
            continue
        ID, score = line.strip().split()
        chr_info, pos_info = ID.split('P')
        try:
            chr = 'chr' + str(int(chr_info[3:]))
        except:
            chr = 'chr' + chr_info[3:]
        pos = int(pos_info[:-2])
        new_pos_list = lo.convert_coordinate(chr, pos)
        if len(new_pos_list) != 1:
            continue
        chr = new_pos_list[0][0]
        if chr != chr_target:
            continue
        pos = new_pos_list[0][1]
        strand = new_pos_list[0][2]
        score = float(score)
        assert ID not in ID_pos
        assert ID not in ID_score
        ID_pos[ID] = pos
        ID_score[ID] = score
    return ID_pos, ID_score

lo = LiftOver("hg38", "hg19")
def read_DamID (fname, chr_choice='chr1'):
    i = 0
    lID_value = {}
    lID_interval = {}
    for line in open(fname):
        cols = line.strip().split()
        chr, st, ed, value = cols
        if chr != chr_choice:
            continue
        st, ed = int(st), int(ed)
        new_st_list = lo.convert_coordinate(chr, st)
        new_ed_list = lo.convert_coordinate(chr, ed)
        #print new_st_list
        #print new_ed_list
        if len(new_st_list) <=0 or len(new_ed_list) <=0:
            continue
        if chr != new_st_list[0][0] or chr != new_ed_list[0][0]:
            continue
        st = new_st_list[0][1]
        ed = new_ed_list[0][1]
        lID_interval[i] = (st, ed)
        lID_value[i] = float(value)
        i +=1
    return lID_value, lID_interval

path = "./data/"

#lID_pos1, lID_score1 = read_data(path + "GSM1612855-28369.txt", "chr1")
#lID_pos2, lID_score2 = read_data(path + "GSM1612856-28372.txt", 'chr1')

#keys = list(set(lID_pos1.keys()) & set(lID_pos2.keys()))

#lID_pos, lID_score = {}, {}
#lID_interval = {}
#for key in keys:
#    assert lID_pos1[key] == lID_pos2[key]
#    pos = lID_pos1[key]
#    lID_pos[key] = pos
#    lID_score[key] = np.mean([lID_score1[key], lID_score2[key]])
#    lID_interval[key] = (pos-500, pos+500)

    
lID_score, lID_interval = read_DamID(path + "HCT116_LMNB1_DamID.bedgraph", chr_choice="chr1")
linterval_dict = Interval_dict.double_hash(lID_interval, 100000, 250000000)

#ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path + "hg19_chr1_NCP_anot.cn")
name_ID_value = load_file.read_tabular_file (path + "hg19_chr1_167win25step_anot.cn", mode='col', jump=10)
ID_pos = name_ID_value['PhysicalPosition']

#names = ['data/sp_spd_tests_detail/sp7', 'ATcontent', 'CpGNumber', 'k27ac', 'k9ac', 'k4me3', 'k36me3', 'k9me2', 'k9me3', 'k27me3']
#names = ['data/sp_spd_tests_detail/sp7', 'ATcontent', 'CpGNumber', 'k27ac', 'k9ac', 'k4me3', 'k36me3', 'k9me2', 'k9me3', 'k27me3']

ID_CG = name_ID_value['CpGNumber']
ID_me = name_ID_value['meGCNumber']
ID_mefrac = {}
for ID in ID_CG:
    CG = ID_CG[ID]
    if CG <= 0:
        continue
    me = ID_me[ID]
    mefrac = float(me) / (2*CG)
    ID_mefrac[ID] = mefrac
name_ID_value['meCpG density'] = ID_mefrac
#names.append('meCpG density')

#names = name_ID_value.keys()
names = ['meGCNumber', 'k27me3', 'k27ac', 'k9me3', 'k9me2', 'ATcontent', 'k36me3', 'meCpG density', 'k4me3', 'k9ac', 'CpGNumber', 'data/sp_spd_tests_detail/sp7']

name_lID_values = {}
for ID in ID_pos:
    pos = ID_pos[ID]
    lIDs = linterval_dict.find(pos)
    if len(lIDs) <= 0:
        continue
    for name in names:
        if name not in name_lID_values:
            name_lID_values[name] = {}
        try:
            value = name_ID_value[name][ID]
        except:
            continue
        if np.isnan(value):
            continue
        for lID in lIDs:
            if lID not in name_lID_values[name]:
                name_lID_values[name][lID] = []
            name_lID_values[name][lID].append(value)
            
name_lID_mean = {}
for name in names:
    if name not in name_lID_mean:
        name_lID_mean[name] = {}
    for lID in name_lID_values[name]:
        assert lID not in name_lID_mean[name]
        name_lID_mean[name][lID] = np.mean(name_lID_values[name][lID])

X_list, Y_list = [], []
for i in range(len(names)):
    name = names[i]
    lID_mean = name_lID_mean[name]
    X, Y = [], []
    for lID in lID_mean:
        score = lID_score[lID]
        mean = lID_mean[lID]
        X.append(score)
        Y.append(mean)
    X_list.append(X)
    Y_list.append(Y)

corr_list = []
for i in range(len(names)):
    name = names[i].split('/')[-1]
    X, Y = X_list[i], Y_list[i]
    corr = get_corr(X, Y)
    corr_list.append(corr)
    print name, corr
    feature_list = [[x] for x in X]
    test_list = [[y] for y in Y]
    reg = linear_model.Ridge(alpha=0.5)
    reg.fit (feature_list, test_list)
    Ypred = reg.predict(feature_list)
    Ypred = [ value[0] for value in Ypred]
    fig = plt.figure()
    plt.plot(X, Y, ',', alpha=0.5)
    plt.plot(X, Ypred, 'r--')
    xloc, yloc = np.mean([np.median(X), max(X)]), np.mean([np.median(Y), max(Y)])
    if name == 'sp7':
        yloc -= 1
    plt.text(xloc, yloc, str(round(corr,3)), fontsize=20, va='center', ha='center')
    plt.xlabel("DamID score")
    if name == 'sp7':
        plt.ylim([-2,2])
    if name == "sp7":
        name = "Condensability (A.U.)"
    plt.ylabel(name)
    plt.title("5kb window")
    #plt.savefig("LAD_scatter_" + name + ".png",bbox_inches='tight')
    #plt.show()
    plt.close()

names[-1] = "Condensability"
fig = plt.figure()
plt.bar(range(len(corr_list)), corr_list, width=0.5, color='g')
plt.xticks(range(len(corr_list)), names, rotation=90)
plt.ylabel("Pearson Correlation")
plt.axhline(y=0, color='k', linestyle='--')
plt.title("Correlation with DamID score")
plt.savefig("bar_corr.png",bbox_inches='tight')
#plt.show()
plt.close()
