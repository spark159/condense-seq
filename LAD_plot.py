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

lo = LiftOver("hg18", "hg19")

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

lID_pos1, lID_score1 = read_data("data/GSM1612855-28369.txt", "chr1")
lID_pos2, lID_score2 = read_data("data/GSM1612856-28372.txt", 'chr1')

keys = list(set(lID_pos1.keys()) & set(lID_pos2.keys()))

lID_pos, lID_score = {}, {}
lID_interval = {}
for key in keys:
    assert lID_pos1[key] == lID_pos2[key]
    pos = lID_pos1[key]
    lID_pos[key] = pos
    lID_score[key] = np.mean([lID_score1[key], lID_score2[key]])
    lID_interval[key] = (pos-500, pos+500)

linterval_dict = Interval_dict.double_hash(lID_interval, 100000, 250000000)

ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("data/hg19_chr1_171_everything_anot.cn")

name_lID_values = {}
for ID in ID_pos:
    pos = ID_pos[ID]
    lIDs = linterval_dict.find(pos)
    if len(lIDs) <= 0:
        continue
    for name in name_ID_value:
        if name not in name_lID_values:
            name_lID_values[name] = {}
        value = name_ID_value[name][ID]
        for lID in lIDs:
            if lID not in name_lID_values[name]:
                name_lID_values[name][lID] = []
            name_lID_values[name][lID].append(value)
            
name_lID_mean = {}
for name in name_lID_values:
    if name not in name_lID_mean:
        name_lID_mean[name] = {}
    for lID in name_lID_values[name]:
        assert lID not in name_lID_mean[name]
        name_lID_mean[name][lID] = np.mean(name_lID_values[name][lID])

names = name_lID_mean.keys()
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
    plt.plot(X, Y, ',', alpha=0.2)
    plt.plot(X, Ypred, 'r--')
    plt.xlabel("DamID score")
    plt.ylabel(name)
    plt.title("1kb window")
    plt.savefig("LAD_scatter_" + name + ".png",bbox_inches='tight')
    #plt.show()
    plt.close()
