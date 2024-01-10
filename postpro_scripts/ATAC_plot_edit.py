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
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn
from matplotlib.colors import LinearSegmentedColormap


pastel_jet = LinearSegmentedColormap.from_list('white_viridis',
                                             [(0, '#ffffff'),
                                              (0.05, 'tab:cyan'),
                                              (0.1, 'tab:blue'),
                                              (0.3, 'tab:green'),
                                              (0.5, 'yellow'),
                                              (0.7, 'tab:orange'),
                                              (0.9, 'tab:red'),
                                              (1, 'darkred')
                                             ], N=256)

def density_scatter(x , y, ax = None, sort = True, bins = 20, density = False, **kwargs )   :
    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d(x, y, bins = bins, density=density )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)
    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0
    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
    img = ax.scatter( x, y, c=z, **kwargs )
    cbar = plt.colorbar(img)
    cbar.ax.tick_params(labelsize=5)
    #cbar = plt.colorbar(cm.ScalarMappable(norm = norm), ax=img)
    #cbar.ax.set_ylabel('Density')
    return ax

def read_ATAC (fname, chr_choice):
    lID_score = {}
    lID_interval = {}
    for line in open(fname):
        cols = line.strip().split()
        chr, st, ed, id, _, _, counts, score, _ = cols
        if chr != chr_choice:
            continue
        st, ed = int(st), int(ed)
        counts = int(counts)
        #score = float(score)
        score = np.log2(counts+1.0)
        lID_score[id] = score
        lID_interval[id] = (st, ed)
    return lID_score, lID_interval

def read_ATAC_foldchange (fname, chr_choice):
    lID_score = {}
    lID_interval = {}
    lID = 0
    for line in open(fname):
        cols = line.strip().split()
        chr, st, ed, foldchange = cols
        if chr != chr_choice:
            continue
        st, ed = int(st), int(ed)
        foldchange = float(foldchange)
        lID_score[lID] = np.log2(1.0+foldchange)
        lID_interval[lID] = (st, ed)
        lID +=1
    return lID_score, lID_interval

def read_binfile (fname, chr_choice):
    ID_score1, ID_score2 = {}, {}
    ID_interval = {}
    for line in open(fname):
        cols = line.strip().split()
        binID, chr, st, ed, count1, count2, count0, GC = cols
        if chr != chr_choice:
            continue
        binID = int(binID)
        st, ed = int(st), int(ed)
        count1 = 1 + float(count1)
        count2 = 1 + float(count2)
        count0 = 1 + float(count0)
        foldchange1 = count1/count0
        foldchange2 = count2/count0
        ID_score1[binID] = -np.log(foldchange1)
        ID_score2[binID] = -np.log(foldchange2)
        ID_interval[binID] = (st, ed)
    return ID_interval, ID_score1, ID_score2
    
path = ""
ID_interval, ID_score1, ID_score2 = read_binfile(path + "H1_NCP_sp_1kb_bin.cn", chr_choice="chr1")
lID_score, lID_interval = read_ATAC_foldchange(path + "H1_ATAC_foldchange.bedgraph", chr_choice="chr1")

interval_dict = Interval_dict.double_hash(ID_interval, 100000, 250000000)

for lID in lID_score:
    rst, red = lID_interval[lID]
    interval_dict.insert_range(rst, red, lID_score[lID])

ID_lscore = interval_dict.get()
del interval_dict

X, Y = [], []
for ID in ID_lscore:
    score2 = ID_score2[ID]
    lscore = ID_lscore[ID] /1000
    X.append(score2)
    Y.append(lscore)

X, Y = np.asarray(X), np.asarray(Y)

fig = plt.figure()
#plt.plot(X, Y, '.', markersize=1, alpha=0.1)
density_scatter(X, Y, bins = [20,20], s=3, cmap=pastel_jet, ax=plt.gca())
plt.xlabel("Condensability")
plt.ylabel("ATAC score")
plt.title("Condensability VS ATAC score (chr1, 1kb)")
plt.xlim([-2.5, 2.5])
plt.ylim([-0.1, 3])
plt.show()
plt.close()
