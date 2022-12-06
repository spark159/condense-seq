import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
import Interval_dict
from scipy import stats
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import interpn
import subprocess

def Z_score (ID_value):
    mean = np.mean(ID_value.values())
    std = np.std(ID_value.values())
    ID_zscore = {}
    for ID in ID_value:
        zscore = float(ID_value[ID] - mean)/std
        ID_zscore[ID] = zscore
    return ID_zscore

def write_cn (ID_value, fname, IDs=None):
    f = open(fname, 'w')
    print >> f, 'SNP\tChromosome\tPhysicalPosition\tValue'
    if IDs == None:
        IDs = ID_value.keys()
    count = 0
    for ID in sorted(IDs):
        chrnum, st, ed = ID
        chr = 'chr' + str(chrnum)
        st, ed = int(st), int(ed)
        pos = int(round(0.5*(st+ed)))
        s = [str(count), chr, str(pos), str(ID_value[ID])]
        print >> f, '\t'.join(s)
        count +=1
    f.close()
        

# "jet-like" colormap with white background
pastel_jet = LinearSegmentedColormap.from_list('white_viridis',
                                             [(0, '#ffffff'),
                                              (0.03, 'tab:cyan'),
                                              (0.1, 'tab:blue'),
                                              (0.3, 'tab:green'),
                                              (0.5, 'yellow'),
                                              (0.7, 'tab:orange'),
                                              (0.9, 'tab:red'),
                                              (1, 'darkred')
                                             ], N=256)


def read_bin_score (fname):
    name_binID_score = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        if First:
            names = cols[4:]
            First = False
            continue
        _, chr, st, ed = cols[:4]
        st, ed = int(st), int(ed)
        scores = [float(score) for score in cols[4:]]
        try:
            chrnum = int(chr[3:])
        except:
            chrnum = chr[3:]
        binID = (chrnum, st, ed)
        for name, score in zip(names, scores):
            if name not in name_binID_score:
                name_binID_score[name] = {}
            name_binID_score[name][binID] = score
    return name_binID_score

            
def density_scatter(x , y, ax = None, sort = True, bins = 20, density = False, **kwargs )   :
    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d(x, y, bins = bins, density=density )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T ,
                 method = "splinef2d", bounds_error = False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        #print idx
        x, y, z = x[idx], y[idx], z[idx]

    img = ax.scatter( x, y, c=z, **kwargs )
    cbar = plt.colorbar(img)
    cbar.ax.tick_params(labelsize=5)
    #cbar = plt.colorbar(cm.ScalarMappable(norm = norm), ax=img)
    #cbar.ax.set_ylabel('Density')

    return ax
        
        
path = '/home/spark159/../../media/spark159/sw/'

field_ID_score1 = read_bin_score(path+"mCD8T_WT-NCP_sp_1kb_score.cn")
field_ID_score2 = read_bin_score(path+"mCD8T_inht-NCP_sp_1kb_score.cn")
field_ID_score3 = read_bin_score(path+"mCD8T_KO-NCP_sp_1kb_score.cn")

ID_score1 = field_ID_score1['mCD8T-WT-NCP-sp-8.bam']
ID_score2 = field_ID_score2['mCD8T-inht-NCP-sp-8.bam']
ID_score3 = field_ID_score3['mCD8T-KO-NCP-sp-8.bam']

IDs = list(set(ID_score1) & set(ID_score2) & set(ID_score3))

# get dscores
dscore_ID1 = []
dscore_ID2 = []
for ID in IDs:
    dscore1 = ID_score2[ID] - ID_score1[ID]
    dscore2 = ID_score3[ID] - ID_score1[ID]
    dscore_ID1.append([dscore1, ID])
    dscore_ID2.append([dscore2, ID])
dscore_ID1 = sorted(dscore_ID1)
dscore_ID2 = sorted(dscore_ID2)

# get dzscore after standardization 
ID_zscore1 = Z_score(ID_score1)
ID_zscore2 = Z_score(ID_score2)
ID_zscore3 = Z_score(ID_score3)

dzscore_ID1 = []
dzscore_ID2 = []
for ID in IDs:
    dzscore1 = ID_zscore2[ID] - ID_zscore1[ID]
    dzscore2 = ID_zscore3[ID] - ID_zscore1[ID]
    dzscore_ID1.append([dzscore1, ID])
    dzscore_ID2.append([dzscore2, ID])
dzscore_ID1 = sorted(dzscore_ID1)
dzscore_ID2 = sorted(dzscore_ID2)

# save top 20%
fname1 = 'Down_dscore_inht.cn'
fname2 = 'Up_dscore_inht.cn'
fname3 = 'Down_dscore_KO.cn'
fname4 = 'Up_dscore_KO.cn'
fname5 = 'Down_dzscore_inht.cn'
fname6 = 'Up_dzscore_inht.cn'
fname7 = 'Down_dzscore_KO.cn'
fname8 = 'Up_dzscore_KO.cn'

fract = 0.01
sample_size = int(fract*len(IDs))

write_cn({ID:dscore for dscore, ID in dscore_ID1[:sample_size]}, fname=fname1)
write_cn({ID:dscore for dscore, ID in dscore_ID1[::-1][:sample_size]}, fname=fname2)
write_cn({ID:dscore for dscore, ID in dscore_ID2[:sample_size]}, fname=fname3)
write_cn({ID:dscore for dscore, ID in dscore_ID2[::-1][:sample_size]}, fname=fname4)
write_cn({ID:dzscore for dzscore, ID in dzscore_ID1[:sample_size]}, fname=fname5)
write_cn({ID:dzscore for dzscore, ID in dzscore_ID1[::-1][:sample_size]}, fname=fname6)
write_cn({ID:dzscore for dzscore, ID in dzscore_ID2[:sample_size]}, fname=fname7)
write_cn({ID:dzscore for dzscore, ID in dzscore_ID2[::-1][:sample_size]}, fname=fname8)

subprocess.call('PATH=$PATH:/home/spark159/IGV_2.15.2/.', shell=True)
fnames = [fname1, fname2, fname3, fname4, fname5, fname6, fname7, fname8]
for fname in fnames:
    subprocess.call(['igvtools', 'toTDF', fname, fname.rsplit('.', 1)[0] + '.tdf', 'mm10'])










"""

fig = plt.figure()
plt.hist(ID_value1.values(), bins=100, alpha=0.5, label='WT')
plt.hist(ID_value2.values(), bins=100, alpha=0.5, label='+inht')
plt.hist(ID_value3.values(), bins=100, alpha=0.5, label='KO')
plt.xlim([-1, 6])
plt.legend()
plt.savefig("comp8.png", bbox_inches='tight')
#plt.show()
plt.close()

X, Y, Z = [], [], []
for ID in IDs:
    X.append(ID_value1[ID])
    Y.append(ID_value2[ID])
    Z.append(ID_value3[ID])

rX = statis.standardization(X)
rY = statis.standardization(Y)
rZ = statis.standardization(Z)

fig = plt.figure()
plt.plot(X, Y, ',')
#plt.show()
plt.close()

fig = plt.figure()
density_scatter(X, Y, bins = [20,20], s=3, cmap=pastel_jet, ax=plt.gca(), sort=False)
plt.plot([0, 5], [0, 5], 'k--', ms=2, alpha=0.7)
#plt.xlim([0,5])
#plt.ylim([0,5])
plt.xlabel("score (WT)")
plt.ylabel("score (+ODC inht)")
plt.title("1kb bins")
plt.savefig("WTVSinht.png", bbox_inches='tight')
#plt.show()
plt.close()

fig = plt.figure()
density_scatter(X, Z, bins = [20,20], s=3, cmap=pastel_jet, ax=plt.gca(), sort=False)
plt.plot([0, 5], [0, 5], 'k--', ms=2, alpha=0.7)
#plt.xlim([0,5])
#plt.ylim([0,5])
plt.xlabel("score (WT)")
plt.ylabel("score (ODC KO)")
plt.title("1kb bins")
plt.savefig("WTVSKO.png", bbox_inches='tight')
#plt.show()
plt.close()
"""
