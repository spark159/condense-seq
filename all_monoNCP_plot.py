import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import math
import Interval_dict
from scipy.stats import gaussian_kde
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn
from matplotlib.colors import LinearSegmentedColormap

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
        x, y, z = x[idx], y[idx], z[idx]

    img = ax.scatter( x, y, c=z, **kwargs )
    cbar = plt.colorbar(img)
    cbar.ax.tick_params(labelsize=5)
    #cbar = plt.colorbar(cm.ScalarMappable(norm = norm), ax=img)
    #cbar.ax.set_ylabel('Density')

    return ax

#path = "/home/spark159/../../media/spark159/sw/sp_spd_tests_detail/"
#path='./data/'
path = ""

#GC/me/PTM analysis
#ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+ "hg19_chr1_1001win501step_anot.cn")
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+ "H1_NCP_sp_chr1_anot.cn")

ID_score1 = name_ID_value["work/2021_06_07_H1_sp_detail/H1-NCP-sp-4"]
ID_score2 = name_ID_value["work/2021_06_07_H1_sp_detail/H1-NCP-sp-8"]
ID_AT = name_ID_value['ATcontent']

for ID in ID_AT:
    ID_AT[ID] = 100*ID_AT[ID]

graphics.Scatter_plot(ID_AT, ID_score1, note='H1-NCP-sp-4', ylim=[-2.5, 3])
graphics.Scatter_plot(ID_AT, ID_score2, note='H1-NCP-sp-8', ylim=[-2.5, 3])

#sys.exit(1)


ID_score = ID_score1
X, Y = [], []
xvalue_scores = {}
for ID in ID_score1.keys():
    xvalue = ID_AT[ID]
    score = ID_score2[ID]
    X.append(xvalue)
    Y.append(score)
    if xvalue not in xvalue_scores:
        xvalue_scores[xvalue] = []
    xvalue_scores[xvalue].append(score)

Xmean, Ymean, Yerr = [], [], []
for xvalue in xvalue_scores:
    if len(xvalue_scores[xvalue]) <= 1:
        continue
    Xmean.append(xvalue)
    Ymean.append(np.mean(xvalue_scores[xvalue]))
    Yerr.append(np.std(xvalue_scores[xvalue]/np.sqrt(len(xvalue_scores[xvalue]))))

X, Y = np.asarray(X), np.asarray(Y)

# plot density scatter
fig = plt.figure(figsize=(3, 2))
density_scatter(X, Y, bins = [20,20], s=3, cmap=pastel_jet, ax=plt.gca())
plt.errorbar(Xmean, Ymean, yerr=Yerr, fmt='k.', markersize=1, elinewidth=0.5)
plt.plot(Xmean, Ymean,'k.', markersize=1)
plt.xlabel("AT content (%)", fontsize=8)
plt.ylabel("Condensability (A.U.)", fontsize=8)
plt.gca().tick_params(axis='both', which='major', labelsize=5)
plt.gca().tick_params(axis='both', which='minor', labelsize=5)
plt.xlim([0, 100])
plt.ylim([-2.5, 3])
plt.savefig("ATvsCondensability_scatter.png", bbox_inches='tight', dpi=500)
#plt.savefig("ATvsCondensability_scatter.svg", format='svg', bbox_inches='tight')
plt.close()
sys.exit(1)


# plot density scatter by using kde
# Calculate the point density
#X, Y = np.asarray(X), np.asarray(Y)
#XY = np.vstack([X,Y])
#Z = gaussian_kde(XY)(XY)

# Sort the points by density, so that the densest points are plotted last
#order = np.argsort(Z)
#X, Y, Z = X[order], Y[order], Z[order]

#fig = plt.figure()
#plt.scatter(X, Y, c=Z, s=5, cmap='jet', edgecolor='', alpha=0.1)
#plt.errorbar(Xmean, Ymean, yerr=Yerr, fmt='k.')
#plt.plot(Xmean, Ymean,'k.')
#plt.xlabel("AT content (%)")
#plt.ylabel("Condensability (A.U.)")
#plt.ylim([-2.5, 3])
#plt.savefig("ATvsCondensability_scatter.png")
#plt.show()
#plt.close()

#with open("X.pickle", "wb") as f:
#    pickle.dump(X, f)
#with open("Y.pickle", "wb") as f:
#    pickle.dump(Y, f)
#with open("Z.pickle", "wb") as f:
#    pickle.dump(Z, f)


sys.exit(1)

#graphics.Scatter_plot (ID_AT, ID_score1, ylim=[-3, 3.5], note='raw')

"""
#divide into domains
gID_field_values, field_gID_values = load_file.read_GTF ("data/Homo_sapiens.GRCh37.87.gtf", "chr1", mode="both")
x
gID_ginterval = {}
for gID in gID_field_values.keys():
    TSS = gID_field_values[gID]['TSS']
    TTS = gID_field_values[gID]['TTS']
    strand = gID_field_values[gID]['strand']
    interval = (TSS-500, TSS+500)
    #if strand == '+':
        #interval = (TSS, TTS)
        #interval = (TSS-250, TSS)
        #interval = (TSS, TSS+2500)
    #else:
        #interval = (TTS, TSS)
        #interval = (TSS, TSS+250)
        #interval = (TSS-2500, TSS)
    gID_ginterval[gID] = interval

ginterval_dict = Interval_dict.double_hash(gID_ginterval, 10000, 250000000)

temp_name_ID_value = {}
for ID in ID_pos:
    pos = ID_pos[ID]
    gIDs = ginterval_dict.find(pos)
    if not gIDs:
        continue
    for name in name_ID_value:
        if name not in temp_name_ID_value:
            temp_name_ID_value[name] = {}
        temp_name_ID_value[name][ID] = name_ID_value[name][ID]

name_ID_value = temp_name_ID_value
"""

#ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("hg19_chr1_171_anot.cn")
ID_score1 = name_ID_value['data/sp_spd_tests_detail/sp7']
ID_score2 = name_ID_value['data/sp_spd_tests_detail/sp8']
ID_AT = name_ID_value['ATcontent']
ID_CpG = name_ID_value['CpGNumber']
ID_me = name_ID_value['meGCNumber']

ID_chip_list, chip_names = [], []
for name in name_ID_value:
    if name.startswith('k'):
        ID_value = name_ID_value[name]
        chip_names.append(name)
        ID_chip_list.append(ID_value)

for ID in ID_AT:
    ID_AT[ID] = ID_AT[ID]*100

"""
ID_meCpGfrac = {}
ID_meGCfrac = {}
for ID in ID_CpG:
    CpG = ID_CpG[ID]
    me = ID_me[ID]
    GC = 100-ID_AT[ID]
    if CpG <= 0:
        value = 0.0
    else:
        float(me)/(171*CpG)
    ID_meCpGfrac[ID] = value
    ID_meGCfrac[ID] = float(me)/GC
"""

#new_ID_score2 = statis.neutralize_score_by_target(ID_score2, ID_AT)

#graphics.Scatter_plot(ID_AT, ID_score1, note='sp6')
#graphics.Scatter_plot(ID_AT, ID_score2, note='sp10')

#frac=[(4**i) for i in range(1,11)]
frac=[1 for i in range(10)]
frac = frac[::-1]


#graphics.PartitionMeanplot(ID_AT, ID_score1, ID_me, frac, note="sp6_me")
#graphics.PartitionMeanplot(ID_AT, ID_score2, ID_me, frac, note="sp7_me")
#graphics.PartitionScatterplot(ID_AT, ID_score1, ID_me, frac, note="sp6_me")
graphics.PartitionBoxplot(ID_score1, ID_CpG, frac, xlabel='CpG number', note="sp6_CpG")

for i in range(len(ID_chip_list)):
    ID_chip = ID_chip_list[i]
    name = chip_names[i]
    #graphics.PartitionMeanplot(ID_AT, ID_score1, ID_chip, frac, note="sp6_" + name)
    #graphics.PartitionMeanplot(ID_AT, ID_score2, ID_chip, frac, note="sp7_" + name)
    #graphics.PartitionScatterplot(ID_AT, ID_score1, ID_chip, frac, note="sp6_" + name)
    graphics.PartitionBoxplot(ID_score1, ID_chip, frac, xlabel=name, note="sp6_" + name)
