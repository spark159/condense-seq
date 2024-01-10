import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
import pickle
from sklearn import linear_model
from sklearn.neighbors import LocalOutlierFactor
from sklearn.covariance import EllipticEnvelope
import seaborn as sns
from scipy import stats
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
        #print idx
        x, y, z = x[idx], y[idx], z[idx]

    img = ax.scatter( x, y, c=z, **kwargs )
    cbar = plt.colorbar(img)
    cbar.ax.tick_params(labelsize=5)
    #cbar = plt.colorbar(cm.ScalarMappable(norm = norm), ax=img)
    #cbar.ax.set_ylabel('Density')

    return ax


## for H1
#name_dict = {"E1":"Polycomb repressed",
#             "E2":"Poised promoter",
#             "E3":"Weak promoter",
#             "E4":"Strong enhancer",
#             "E5":"Active promoter",
#             "E6":"Weak enhancer",
#             "E7":"Quiescence1",
#             "E8":"Quiescence2",
#             "E9":"Heterochromatin",
#             "E10":"Tx elongation",
#             "E11":"Weak Tx",
#             "E12":"Insulator"}

# state for H1
#states = ["Active promoter", "Weak promoter", "Poised promoter", "Strong enhancer", "Weak enhancer", "Tx elongation", "Weak Tx", "Insulator", "Polycomb repressed", "Heterochromatin", "Quiescence1", "Quiescence2"]


#state_intervals = read_chromHMM("H1_12_segments.bed", chr_target='chr1', change=name_dict)

# for GM12878
#name_dict = {"E1":"Polycomb repressed",
#             "E2":"Quiescence",
#             "E3":"Heterochromatin",
#             "E4":"Weak Tx",
#             "E5":"Tx elongation",
#             "E6":"Weak enhancer",
#             "E7":"Active enhancer",
#             "E8":"Strong enhancer",
#             "E9":"Active promoter",
#             "E10":"Weak promoter",
#             "E11":"Poised promoter",
#             "E12":"Insulator"}

# state for GM
#states = ["Active promoter", "Weak promoter", "Poised promoter", "Strong enhancer", "Active enhancer", "Weak enhancer", "Tx elongation", "Weak Tx", "Insulator", "Polycomb repressed", "Heterochromatin", "Quiescence"]

#state_intervals = read_chromHMM("GM12878_12_segments.bed", chr_target='chr1', change=name_dict)

# for mouse CD8 T cell
name_dict = {"E1":"Weak Tx",
             "E2":"Tx elongation",
             "E3":"Weak enhancer2",
             "E4":"Strong enhancer2",
             "E5":"Strong enhancer1",
             "E6":"Weak enhancer1",
             "E7":"Active promoter",
             "E8":"Poised promoter",
             "E9":"Polycomb repressed1",
             "E10":"Polycomb repressed2",
             "E11":"Quiescence",
             "E12":"Heterochromatin"}

# state for mouse CD8 T cell
states = ["Active promoter", "Poised promoter", "Strong enhancer1", "Strong enhancer2", "Weak enhancer1", "Weak enhancer2", "Tx elongation", "Weak Tx", "Polycomb repressed1", "Polycomb repressed2", "Heterochromatin", "Quiescence"]



def tuple_cmp (a, b):
    if a[0] <= b[0]:
        return -1
    return 1

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

def read_chromHMM(fname, chr_list=None, change=False):
    chr_state_intervals = {}
    for line in open(fname):
        cols = line.strip().split()
        chr, st, ed, state = cols[:4]
        if chr_list!=None and chr not in chr_list:
            continue
        if chr not in chr_state_intervals:
            chr_state_intervals[chr] = {}
        st, ed = int(st), int(ed)
        if change:
            state = change[state]
        if state not in chr_state_intervals[chr]:
            chr_state_intervals[chr][state] = []
        chr_state_intervals[chr][state].append((st,ed))
    return chr_state_intervals

# read binnum file and convert it to score
def read_binnum_file (fname):
    chr_binID_range = {}
    chr_binID_score = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        if First:
            binID = 0
            First = False
            continue
        cols = line.split()
        ID, chr, st, ed = cols[:4]

        if chr not in chr_binID_range:
            chr_binID_range[chr] = {}
            binID = 0
        chr_binID_range[chr][binID] = (int(st), int(ed))

        if cols[-1] == '*':
            binID +=1
            continue

        counts = [int(count) for count in cols[4:]]
        sup1, sup2 = counts[:2]
        control = counts[-1]

        if sup2 <=0 or control <=0:
            binID +=1
            continue
        
        score = -np.log(float(sup2)/control)

        if chr not in chr_binID_score:
            chr_binID_score[chr] = {}
            binID = 0
        chr_binID_score[chr][binID] = score
        binID +=1

    return chr_binID_range, chr_binID_score
        

# parameters
path = "/home/spark159/../../media/spark159/sw/"

#cell1, cell2 = "H1", "GM"
cell1, cell2 = "mCD8T","mCD8T"
#sample1, sample2 = "NCP", "NCP"
sample1, sample2 = "WT-NCP", "inht-NCP"
agent = "sp"
bin_size = 1000

bincount_fname1 = "_".join([cell1, sample1, agent, str(int(bin_size/1000.0)) + 'kb']) + "_num.cn"
bincount_fname2 = "_".join([cell2, sample2, agent, str(int(bin_size/1000.0)) + 'kb']) + "_num.cn"

# read bincount files
chr_binID_range1, chr_binID_score1 = read_binnum_file(path + bincount_fname1)
chr_binID_range2, chr_binID_score2 = read_binnum_file(path + bincount_fname2)

# check range match
assert chr_binID_range1 == chr_binID_range2
chr_binID_range = copy.deepcopy(chr_binID_range1)
del chr_binID_range1
del chr_binID_range2

# finalize common gene set
chr_binIDs = {}
total_count = 0
for chr in chr_binID_score1:
    binIDs1 = set(chr_binID_score1[chr])
    binIDs2 = set(chr_binID_score2[chr])
    binIDs = sorted(list(binIDs1 & binIDs2))
    chr_binIDs[chr] = binIDs
    total_count += len(binIDs)
print 'Total bin count:' + str(total_count)

# plot scatter plot
X, Y = [], []
for chr in chr_binIDs:
    for binID in chr_binIDs[chr]:
        X.append(chr_binID_score1[chr][binID])
        Y.append(chr_binID_score2[chr][binID])

X = statis.standardization(X)
Y = statis.standardization(Y)

fig = plt.figure()
#plt.plot(X, Y, 'k,', alpha=0.2, markersize=1.5)
density_scatter(X, Y, bins = [20,20], s=3, cmap=pastel_jet, ax=plt.gca(), sort=False)
#plt.plot([0, 7], [0, 7], 'k--', alpha=0.7)
plt.plot([-5, 5], [-5, 5], 'k--', alpha=0.7)
plt.title("Condensability (1kb)")
plt.xlabel(" ".join([cell1, sample1]))
plt.ylabel(" ".join([cell2, sample2]))
plt.savefig("VS".join([cell1+sample1, cell2+sample2]) + '.png')
#plt.show()
plt.close()

sys.exit(1)

# categorize according to the chromatin state
# read chromHMM file
chr_state_intervals = read_chromHMM(chromHMM_fname, change=name_dict)

chr_list = ['chr1']
chr_binID_state = {}
for chr in chr_list:
    print "categorizing %s data" % (chr)
    binID_states = {}
    state_intervals = chr_state_intervals[chr]
    for state in state_intervals:
        for interval in state_intervals[state]:
            st, ed = interval
            st_binID = int(st) / bin_size
            ed_binID = int(ed) / bin_size
            ed_binID = min(ed_binID, max(chr_binID_range[chr].keys()))
            for k in range(st_binID, ed_binID+1):
                bin_st = k*bin_size
                bin_ed = bin_st + bin_size
                value = min(bin_ed, ed) - max(bin_st, st)
                if value <= 0:
                    continue
                if k not in binID_states:
                    binID_states[k] = []
                binID_states[k].append((value, state))

    binID_state = {}
    for binID in binID_states:
        state = sorted(binID_states[binID], cmp=tuple_cmp, reverse=True)[0][1]
        binID_state[binID] = state

    if chr not in chr_binID_state:
        chr_binID_state[chr] = {}
    chr_binID_state[chr].update(binID_state)


state_tuples = {}
for chr in chr_list:
    for binID in chr_binID_state[chr]:
        state = chr_binID_state[chr][binID]
        try:
            x = chr_binID_score1[chr][binID]
            y = chr_binID_score2[chr][binID]
        except:
            continue
        if state not in state_tuples:
            state_tuples[state] = [[], []]
        state_tuples[state][0].append(x)
        state_tuples[state][1].append(y)


fig = plt.figure()
for state in states:
    X, Y = state_tuples[state][0], state_tuples[state][1]
    plt.plot(X, Y, '.', alpha=0.2, markersize=1.5, label=state)
plt.plot([0, 7], [0, 7], 'k--', alpha=0.7)
plt.title("Condensability (chromosome1, 1kb)")
plt.xlabel('Mouse CD8 T cell (WT)')
plt.ylabel('Mouse CD8 T cell (+ODC inhibitor)')
leg = plt.legend(loc='best', numpoints=1, prop={'size': 6})
for lh in leg.legendHandles:
    lh._legmarker.set_markersize(15)
    lh._legmarker.set_alpha(1)
plt.savefig('WTVSODCinht_chromHMM.png')
#plt.show()
plt.close()

for state in states:
    fig = plt.figure()
    X, Y = state_tuples[state][0], state_tuples[state][1]
    plt.plot(X, Y, '.', alpha=0.5, markersize=1.5)
    plt.plot([0, 7], [0, 7], 'k--', alpha=0.7)
    plt.title("Condensability (chromosome1, 1kb)")
    plt.xlabel('Mouse CD8 T cell (WT)')
    plt.ylabel('Mouse CD8 T cell (+ODC inhibitor)')
    plt.savefig('WTVSODCinht_' + state + '.png')
    #plt.show()
    plt.close()




sys.exit(1)


"""
# finalize common gene set
gIDs = list(set(gID_mscore1) & set(gID_mscore2) & set(gID_FPKM1.keys()) & set(gID_FPKM2.keys()))
print 'Total gene count:' + str(len(gIDs))

## standardization of scores
#scores1 = [gID_mscore1[gID] for gID in gIDs]
#scores2 = [gID_mscore2[gID] for gID in gIDs]
#score_mean1, score_std1 = np.mean(scores1), np.std(scores1)
#score_mean2, score_std2 = np.mean(scores2), np.std(scores2)

#for gID in gIDs:
#    gID_mscore1[gID] = float(gID_mscore1[gID] - score_mean1) / score_std1
#    gID_mscore2[gID] = float(gID_mscore2[gID] - score_mean2) / score_std2


# find Ensemble Gene ID for ESC marker genes / PC marker genes
ESC_gname_gIDs = {gname :[] for gname in ESC_gnames}
PC_gname_gIDs = {gname:[] for gname in PC_gnames}
mCD4Tcell_gname_gIDs = {gname:[] for gname in mCD4Tcell_gnames}
HOX_gname_gID = {}

for gID in gIDs:
    gname = gID_field_values[gID]['geneName'].upper()
    if gname.startswith('HOX'):
        HOX_gname_gID[gname] = gID
    try:
        ESC_gname_gIDs[gname].append(gID)
    except:
        pass
    try:
        PC_gname_gIDs[gname].append(gID)
    except:
        pass
    try:
        mCD4Tcell_gname_gIDs[gname].append(gID)
    except:
        pass
    
ESC_gID_gname = {}
for gname in ESC_gname_gIDs:
    if len(ESC_gname_gIDs[gname]) == 1:
        gID = ESC_gname_gIDs[gname][0]
        assert gID not in ESC_gID_gname
        ESC_gID_gname[gID] = gname

PC_gID_gname = {}
for gname in PC_gname_gIDs:
    if len(PC_gname_gIDs[gname]) == 1:
        gID = PC_gname_gIDs[gname][0]
        assert gID not in PC_gID_gname
        PC_gID_gname[gID] = gname

mCD4Tcell_gID_gname = {}
for gname in mCD4Tcell_gname_gIDs:
    if len(mCD4Tcell_gname_gIDs[gname]) == 1:
        gID = mCD4Tcell_gname_gIDs[gname][0]
        assert gID not in mCD4Tcell_gID_gname
        mCD4Tcell_gID_gname[gID] = gname


HOX_gID_gname = {}
for gname in HOX_gname_gID:
    gID = HOX_gname_gID[gname]
    HOX_gID_gname[gID] = gname

# find bivalent genes
gname_binfo, btype_gnames = read_bivalent("Bivalent_info.csv")
btype_gID_gname = {}
for gID in gIDs:
    gname = gID_field_values[gID]['geneName']
    try:
        btype = gname_binfo[gname]['btype']
        if btype not in btype_gID_gname:
            btype_gID_gname[btype] = {}
        assert gID not in btype_gID_gname[btype]
        btype_gID_gname[btype][gID] = gname
    except:
        continue


# set in-IDs and out-IDs
in_gIDs = list(set(gIDs) - set(ESC_gID_gname.keys())) # all genes except ESC marker
out_gIDs = ESC_gID_gname.keys() # ESC marker genes


# H1 FPKM vs GM FPKM
X, Y = [], []
for gID in gIDs:
    X.append(np.log2(1+gID_FPKM1[gID]))
    Y.append(np.log2(1+gID_FPKM2[gID]))
fig = plt.figure()
plt.plot(X, Y, '.')
plt.title('RNA-seq data comparison')
plt.xlabel('H1 hESC logFPKM')
#plt.ylabel('A38-41 hPanc logFPKM')
plt.ylabel('GM12878 logFPKM')
#plt.show()
plt.close()

# H1 FPKM vs H1 score
inX = [np.log2(1+gID_FPKM1[gID]) for gID in in_gIDs]
inY = [gID_mscore1[gID] for gID in in_gIDs]
outX = [np.log2(1+gID_FPKM1[gID]) for gID in out_gIDs]
outY = [gID_mscore1[gID] for gID in out_gIDs]
gnames = [ESC_gID_gname[gID] for gID in out_gIDs]

# get correaltion
X, Y = inX + outX, inY + outY
print 'Spearman corr: ', statis.get_spearman_corr(X, Y)

# polynomial fitting
#feature_list = [[x, x**2, x**3] for x in X]
#test_list = [[y] for y in Y]
#reg = linear_model.Ridge(alpha=0.5)
#reg.fit(feature_list, test_list)
#Xrange = np.linspace(min(X), max(X), num=100)
#Ypred = reg.predict([[x, x**2, x**3] for x in Xrange])
#Ypred = [value[0] for value in Ypred]

group_gIDs = statis.partition({gID:np.log2(1+gID_FPKM1[gID]) for gID in gIDs}, 500)
meanX, meanY = [], []
for i in range(len(group_gIDs)):
    meanX.append(np.median([np.log2(1+gID_FPKM1[gID]) for gID in group_gIDs[i]]))
    meanY.append(np.median([gID_mscore1[gID] for gID in group_gIDs[i]]))

#fig = plt.figure(figsize=(2.8,2.4))
fig = plt.figure()
plt.plot(inX, inY, color='tab:blue', marker=',', linestyle='none', markersize=3, alpha=0.3)
for x, y, gname in zip(outX, outY, gnames):
    #plt.plot(x, y, 'r.', markersize=3, alpha=1, zorder=10)
    #if gname in ESC_tf_cores:
    if gname in ESC_gnames:
        if gname in ESC_tf_cores:
            plt.plot(x, y, 'r.', markersize=3, alpha=1, zorder=10)
            plt.annotate(gname, (x, y), color='red', zorder=40, size=8, weight='bold')
        else:
            plt.plot(x, y, 'k.', markersize=3, alpha=1, zorder=10)
            plt.annotate(gname, (x, y), color='black', zorder=40, size=8)
#sns.kdeplot(data=[[X[i], Y[i]] for i in range(len(X))])
#plt.plot(Xrange, Ypred, 'k--', alpha=0.7)
plt.plot(meanX, meanY, 'k--', alpha=0.7)
plt.title('H1-hESC (Near TSS)', fontsize=12)
plt.xlabel('Gene expression (logFPKM)', fontsize=12)
plt.ylabel('Condensabiltiy (A.U.)', fontsize=12)
plt.gca().tick_params(axis='both', which='major', labelsize=8)
plt.gca().tick_params(axis='both', which='minor', labelsize=8)
plt.ylim([-2.5, 2.5])
#plt.savefig('nearTSScondvsExpreesion_hESC.svg', format='svg', bbox_inches='tight')
plt.savefig('nearTSScondvsExpreesion_hESC.png', bbox_inches='tight')
#plt.show()
plt.close()

# Panc FPKM vs GM score
inX = [np.log2(1+gID_FPKM2[gID]) for gID in in_gIDs]
inY = [gID_mscore2[gID] for gID in in_gIDs]
outX = [np.log2(1+gID_FPKM2[gID]) for gID in out_gIDs]
outY = [gID_mscore2[gID] for gID in out_gIDs]
gnames = [ESC_gID_gname[gID] for gID in out_gIDs]

X, Y = inX + outX, inY + outY
print 'Spearman corr: ', statis.get_spearman_corr(X, Y)

# polynomial fitting
#feature_list = [[x, x**2, x**3] for x in X]
#test_list = [[y] for y in Y]
#reg = linear_model.Ridge(alpha=0.5)
#reg.fit(feature_list, test_list)
#Xrange = np.linspace(min(X), max(X), num=100)
#Ypred = reg.predict([[x, x**2, x**3] for x in Xrange])
#Ypred = [value[0] for value in Ypred]

group_gIDs = statis.partition({gID:np.log2(1+gID_FPKM2[gID]) for gID in gIDs}, 500)
meanX, meanY = [], []
for i in range(len(group_gIDs)):
    meanX.append(np.median([np.log2(1+gID_FPKM2[gID]) for gID in group_gIDs[i]]))
    meanY.append(np.median([gID_mscore2[gID] for gID in group_gIDs[i]]))

fig = plt.figure()
plt.plot(inX, inY, color='tab:blue', marker=',', linestyle='none', markersize=3, alpha=0.3)
for x, y, gname in zip(outX, outY, gnames):
    #plt.plot(x, y, 'r.', markersize=3, alpha=1, zorder=10)
    #if gname in ESC_tf_cores:
    if gname in ESC_gnames:
        if gname in ESC_tf_cores:
            plt.plot(x, y, 'r.', markersize=3, alpha=1, zorder=10)
            plt.annotate(gname, (x, y), color='red', zorder=40, size=8, weight='bold')
        else:
            plt.plot(x, y, 'k.', markersize=3, alpha=1, zorder=10)
            plt.annotate(gname, (x, y), color='black', zorder=40, size=8)
#plt.plot(Xrange, Ypred, 'k--', alpha=0.7)
plt.plot(meanX, meanY, 'k--', alpha=0.7)
#plt.title('A38-41 hPanc')
#plt.title('A38-41 hPanc (Near TSS)', fontsize=12)
plt.title('GM12878 (Near TSS)', fontsize=12)
plt.xlabel('Gene expression (logFPKM)', fontsize=12)
plt.ylabel('Condensabiltiy (A.U.)', fontsize=12)
plt.gca().tick_params(axis='both', which='major', labelsize=8)
plt.gca().tick_params(axis='both', which='minor', labelsize=8)
plt.ylim([-3.5, 3.5])
#plt.savefig('nearTSScondvsExpreesion_hPanc.svg', format='svg', bbox_inches='tight')
#plt.savefig('nearTSScondvsExpreesion_GM.svg', format='svg', bbox_inches='tight')
plt.savefig('nearTSScondvsExpreesion_GM.png', bbox_inches='tight')
#plt.show()
plt.close()


## H1 score VS GM score
X, Y = [], []
C = []
for gID in gIDs:
    X.append(gID_mscore1[gID])
    Y.append(gID_mscore2[gID])
    C.append(np.log2(1+gID_FPKM2[gID]) - np.log2(1+gID_FPKM1[gID]))
C = stats.zscore(C)

outX = [gID_mscore1[gID] for gID in out_gIDs]
outY = [gID_mscore2[gID] for gID in out_gIDs]
gnames = [ESC_gID_gname[gID] for gID in out_gIDs]


# draw all genes

# custom diverging colormap with white background
pastel_jet = LinearSegmentedColormap.from_list('white_viridis',
                                             [(0, 'darkblue'),
                                              (0.1, 'blue'),
                                              (0.2, 'tab:blue'),
                                              (0.4, 'tab:cyan'),
                                              (0.5, 'ivory'),
                                              (0.6, 'tab:orange'),
                                              (0.8, 'tab:red'),
                                              (0.9, 'red'),
                                              (1, 'darkred')
                                             ], N=256)


fig = plt.figure()
#plt.scatter(X, Y, c=C, cmap='Spectral', vmin=-3, vmax=3, alpha=0.3, s=3)
#plt.scatter(X, Y, c=C, cmap=pastel_jet, vmin=-5, vmax=5, alpha=0.3, s=2)
plt.plot(X, Y, 'k.', alpha=0.2, markersize=1.5)
        
for gID in mCD4Tcell_gID_gname:
    gname = mCD4Tcell_gID_gname[gID]
    x, y = gID_mscore1[gID], gID_mscore2[gID]
    #plt.plot(x, y, 'kx', markersize=5, alpha=1, zorder=10, mew=1.5)
    #plt.annotate(gname, (x, y), color='black', zorder=40, size=8, weight='bold')

#plt.plot([min(X), max(X)], [min(Y), max(Y)], 'k--', alpha=0.7)
#plt.plot([-2.5, 2], [-3.5, 3.5], 'k--', alpha=0.5)
#plt.plot([-4, 3.5], [-4, 3.5], 'r--', alpha=0.7)
plt.plot([0, 7], [0, 7], 'r--', alpha=0.7)
plt.title("Condensability near TSS (5kb)")
plt.xlabel('Mouse CD8 T cell (WT)')
#plt.ylabel('A38-41 hPanc')
plt.ylabel('Mouse CD8 T cell (+ODC inhibitor)')
#plt.ylim([-2.5, 2.5])
#plt.xlim([-3, 2.5])
#plt.ylim([-3, 2.5])
#plt.xlim([-4, 3.5])
#plt.ylim([-4, 3.5])
plt.xlim([0, 7])
plt.ylim([0, 7])

#cbar = plt.colorbar()
#cbar.ax.set_ylabel('Gene expression (A38-41 hPanc - H1 hESC)', rotation=-90, va="bottom")
#cbar.ax.set_ylabel('Gene expression (GM12878 - H1 hESC)', rotation=-90, va="bottom")
#cbar.ax.set_ylabel('Gene expression (+inht - WT)', rotation=-90, va="bottom")
#plt.savefig('A38VSH1hESC_all.png')
plt.savefig('WTVSODCinht_all.png')
plt.show()
plt.close()


# score difference vs mean gene expression (MA-plot)
X, Y = [] ,[]
gID_dscore = {}
for gID in gIDs:
    #mean_FPKM = 0.5*(np.log(1+gID_FPKM1[gID]) + np.log(1+gID_FPKM2[gID]))
    #mean_FPKM = 0.5*(gID_FPKM1[gID] + gID_FPKM2[gID])
    mean_score = 0.5*(gID_mscore1[gID] + gID_mscore2[gID])
    score_diff = gID_mscore2[gID] - gID_mscore1[gID]
    X.append(mean_score)
    Y.append(score_diff)
    gID_dscore[gID] = score_diff

fig = plt.figure()
plt.plot(X, Y, 'k.', markersize=1.5, alpha=0.2)
#plt.scatter(X, Y, c=C, cmap=pastel_jet, vmin=-5, vmax=5, alpha=0.3, s=2)
plt.axhline(y=0, linestyle='--', color='r')
plt.title("Condensability mean vs difference")
plt.xlabel('mean condensability')
#plt.ylabel('A38-41 hPanc - H1 hESC')
plt.ylabel('+ODC inht - WT')
#cbar = plt.colorbar()
#cbar.ax.set_ylabel('Gene expression (+inht - WT)', rotation=-90, va="bottom")
plt.savefig('MA_plot_all.png')
#plt.show()
plt.close()



# draw only bivalent genes
for btype in btype_gID_gname:
    Bt_gID_gname = btype_gID_gname[btype]
    
    fig = plt.figure()
    plt.scatter(X, Y, c='lightgrey', alpha=0.05, s=1)

    Bt_X, Bt_Y = [], []
    Bt_C = []
    for gID in btype_gID_gname[btype]:
        gname = btype_gID_gname[btype][gID]
        x, y = gID_mscore1[gID], gID_mscore2[gID]
        Bt_X.append(x)
        Bt_Y.append(y)
        Bt_C.append(np.log2(1+gID_FPKM2[gID]) - np.log2(1+gID_FPKM1[gID]))
    plt.scatter(Bt_X, Bt_Y, c=Bt_C, s=3, alpha=0.5, zorder=10, vmin=-3, vmax=3,  cmap='bwr')

    #plt.plot([min(X), max(X)], [min(Y), max(Y)], 'k--', alpha=0.7)
    plt.plot([-2.5, 2], [-3.5, 3.5], 'k--', alpha=0.5)
    plt.title("Condensability near TSS (5kb) %s" % (btype))
    plt.xlabel('H1 hESC')
    plt.ylabel('GM12878')
    plt.xlim([-2.5, 2])
    plt.ylim([-3.5, 3.5])
    #plt.xlim([-2.5, 2.5])
    #plt.ylim([-2.5, 2.5])
    cbar = plt.colorbar()
    #cbar.ax.set_ylabel('Gene expression (A38-41 hPanc - H1 hESC)', rotation=-90, va="bottom")
    cbar.ax.set_ylabel('Gene expression (GM12878 - H1 hESC)', rotation=-90, va="bottom")
    #plt.savefig('A38VSH1hESC_%s.png' % (btype))
    plt.savefig('GMVSH1hESC_%s.png' % (btype))
    #plt.show()
    plt.close()


# draw only HOX genes
fig = plt.figure()
plt.scatter(X, Y, c='lightgrey', alpha=0.05, s=1)

HOX_X, HOX_Y = [], []
HOX_C = []
for gID in HOX_gID_gname:
    gname = HOX_gID_gname[gID]
    x, y = gID_mscore1[gID], gID_mscore2[gID]
    HOX_X.append(x)
    HOX_Y.append(y)
    HOX_C.append(np.log2(1+gID_FPKM2[gID]) - np.log2(1+gID_FPKM1[gID]))
plt.scatter(HOX_X, HOX_Y, c=HOX_C, s=12, alpha=1, zorder=10, vmin=-3, vmax=3, edgecolor='k', linewidth=0.5, cmap='bwr')

#plt.plot([min(X), max(X)], [min(Y), max(Y)], 'k--', alpha=0.7)
plt.plot([-2.5, 2], [-3.5, 3.5], 'k--', alpha=0.5)
plt.title("Condensability near TSS (5kb) only HOX genes")
plt.xlabel('H1 hESC')
#plt.ylabel('A38-41 hPanc')
plt.ylabel('GM12878')
#plt.xlim([-2.5, 2.5])
#plt.ylim([-2.5, 2.5])
plt.xlim([-2.5, 2])
plt.ylim([-3.5, 3.5])
cbar = plt.colorbar()
#cbar.ax.set_ylabel('Gene expression (A38-41 hPanc - H1 hESC)', rotation=-90, va="bottom")
#plt.savefig('A38VSH1hESC_HOX.png')
cbar.ax.set_ylabel('Gene expression (GM12878 - H1 hESC)', rotation=-90, va="bottom")
plt.savefig('GMVSH1hESC_HOX.png')
#plt.show()
plt.close()


# draw only ES markers
fig = plt.figure()
plt.scatter(X, Y, c='lightgrey', alpha=0.05, s=1)

ESC_X, ESC_Y = [], []
ESC_C = []
for gID in ESC_gID_gname:
    gname = ESC_gID_gname[gID]
    x, y = gID_mscore1[gID], gID_mscore2[gID]
    ESC_X.append(x)
    ESC_Y.append(y)
    ESC_C.append(np.log2(1+gID_FPKM2[gID]) - np.log2(1+gID_FPKM1[gID]))
plt.scatter(ESC_X, ESC_Y, c=ESC_C, s=12, alpha=1, zorder=10, vmin=-3, vmax=3, edgecolor='k', linewidth=0.5, cmap='bwr')

#plt.plot([min(X), max(X)], [min(Y), max(Y)], 'k--', alpha=0.7)
plt.plot([-2.5, 2], [-3.5, 3.5], 'k--', alpha=0.5)
plt.title("Condensability near TSS (5kb) only ES markers")
plt.xlabel('H1 hESC')
#plt.ylabel('A38-41 hPanc')
plt.ylabel('GM12878')
#plt.xlim([-2.5, 2.5])
#plt.ylim([-2.5, 2.5])
plt.xlim([-2.5, 2])
plt.ylim([-3.5, 3.5])
cbar = plt.colorbar()
#cbar.ax.set_ylabel('Gene expression (A38-41 hPanc - H1 hESC)', rotation=-90, va="bottom")
#plt.savefig('A38VSH1hESC_ESm.png')
cbar.ax.set_ylabel('Gene expression (GM12878 - H1 hESC)', rotation=-90, va="bottom")
plt.savefig('GMVSH1hESC_ESm.png')
#plt.show()
plt.close()





### H1 score VS Panc score
#inX = [gID_mscore1[gID] for gID in in_gIDs]
#inY = [gID_mscore2[gID] for gID in in_gIDs]
#outX = [gID_mscore1[gID] for gID in out_gIDs]
#outY = [gID_mscore2[gID] for gID in out_gIDs]
#gnames = [ESC_gID_gname[gID] for gID in out_gIDs]
#
#X = inX + outX
#Y = inY + outY
#
#fig = plt.figure()
#plt.plot(inX, inY, color='tab:blue', marker=',', linestyle='none', markersize=3, alpha=0.3)
#for x, y, gname in zip(outX, outY, gnames):
#    plt.plot(x, y, 'r.', markersize=3, alpha=1, zorder=10)
#    if gname in ESC_tf_cores:
#        plt.annotate(gname, (x, y), color='red', zorder=40, size=8, weight='bold')
#plt.plot([min(X), max(X)], [min(Y), max(Y)], 'k--', alpha=0.7)
#plt.title("Condensability near TSS (5kb)")
#plt.xlabel('H1 hESC')
#plt.ylabel('A38-41 hPanc')
#plt.xlim([-2.5, 2.5])
#plt.ylim([-2.5, 2.5])
#plt.show()
#plt.close()

#sys.exit(1)

#data_list = []
#for gID in gIDs:
#    data = [gID_mscore1[gID], gID_mscore2[gID]]
#    data_list.append(data)

#clf = EllipticEnvelope()
#outcheck = clf.fit_predict(data_list)


#fig = plt.figure()
#for i in range(len(data_list)):
#    data = data_list[i]
#    if outcheck[i] < 0:
#        #mark = 'r.'
#        mark = 'k.'
#    else:
#        mark = 'k.'
#    plt.plot([data[0]], [data[1]], mark, markersize=3, alpha=0.5)

#plt.title("Condensability near TSS (5kb)")
#plt.xlabel('H1 hESC')
#plt.ylabel('A38-41 hPanc')
#plt.show()
#plt.close()

#sys.exit(1)
    

#inX, inY = [], []
#outX, outY = [], []
#for i in range(len(X), len(Y)):
#    if outcheck[i] < 0:
#        outX.append(X[i])
#        outY.append(Y[i])
#    else:
#        inX.append(X[i])
#        inY.append(Y[i])

#reg = linear_model.Ridge(alpha=0.5)
#reg.fit (feature_list, test_list)
#Ypred = reg.predict(feature_list)
#Ypred = [ value[0] for value in Ypred]


#fig = plt.figure()
#plt.plot(X, Y, 'k.', markersize=3, alpha=0.5)
#plt.plot(X, Ypred, 'r--')
#plt.title("Condensability near TSS (5kb)")
#plt.xlabel('H1 hESC')
#plt.ylabel('A38-41 hPanc')
#plt.show()
#plt.close()

# gene expression difference VS score differences
inX = [np.log2(1+gID_FPKM2[gID]) - np.log2(1+gID_FPKM1[gID]) for gID in in_gIDs]
inY = [gID_mscore2[gID] - gID_mscore1[gID] for gID in in_gIDs]
outX = [np.log2(1+gID_FPKM2[gID]) - np.log2(1+gID_FPKM1[gID]) for gID in out_gIDs]
outY = [gID_mscore2[gID] - gID_mscore1[gID] for gID in out_gIDs]
gnames = [ESC_gID_gname[gID] for gID in out_gIDs]

fig = plt.figure()
plt.plot(inX, inY, color='tab:blue', marker='.', linestyle='none', markersize=2, alpha=0.2)
for x, y, gname in zip(outX, outY, gnames):
    plt.plot(x, y, 'r.', markersize=2, alpha=1, zorder=10)
    if gname in ESC_tf_cores:
        plt.annotate(gname, (x, y), color='red', zorder=40, size=8, weight='bold')
#plt.xlabel('logFPKM (A38-41 hPanc - H1 hESC)')
#plt.ylabel('Condensability (A38-41 hPanc - H1 hESC)')
plt.xlabel('logFPKM (GM12878 - H1 hESC)')
plt.ylabel('Condensability (GM12878 - H1 hESC)')
plt.axvline(x=0, color='k', linestyle='--', alpha=0.7)
plt.axhline(y=0, color='k', linestyle='--', alpha=0.7)
#plt.xlim([-2.5, 2.5])
#plt.ylim([-2.5, 2.5])
#plt.show()
plt.close()

sys.exit(1)

# score difference vs mean gene expression (MA-plot)
X, Y = [] ,[]
gID_dscore = {}
for gID in gIDs:
    #mean_FPKM = 0.5*(np.log(1+gID_FPKM1[gID]) + np.log(1+gID_FPKM2[gID]))
    mean_FPKM = 0.5*(gID_FPKM1[gID] + gID_FPKM2[gID])
    score_diff = gID_mscore2[gID] - gID_mscore1[gID]
    X.append(np.log2(1+mean_FPKM))
    Y.append(score_diff)
    gID_dscore[gID] = score_diff

fig = plt.figure()
plt.plot(X, Y, 'k.', markersize=3, alpha=0.3)
plt.axhline(y=0, linestyle='--', color='r')
plt.title("Condensability difference vs mean gene expression")
plt.xlabel('mean log2(1+FPKM))')
#plt.ylabel('A38-41 hPanc - H1 hESC')
plt.ylabel('GM12878 - H1 hESC')
plt.show()
plt.close()


# score difference histogram and partitions
# Partition by score
med = np.median(gID_dscore.values())
std = np.std(gID_dscore.values())
lines = [med-0.5*std-i*std for i in range(3)] + [med+0.5*std+i*std for i in range(3)]
lines = sorted(lines)
p_num = len(lines)+1
p_range = []
for i in range(p_num):
    if i == 0:
        st = -np.inf
        ed = lines[i]
    elif i == p_num-1:
        st = ed
        ed = np.inf
    else:
        st = ed
        ed = lines[i]
    p_range.append((st, ed))
fig = plt.figure()
plt.hist(gID_dscore.values(), bins=100)
for line in lines:
    plt.axvline(x=line, color='k', linestyle='--')
num_rom = {1:'I', 2:'II', 3:'III', 4:'IV', 5:'V', 6:'VI', 7:'VII'}
for i in range(p_num):
    st, ed = p_range[i]
    if i == 0:
        x = np.mean([-3, ed])
    elif i == p_num - 1:
        x = np.mean([st, 3])
    else:
        x = np.mean([st, ed])
    plt.text(x, 10000, num_rom[i+1], fontsize=20, va='center', ha='center')
plt.xlim([-3,3])
#plt.title("Chromosome1")
plt.title("All genes")
#plt.xlabel("Condensability difference (A38-41 hPanc - H1 hESC)")
plt.xlabel("Condensability difference (GM12878 - H1 hESC)")
plt.ylabel("Gene Counts")
#plt.savefig("partition_hist.png")
plt.show()
plt.close()

p_IDs = [[] for i in range(p_num)]
for ID in gID_dscore:
    dscore = gID_dscore[ID]
    for i in range(p_num):
        st, ed = p_range[i]
        if dscore >= st and dscore < ed:
            break
    p_IDs[i].append(ID)

for i in range(len(p_IDs)):
    f = open("output_" + str(i) + ".txt", "w")
    for ID in p_IDs[i]:
        print >> f, ID
    f.close()
"""
