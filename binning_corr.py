import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
import random
import matplotlib.cm as cm
from scipy.optimize import curve_fit
from sklearn import linear_model
from scipy.special import expit

# read annotation file
path = "./data/"
field_ID_value = load_file.read_tabular_file (path + "hg19_chr1_167win25step_anot.cn", mode='col', jump=10)

ID_pos = field_ID_value['PhysicalPosition']
ID_AT = field_ID_value['ATcontent']
ID_k27ac = field_ID_value['k27ac']
ID_score1 = field_ID_value["data/sp_spd_tests_detail/sp7"]
ID_score2 = field_ID_value["data/sp_spd_tests_detail/sp8"]

# binning parameters and choice of data
#names = ['data/sp_spd_tests_detail/sp7', 'ATcontent', 'CpGNumber', 'k27ac', 'k9ac', 'k4me3', 'k36me3', 'k9me2', 'k9me3', 'k27me3']
names = ['data/sp_spd_tests_detail/sp7', 'Compartments']

#step_size = 50000 # sampling step size

# temporal binning
#i = 20
#bin_size = int(step_size / i)
#blur_win = int(4*i + 1)
#bin_num = int(max(ID_pos.values()))/bin_size + 1

bin_size = int(50000)
bin_num = int(max(ID_pos.values()))/bin_size + 1
#blur_win = int(1)


# binning the annotation data
name_binID_mean = {}
name_binID_count = {}
for ID in ID_pos:
    binID = int(ID_pos[ID]) / int(bin_size)
    for name in names:
        if name not in name_binID_mean:
            name_binID_mean[name] = {}
        if name not in name_binID_count:
            name_binID_count[name] = {}
        try:
            value = field_ID_value[name][ID]
        except:
            continue
        if np.isnan(value):
            continue
        if binID not in name_binID_mean[name]:
            name_binID_mean[name][binID] = 0.0
        name_binID_mean[name][binID] += value
        if binID not in name_binID_count[name]:
            name_binID_count[name][binID] = 0
        name_binID_count[name][binID] += 1

for name in name_binID_mean:
    for binID in name_binID_mean[name]:
        name_binID_mean[name][binID] = float(name_binID_mean[name][binID]) / name_binID_count[name][binID]

# read RPKM file
if 'Gene density' in names or 'Gene activity' in names:
    geneID_field_values, field_geneID_values = load_file.read_GTF (path + "Homo_sapiens.GRCh37.87.gtf", chr_list=['chr1'], mode="both")
    geneID_RPKM = load_file.read_RPKM (path+"GSE63124_all_gene_raw_readcounts.txt", path+"Homo_sapiens.GRCh37.87.gtf", "chr1")

    geneID_pos = {}
    for geneID in geneID_field_values:
        pos = geneID_field_values[geneID]['TSS']
        geneID_pos[geneID] = pos

    if 'Gene density' in names:
        name_binID_mean['Gene density'] = {binID:0 for binID in range(bin_num)}
        name_binID_count['Gene density'] = {binID:0 for binID in range(bin_num)}

    if 'Gene activity' in names:
        name_binID_mean['Gene activity'] = {binID:0 for binID in range(bin_num)}
        name_binID_count['Gene activity'] = {binID:0 for binID in range(bin_num)}

    min_RPKM = min(geneID_RPKM.values())
    for geneID in geneID_pos:
        binID = int(geneID_pos[geneID]) / int(bin_size)
        if 'Gene density' in names:
            name_binID_mean['Gene density'][binID] += 1.0
        if 'Gene activity' in names:
            try:
                #name_binID_mean['Gene activity'][binID] += geneID_RPKM[geneID]
                name_binID_mean['Gene activity'][binID] += np.log2(geneID_RPKM[geneID] - min_RPKM + 1)
                #name_binID_mean['Gene activity'][binID] += np.log2(geneID_RPKM[geneID])
            except:
                name_binID_mean['Gene activity'][binID] += 0.0
                #name_binID_mean['Gene activity'][binID] += np.nan
            name_binID_count['Gene activity'][binID] += 1

    if 'Gene activity' in names:
        for binID in name_binID_mean['Gene activity']:
            if name_binID_count['Gene activity'][binID] <= 0:
                continue
            name_binID_mean['Gene activity'][binID] /= name_binID_count['Gene activity'][binID] 

# read A/B compartment annotation
if 'Compartments' in names:
    def read_eigenfile (fname, bin_size=1000000):
        eigen_list = []
        interval_list = []
        i = 0
        for line in open(fname):
            line = line.strip()
            if not line:
                continue
            try:
                value = float(line)
            except:
                value = np.nan
            st = i*bin_size
            ed = (i+1)*bin_size
            eigen_list.append(value)
            interval_list.append((st,ed))
            i +=1
        return eigen_list, interval_list

    fnames = ['eigen_WT_50kb.txt', 'eigen_CohesinKO_50kb.txt'] 
    for i in range(len(fnames)):
        fname = fnames[i]
        eigen_list, interval_list = read_eigenfile(path+fname, bin_size=50000)

        binID_eigens = {}
        for value, interval in zip(eigen_list, interval_list):
            st, ed = interval
            st_binID, ed_binID = st / bin_size, ed / bin_size
            for binID in range(st_binID, ed_binID):
                if binID not in binID_eigens:
                    binID_eigens[binID] = []
                binID_eigens[binID].append(value)

        binID_eigen = {binID:0 for binID in range(bin_num)}
        for binID in binID_eigen:
            try:
                binID_eigen[binID] = np.mean(binID_eigens[binID])
            except:
                continue
                #binID_eigen[binID] = np.nan
        name_binID_mean['Compartments' + str(i)] = binID_eigen

# get linear data array
name_sig = {}
for name in name_binID_mean:
    binID_mean = name_binID_mean[name]
    sig = []
    for binID in range(bin_num):
        try:
            sig.append(binID_mean[binID])
        except:
            sig.append(np.nan)
    #if name != 'Compartments':
    #    sig = statis.slow_moving_average2(sig, blur_win)
    name_sig[name] = sig

# plot A/B compartment and condensability
X = []
Y1, Y2 = [], []
for i in range(bin_num):
    x = name_sig['data/sp_spd_tests_detail/sp7'][i]
    y1, y2 = name_sig['Compartments0'][i], name_sig['Compartments1'][i]
    if np.isnan(x) or np.isnan(y1) or np.isnan(y2):
        continue
    X.append(x)
    Y1.append(y1)
    Y2.append(y2)

def func(x, a, b, c):
    return float(a) / (1 + np.exp(-b*(x-c))) - a*0.5

X_test = np.linspace(-1, 1, num=1000)
popt1, pcov1 = curve_fit(func, X, Y1)
popt2, pcov2 = curve_fit(func, X, Y2)
Y1_predict = func(X_test, *popt1)
Y2_predict = func(X_test, *popt2)

fig = plt.figure()
plt.plot(X, Y1, 'o', markersize=2, alpha=0.2)
plt.plot(X_test, Y1_predict, 'b', linewidth=2, label='WT')
plt.plot(X, Y2, 'o', markersize=2, alpha=0.2)
plt.plot(X_test, Y2_predict, 'r', linewidth=2, label='Cohesin KO')
leg = plt.legend()
#for lh in leg.legendHandles:
#    lh._legmarker.set_markersize(15)
#    lh._legmarker.set_alpha(1)
plt.xlim([-1,1])
plt.ylim([-0.025,0.03])
plt.xlabel("Condensability (A.U.)")
plt.ylabel("Eigenvector")
plt.title("50kb window")
plt.show()
plt.close()    

# make bedgraph file
def make_bedgraph (fname, binID_value, header=None):
    f = open(fname, 'w')
    if header == None:
        header = 'track type=bedGraph name="BedGraph Format" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20'
        print >> f, header
    for i in range(len(binID_value)):
        st, ed = i*bin_size, (i+1)*bin_size
        value = float(binID_value[i])
        if np.isnan(value):
            continue
        print >> f, "chr1\t" + "%d\t%d\t%f" % (st, ed, value)
    f.close()

binID_value = name_sig["data/sp_spd_tests_detail/sp7"]
make_bedgraph("condensabiltiy.bedgraph", binID_value)


    
