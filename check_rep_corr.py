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
import scipy
import seaborn as sns
import matplotlib as mpl
from scipy.stats import gaussian_kde
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn
from matplotlib.colors import LinearSegmentedColormap

def read_bin_num (fname, skip_star=False, chr_choices=None):
    First = True
    ID_pos = {}
    name_ID_count = {}
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        if First:
            names = []
            for name in cols[4:-1]:
                #agent = name.rsplit('.', 1)[0].split('-')[-2]
                #tnum = int(name.rsplit('.', 1)[0].split('-')[-1])
                names.append(name.rsplit('.', 1)[0]) 
            First = False
            continue
        if cols[-1] == '*':
            if skip_star:
                continue
            cols = cols[:-1]
        ID, chr, start, end = cols[:4]
        start, end = int(start), int(end)
        if chr_choices!=None and chr not in chr_choices:
            continue
        pos = int(float(start + end)/2)
        assert ID not in ID_pos
        ID_pos[ID] = pos
        counts = [int(count) for count in cols[4:]]
        for name, count in zip(names, counts):
            if name not in name_ID_count:
                name_ID_count[name] = {}
            assert ID not in name_ID_count[name]
            name_ID_count[name][ID] = count
    return ID_pos, name_ID_count



def read_bin_file (fname, chr_choices=None):
    First = True
    ID_pos = {}
    name_ID_value = {}
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        if First:
            names = []
            for name in cols[4:]:
                #agent = name.rsplit('.', 1)[0].split('-')[-2]
                #tnum = int(name.rsplit('.', 1)[0].split('-')[-1])
                #names.append(tnum)
                names.append(name.rsplit('.', 1)[0])
            First = False
            continue
        ID, chr, start, end = cols[:4]
        start, end = int(start), int(end)
        if chr_choices!=None and chr not in chr_choices:
            continue
        pos = int(float(start + end)/2)
        ID = (chr, start, end)
        assert ID not in ID_pos
        ID_pos[ID] = pos
        values = [float(value) for value in cols[4:]]
        for name, value in zip(names, values):
            if name not in name_ID_value:
                name_ID_value[name] = {}
            assert ID not in name_ID_value[name]
            name_ID_value[name][ID] = value
    return ID_pos, name_ID_value


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

def density_scatter(x , y, ax = None, sort = False, bins = 20, density = False, **kwargs )   :
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





# parameters
# set path
#path = ""
#path = "/home/spark159/../../media/spark159/sw/"
path = "/home/spark159/../../storage/"

# pathway for each replicates
rep_path = {1:"/home/spark159/../../storage/",
            2:"/home/spark159/../../storage/replicates/"}

agent_fullname = {'sp':'Spermine(4+)',
                  'spd':'Spermidine(3+)',
                  'CoH':'Cobalt Hexammine(3+)',
                  'PEG':'PEG 8000',
                  'HP1a':'HP1$\\alpha$',
                  'HP1bSUV':'HP1$\\beta$+tSUV',
                  'LKH':'Linker histone1',
                  'Ki67':'Ki67',
                  'FUS':'FUS',
                  'Mg':'Magnesium',
                  'Ca':'Calcium'}

cell_fullname = {'H1':'H1-hESC',
                 'GM':'GM12878',
                 'mCD8T':'Mouse CD8+ T cell'}

agent_charge = {'sp':'4+',
                'spd':'3+',
                'CoH':'3+',
                'PEG':'',
                'Ca':'2+',
                'Mg':'2+',
                'HP1a':'',
                'HP1bSUV':'',
                'LKH':'',
                'Ki67':'',
                'FUS':''}


# experiment list (cell, sample, agent, tnum, rep)
#exp_list = [('H1', 'NCP', 'sp'),
#            ('H1', 'NCP', 'spd'),
#            ('H1', 'NCP', 'CoH'),
#            ('H1', 'NCP', 'PEG'),
#            ('H1', 'NCP', 'Ca'),
#            ('H1', 'NCP', 'HP1a'),
#            ('H1', 'NCP', 'LKH'),
#            ('H1', 'NCP', 'Ki67'),
#            ('H1', 'DNA', 'HP1a'),
#            ('H1', 'DNA', 'LKH'),
#            ('H1', 'DNA', 'Ki67'),
#            ('mCD8T', 'WT-NCP', 'sp'),
#            ('mCD8T', 'inht-NCP', 'sp'),
#            ('mCD8T', 'KO-NCP', 'sp')]

#exp_list = [('H1', 'NCP', 'sp', 8, 1),
#            ('H1', 'NCP', 'sp', 8, 2),
#            ('H1', 'NCP', 'spd', 6, 1),
#            ('H1', 'NCP', 'spd', 6, 2),
#            ('H1', 'NCP', 'CoH', 5, 1),
#            ('H1', 'NCP', 'CoH', 5, 2),
#            ('H1', 'NCP', 'PEG', 6, 1),
#            ('H1', 'NCP', 'PEG', 6, 2),
#            ('H1', 'NCP', 'Ca', 5, 1),
#            ('H1', 'NCP', 'Ca', 5, 2),
#            ('H1', 'NCP', 'HP1a', 3, 1),
#            ('H1', 'NCP', 'HP1a', 3, 2),
#            ('H1', 'NCP', 'LKH', 3, 1),
#            ('H1', 'NCP', 'LKH', 3, 2),
#            ('H1', 'NCP', 'Ki67', 4, 1),
#            ('H1', 'NCP', 'Ki67', 4, 2)]

#exp_list = [('H1', 'NCP', 'sp', 8, 1),
#            ('H1', 'NCP', 'sp', 8, 2),
#            ('H1', 'NCP', 'spd', 6, 1),
#            ('H1', 'NCP', 'spd', 6, 2),
#            ('H1', 'NCP', 'CoH', 5, 1),
#            ('H1', 'NCP', 'CoH', 5, 2),
#            ('H1', 'NCP', 'PEG', 6, 1),
#            ('H1', 'NCP', 'PEG', 6, 2),
#            ('H1', 'NCP', 'Ca', 5, 1),
#            ('H1', 'NCP', 'Ca', 5, 2)]

#exp_list = [('H1', 'NCP', 'HP1a', 3, 1),
#            ('H1', 'NCP', 'HP1a', 3, 2),
#            ('H1', 'NCP', 'LKH', 4, 1),
#            ('H1', 'NCP', 'LKH', 4, 2),
#            ('H1', 'NCP', 'Ki67', 4, 1),
#            ('H1', 'NCP', 'Ki67', 4, 2)]

#exp_list = [('H1', 'DNA', 'HP1a', 3, 1),
#            ('H1', 'DNA', 'HP1a', 3, 2),
#            ('H1', 'DNA', 'LKH', 4, 1),
#            ('H1', 'DNA', 'LKH', 4, 2),
#            ('H1', 'DNA', 'Ki67', 4, 1),
#            ('H1', 'DNA', 'Ki67', 4, 2)]

#exp_list = [('GM', 'NCP', 'sp', 8, 1),
#            ('GM', 'NCP', 'sp', 8, 2)]

#exp_list = [('mCD8T', 'WT-NCP', 'sp', 8, 1),
#            ('mCD8T', 'WT-NCP', 'sp', 8, 2),
#            ('mCD8T', 'inht-NCP', 'sp', 8, 1),
#            ('mCD8T', 'inht-NCP', 'sp', 8, 2),
#            ('mCD8T', 'KO-NCP', 'sp', 8, 1),
#            ('mCD8T', 'KO-NCP', 'sp', 8, 2)]



exp_list = [('H1', 'NCP', 'sp', 8, 1),
            ('H1', 'NCP', 'sp', 8, 2),
            ('H1', 'NCP', 'spd', 6, 1),
            ('H1', 'NCP', 'spd', 6, 2),
            ('H1', 'NCP', 'CoH', 5, 1),
            ('H1', 'NCP', 'CoH', 5, 2),
            ('H1', 'NCP', 'PEG', 6, 1),
            ('H1', 'NCP', 'PEG', 6, 2),
            ('H1', 'NCP', 'Ca', 5, 1),
            ('H1', 'NCP', 'Ca', 5, 2),
            ('H1', 'NCP', 'HP1a', 3, 1),
            ('H1', 'NCP', 'HP1a', 3, 2),
            ('H1', 'NCP', 'LKH', 4, 1),
            ('H1', 'NCP', 'LKH', 4, 2),
            ('H1', 'NCP', 'Ki67', 4, 1),
            ('H1', 'NCP', 'Ki67', 4, 2),
            ('H1', 'DNA', 'HP1a', 3, 1),
            ('H1', 'DNA', 'HP1a', 3, 2),
            ('H1', 'DNA', 'LKH', 4, 1),
            ('H1', 'DNA', 'LKH', 4, 2),
            ('H1', 'DNA', 'Ki67', 4, 1),
            ('H1', 'DNA', 'Ki67', 4, 2),
            ('GM', 'NCP', 'sp', 8, 1),
            ('GM', 'NCP', 'sp', 8, 2),
            ('mCD8T', 'WT-NCP', 'sp', 8, 1),
            ('mCD8T', 'WT-NCP', 'sp', 8, 2),
            ('mCD8T', 'inht-NCP', 'sp', 8, 1),
            ('mCD8T', 'inht-NCP', 'sp', 8, 2),
            ('mCD8T', 'KO-NCP', 'sp', 8, 1),
            ('mCD8T', 'KO-NCP', 'sp', 8, 2)]


exp_list = [('H1', 'NCP', 'sp', 4, 1),
            ('H1', 'NCP', 'sp', 4, 2),
            ('H1', 'NCP', 'sp', 8, 1),
            ('H1', 'NCP', 'sp', 8, 2)]



exp_list = [('GM', 'NCP', 'sp', 4, 1),
            ('GM', 'NCP', 'sp', 4, 2),
            ('GM', 'NCP', 'sp', 8, 1),
            ('GM', 'NCP', 'sp', 8, 2)]


exp_list = [('mCD8T', 'WT-NCP', 'sp', 4, 1),
            ('mCD8T', 'WT-NCP', 'sp', 4, 2),
            ('mCD8T', 'WT-NCP', 'sp', 8, 1),
            ('mCD8T', 'WT-NCP', 'sp', 8, 2),
            ('mCD8T', 'inht-NCP', 'sp', 4, 1),
            ('mCD8T', 'inht-NCP', 'sp', 4, 2),
            ('mCD8T', 'inht-NCP', 'sp', 8, 1),
            ('mCD8T', 'inht-NCP', 'sp', 8, 2),
            ('mCD8T', 'KO-NCP', 'sp', 4, 1),
            ('mCD8T', 'KO-NCP', 'sp', 4, 2),
            ('mCD8T', 'KO-NCP', 'sp', 8, 1),
            ('mCD8T', 'KO-NCP', 'sp', 8, 2)]

exp_list = [('H1', 'NCP', 'sp', 0, 1),
            ('H1', 'NCP', 'sp', 0, 2),
            ('H1', 'NCP', 'spd', 0, 1),
            ('H1', 'NCP', 'spd', 0, 2),
            ('H1', 'NCP', 'CoH', 0, 1),
            ('H1', 'NCP', 'CoH', 0, 2),
            ('H1', 'NCP', 'PEG', 0, 1),
            ('H1', 'NCP', 'PEG', 0, 2),
            ('H1', 'NCP', 'Ca', 0, 1),
            ('H1', 'NCP', 'Ca', 0, 2),
            ('H1', 'NCP', 'HP1a', 0, 1),
            ('H1', 'NCP', 'HP1a', 0, 2),
            ('H1', 'NCP', 'LKH', 0, 1),
            ('H1', 'NCP', 'LKH', 0, 2),
            ('H1', 'NCP', 'Ki67', 0, 1),
            ('H1', 'NCP', 'Ki67', 0, 2),
            ('H1', 'DNA', 'HP1a', 0, 1),
            ('H1', 'DNA', 'HP1a', 0, 2),
            ('H1', 'DNA', 'LKH', 0, 1),
            ('H1', 'DNA', 'LKH', 0, 2),
            ('H1', 'DNA', 'Ki67', 0, 1),
            ('H1', 'DNA', 'Ki67', 0, 2),
            ('GM', 'NCP', 'sp', 0, 1),
            ('GM', 'NCP', 'sp', 0, 2),
            ('mCD8T', 'WT-NCP', 'sp', 0, 1),
            ('mCD8T', 'WT-NCP', 'sp', 0, 2),
            ('mCD8T', 'inht-NCP', 'sp', 0, 1),
            ('mCD8T', 'inht-NCP', 'sp', 0, 2),
            ('mCD8T', 'KO-NCP', 'sp', 0, 1),
            ('mCD8T', 'KO-NCP', 'sp', 0, 2)]


#exp_list = [('H1', 'NCP', 'sp', i, 1) for i in range(0, 10)]
#exp_list += [('H1', 'NCP', 'sp', i, 2) for i in range(0, 10)]

#exp_list = [('GM', 'NCP', 'sp', i, 1) for i in range(1, 10)]
#exp_list += [('GM', 'NCP', 'sp', i, 2) for i in range(1, 10)]

#exp_list = [('mCD8T', 'WT-NCP', 'sp', i, 1) for i in range(1, 10)]
#exp_list += [('mCD8T', 'WT-NCP', 'sp', i, 2) for i in range(1, 10)]

#exp_list = [('mCD8T', 'inht-NCP', 'sp', i, 1) for i in range(1, 10)]
#exp_list += [('mCD8T', 'inht-NCP', 'sp', i, 2) for i in range(1, 10)]

exp_list = [('mCD8T', 'KO-NCP', 'sp', i, 1) for i in range(1, 10)]
exp_list += [('mCD8T', 'KO-NCP', 'sp', i, 2) for i in range(1, 10)]

#exp_list = [('H1', 'NCP', 'HP1a', i, 1) for i in range(1, 6)]
#exp_list += [('H1', 'NCP', 'HP1a', i, 2) for i in range(1, 6)]

#exp_list = [('H1', 'NCP', 'LKH', i, 1) for i in range(1, 6)]
#exp_list += [('H1', 'NCP', 'LKH', i, 2) for i in range(1, 6)]

#exp_list = [('H1', 'NCP', 'Ki67', i, 1) for i in range(1, 6)]
#exp_list += [('H1', 'NCP', 'Ki67', i, 2) for i in range(1, 6)]












# bin size
bin_size = 10000
#bin_size = 5000
#bin_size = 1000

# chromosome choice
chr_choices = ['chr1']

# data type
#dtype = 'num'
dtype = 'zscore'
#dtype = 'zChalf'
#dtype = 'Chalf'

# read data
exp_ID_pos = {}
exp_ID_score = {}
for exp in exp_list:
    cell, sample, agent, tnum, rep = exp
    path = rep_path[rep]

    fname = '_'.join([cell, sample, agent,
                      str(int(bin_size/1000.0)) + 'kb',
                      dtype]) + '.cn'

    if dtype in ['score', 'zscore', 'num']:
        field_name = '-'.join([cell, sample, agent, str(tnum)])
        
    elif dtype in ['Chalf', 'zChalf']:
        field_name = 'Chalf'

    if dtype in ['num']:
        ID_pos, field_ID_value = read_bin_num(path + fname, chr_choices=chr_choices) 
    else:
        ID_pos, field_ID_value = read_bin_file(path + fname, chr_choices=chr_choices) 

    ID_score = field_ID_value[field_name]
        
    exp_ID_pos[exp] = ID_pos
    exp_ID_score[exp] = ID_score

# pairing by replicates
rep_pairs = []
for i in range(len(exp_list)-1):
    for j in range(i+1, len(exp_list)):
        exp1, exp2 = exp_list[i], exp_list[j]
        if exp1[:-1] == exp2[:-1]: # same (cell, sample, agent, tnum)
        #if exp1[:-2] == exp2[:-2]: # same (cell, sample, agent)
            rep_pairs.append((exp1, exp2))


# get correlation between replicates
pair_corr = {}
for rep_pair in rep_pairs:
    exp1, exp2 = rep_pair
    
    ID_list = set(exp_ID_pos[exp1].keys()) & set(exp_ID_pos[exp2].keys())
    ID_list = list(ID_list)
        
    cell, sample, agent, tnum, rep = exp1

    X, Y = [], []
    for ID in ID_list:
        X.append(exp_ID_score[exp1][ID])
        Y.append(exp_ID_score[exp2][ID])

    #corr = scipy.stats.spearmanr(X, Y)[0]
    corr = scipy.stats.pearsonr(X, Y)[0]
    #print ("%s VS %s: %1.2f" % (agent1, agent2, corr))

    pair_corr[(exp1, exp2)] = corr

    """
    fig = plt.figure()
    plt.plot(X, Y, 'k.', alpha=0.2)
    #plt.annotate("Spearman %1.2f" % (corr), xy=(0.2, 0.75), fontsize=12, xycoords='axes fraction')
    plt.annotate("Pearson %1.2f" % (corr), xy=(0.2, 0.75), fontsize=12, xycoords='axes fraction')
    #plt.title("%s VS %s" % (agent, agent2))
    #plt.xlabel("fold change"  + str(exp1))
    #plt.ylabel("fold change"  + str(exp2))
    plt.xlabel("replicate #1")
    plt.ylabel("replicate #2")
    plt.xscale('log', basex=2)
    plt.yscale('log', basey=2)
    plt.savefig("%s-%s-%s-%d_rep_corr.png" % (cell, sample ,agent, tnum), dpi=100, bbox_inches='tight')
    #plt.show()
    plt.close()
    """

    fig = plt.figure()
    density_scatter(X, Y, bins = [20,20], s=3, cmap=pastel_jet, ax=plt.gca())
    plt.plot([-10, 10], [-10, 10], 'k--', alpha=0.5)
    plt.annotate("Pearson %1.2f" % (corr), xy=(0.2, 0.75), fontsize=12, xycoords='axes fraction')
    plt.xlim([-5, 5])
    plt.ylim([-5, 5])
    plt.xlabel("replicate #1")
    plt.ylabel("replicate #2")                 
    plt.title('%s %s %s t#%d' % (cell_fullname[cell], sample, agent_fullname[agent], tnum))
    #plt.savefig("%s-%s-%s_rep_corr.svg" % (cell, sample ,agent), format='svg', bbox_inches='tight')
    plt.savefig("%s-%s-%s-%d_rep_corr.png" % (cell, sample ,agent, tnum), dpi=100, bbox_inches='tight')
    plt.close()
