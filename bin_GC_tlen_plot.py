import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import copy
import math
import scipy
from scipy.stats import gaussian_kde
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn
from matplotlib.colors import LinearSegmentedColormap

agent_fullname = {'sp':'Spermine(4+)',
                  'spd':'Spermidine(3+)',
                  'CoH':'Cobalt Hexammine(3+)',
                  'PEG':'PEG 8000',
                  'HP1a':'HP1 $\\alpha$',
                  'LKH':'Linker histone H1',
                  'Ki67':'Ki67',
                  'Mg':'Magnesium',
                  'Ca':'Calcium'}

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

def total (chr_binID_count):
    total = 0
    for chr in chr_binID_count:
        total += sum(chr_binID_count[chr])
    return total

def GC_content(seq):
    seq = seq.upper()
    output = 0.0
    for nt in seq:
        if nt in "GC":
            output += 1.0
        elif nt in 'N':
            output += 0.5
    return output/len(seq)

# read bin score file
def read_bin_score (fname, bin_size):
    name_chr_binID_score = {}
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
        #binID = st / int(bin_size)
        binID = (st, ed)
        for name, score in zip(names, scores):
            if name not in name_chr_binID_score:
                name_chr_binID_score[name] = {}
            if chr not in name_chr_binID_score[name]:
                name_chr_binID_score[name][chr] = {}
            name_chr_binID_score[name][chr][binID] = score
    return name_chr_binID_score



#### parameters
path = "/home/spark159/../../storage/"

# experiment list (cell, sample, agent, tnum)

#exp_list = [(cell, 'NCP', 'sp', 8),
#            (cell, 'NCP', 'spd', 6),
#            (cell, 'NCP', 'CoH', 5),
#            (cell, 'NCP', 'PEG', 6),
#            (cell, 'NCP', 'Ca', 5),
#            (cell, 'NCP', 'Mg', 5),
#            (cell, 'NCP', 'HP1a', 3),
#            (cell, 'NCP', 'HP1bSUV', 4),
#            (cell, 'NCP', 'LKH', 3),
#            (cell, 'NCP', 'Ki67', 4),
#            (cell, 'NCP', 'FUS', 5)]

exp_list = [('H1', 'NCP', 'sp', 8),
            ('H1', 'NCP', 'HP1a', 3),
            ('H1', 'NCP', 'LKH', 3),
            ('H1', 'NCP', 'Ki67', 4)]

#exp_list = [('H1', 'NCP', 'sp', 8),
#            ('H1', 'NCP', 'HP1a', 3)]



# binsize of input data
#bin_size = 10000
#bin_size = 5000
bin_size = 1000


# other parameters
dtype = 'score'
note = ''

if dtype == 'score':
    ylabel = 'Condensability (A.U.)'
elif dtype == 'zscore':
    ylabel = 'Z-score'


### get data and correlation
exp_datas = {}
exp_pcorr = {}
exp_lims = {}
for exp in exp_list:

    cell, sample, agent, tnum = exp

    print cell, sample, agent, tnum
    
    # load score bin file
    fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', dtype]) + '.cn'

    if dtype in ['Chalf', 'zChalf']:
        field_name = 'Chalf'
    else:
        field_name = "%s-%s-%s-%d.bam" % (cell, sample, agent, tnum)
    
    chr_binID_score = read_bin_score(fname, bin_size)[field_name]

    # load bin file
    fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', 'bin']) + '.cn'
    name_chr_binID_score = read_bin_score(fname, bin_size)
    chr_binID_GC = name_chr_binID_score['GCcontent']
    chr_binID_tlen = name_chr_binID_score['Meantlen']
    del name_chr_binID_score

    # get data
    name_data = {}
    for chr in chr_binID_score:
        binID_score = chr_binID_score[chr]
        binID_GC = chr_binID_GC[chr]
        binID_tlen = chr_binID_tlen[chr]
        binIDs = list(set(binID_score) & set(binID_GC) & set(binID_tlen))
        for binID in binIDs:
            score = binID_score[binID]
            AT = (1.0 - binID_GC[binID])*100
            tlen = binID_tlen[binID]
            if 'score' not in name_data:
                name_data['score'] = []
            if 'AT' not in name_data:
                name_data['AT'] = []
            if 'tlen' not in name_data:
                name_data['tlen'] = []
            name_data['score'].append(score)
            name_data['AT'].append(AT)
            name_data['tlen'].append(tlen)

    exp_datas[exp] = name_data

    names = ['AT', 'tlen', 'score']

    # check correlation
    pair_corr = {}
    for i in range(len(names)-1):
        for j in range(i+1, len(names)):
            name1, name2 = names[i], names[j]
            corr = scipy.stats.spearmanr(name_data[name1], name_data[name2])[0]
            pair_corr[(name1, name2)] = corr
            pair_corr[(name2, name1)] = corr
            print '%s VS %s' % (name1, name2), corr

    exp_pcorr[exp] = pair_corr


    # set ranges for plots
    name_lims = {}
    for name in names:
        data = name_data[name]
        median = np.median(data)
        std = np.std(data)
        name_lims[name] = [median-5*std, median+5*std]

    exp_lims[exp] = name_lims

    print


### draw density scatter plot
if False:    
    for exp in exp_list:
        cell, sample, agent, tnum = exp

        name_data = exp_datas[exp]
        pair_corr = exp_pcorr[exp]
        name_lims = exp_lims[exp]

        # setgraph parameters
        name_label = {'AT':'AT content (%)',
                  'score':'Condensabiltiy (A.U.)',
                  'tlen':'Mean DNA length (bp)'}

        # plot density scatter plot
        for i in range(len(names)-1):
            for j in range(i+1, len(names)):
                name1, name2 = names[i], names[j]
                data1, data2 = name_data[name1], name_data[name2]

                hist_bins = []
                for name in [name1, name2]:
                    if name == 'tlen':
                        hist_bins.append(300)
                    else:
                        hist_bins.append(30)

                fig = plt.figure()
                density_scatter(data1,
                                data2,
                                bins = hist_bins,
                                s=3,
                                cmap=pastel_jet,
                                ax=plt.gca(),
                                sort=False)
                plt.xlabel(name_label[name1])
                plt.ylabel(name_label[name2])
                plt.xlim(name_lims[name1])
                plt.ylim(name_lims[name2])
                plt.savefig("%sVS%s_%s_%s_%s_%s.png" % (name1, name2, cell, sample, agent, dtype),
                            bbox_inches='tight')
                plt.close()


### draw correlation matrix
if False:
    for exp in exp_list:
        cell, sample, agent, tnum = exp

        name_data = exp_datas[exp]
        pair_corr = exp_pcorr[exp]
        name_lims = exp_lims[exp]

        # setgraph parameters
        name_label = {'AT':'AT content',
                  'score':'Condensabiltiy',
                  'tlen':'Mean\nDNA length'}
        
        fig, axes = plt.subplots(figsize=(4, 3),
                                 nrows=len(names),
                                 ncols=len(names))

        for i in range(len(names)):
            for j in range(len(names)):
                name1, name2 = names[i], names[j]
                if i > j:
                    data1, data2 = name_data[name1], name_data[name2]
                    axes[i,j].plot(data1, data2, 'k.', ms=0.5, mfc='k', mec='k', alpha=0.05)
                    axes[i,j].set_xlim(name_lims[name1])
                    axes[i,j].set_ylim(name_lims[name2])

                    if j > 0 and i < len(names) -1:
                        axes[i,j].tick_params(axis='both',
                                              which='both',
                                              labelbottom=False,
                                              labelleft=False)
                    if j == 0 and i < len(names) -1:
                        axes[i,j].tick_params(axis='x',
                                              which='both',
                                              labelbottom=False)
                    if j > 0 and i == len(names) - 1:
                        axes[i,j].tick_params(axis='y',
                                              which='both',
                                              labelleft=False)

                elif i == j:
                    matrix = np.zeros((len(names), len(names)))
                    matrix[:] = np.nan
                    axes[i,j].imshow(matrix, origin='lower')
                    label = name_label[name1]
                    axes[i,j].text(len(names)/2,
                                   len(names)/2,
                                   label,
                                   ha="center",
                                   va="center",
                                   fontsize=6,
                                   weight='bold')
                    axes[i,j].set_xlim([0, len(names)-1])
                    axes[i,j].set_ylim([0, len(names)-1])
                    axes[i,j].set_axis_off()
                else:
                    assert i < j
                    corr = pair_corr[(name1, name2)]
                    matrix = np.zeros((len(names), len(names)))
                    matrix[:] = corr
                    img = axes[i,j].imshow(matrix, cmap="bwr", vmin=-0.8, vmax=0.8, origin='lower')
                    if abs(corr) > 0.6:
                        color = "white"
                    else:
                        color = "black"
                    axes[i,j].text(len(names)/2,
                                   len(names)/2,
                                   str(round(corr,2)),
                                   ha="center",
                                   va="center",
                                   fontsize=10,
                                   color=color,
                                   weight='bold')
                    axes[i,j].set_xlim([0, len(names)-1])
                    axes[i,j].set_ylim([0, len(names)-1])
                    axes[i,j].tick_params(axis='both',
                                          which='both',
                                          bottom=False,
                                          top=False,
                                          left=False,
                                          right=False,
                                          labelbottom=False,
                                          labeltop=False,
                                          labelleft=False,
                                          labelright=False)

        plt.subplots_adjust(wspace=0.1, hspace=0.1)
        cbar=fig.colorbar(img, ax=axes, location='right', shrink=0.8)
        cbar.ax.set_ylabel('Spearman correlation', rotation=-90, va="bottom")
        #plt.suptitle("Correlation betwen condensing agents (10 kb bin)")
        plt.savefig("Corr_GC_tlen_%s_%s_%s_%s.png" % (cell, sample, agent, dtype),
                    dpi=1000,
                    bbox_inches='tight')
        #plt.show()
        plt.close()

# plot bargraph of correations
if True:
    # setgraph parameters
    name_label = {'AT':'AT content',
                  'score':'Condensabiltiy',
                  'tlen':'Mean DNA length'}

    fig, axes = plt.subplots(nrows=1, ncols=len(exp_list), sharey=True)

    for i in range(len(exp_list)):
        exp = exp_list[i]
        cell, sample, agent, tnum = exp
        
        pair_corr = exp_pcorr[exp]

        Y = [pair_corr[('AT','score')], pair_corr[('tlen','score')]]
        xlabels = [name_label['AT'], name_label['tlen']]

        barlist = axes[i].bar(range(len(Y)), Y, width=0.5)
        axes[i].set_xticks(range(len(Y)))
        axes[i].set_xticklabels(xlabels, rotation=45, ha='right', rotation_mode='anchor')
        axes[i].set_title(agent_fullname[agent])

        if i > 0:
            axes[i].spines['left'].set_visible(False)
            axes[i].tick_params(left=False)
            
        axes[i].spines['top'].set_visible(False)
        axes[i].spines['right'].set_visible(False)
        axes[i].axhline(y=0, linestyle='-', color='k', alpha=1.0, lw=0.7)

        for bar, corr in zip(barlist, Y):
            if corr >=0:
                color = 'tab:blue'
            else:
                color = 'tab:red'
            bar.set_color(color)

    axes[0].set_ylabel('Spearman correlation')
    plt.savefig("bar.svg", format='svg', bbox_inches='tight')
    plt.close()
