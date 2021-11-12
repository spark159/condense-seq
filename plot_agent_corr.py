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
import random
import matplotlib.cm as cm
import scipy

# rescale the data in old range (old_st, old_ed) into new range (new_st, new_ed)
def rescale (value_list, old_st, old_ed, new_st, new_ed):
    output = []
    for value in value_list:
        if value >= old_st and value <= old_ed:
            new_value = new_st + (new_ed - new_st)*float(value-old_st)/(old_ed-old_st)
        else:
            new_value = np.nan
        output.append(new_value)
    return output

def read_Gband (fname, chr_choices=None):
    chr_ID_Gband = {}
    ID = 0
    for line in open(fname):
        if line.startswith("#"):
            continue
        cols = line.strip().split()
        chr, start, end, name, stain = cols
        if chr_choices !=None and chr not in chr_choices:
            continue
        start, end = int(start), int(end)
        if stain.startswith('g'):
            type = stain[1:4]
            if type == 'neg':
                value = 0
            elif type == 'pos':
                value = int(stain[4:])
            else:
                assert type == 'var'
                value = np.nan
        else:
            type = stain
            value = np.nan
        if chr not in chr_ID_Gband:
            chr_ID_Gband[chr] = {}
        assert ID not in chr_ID_Gband[chr]
        chr_ID_Gband[chr][ID] = {}
        chr_ID_Gband[chr][ID]['interval'] = (start, end)
        chr_ID_Gband[chr][ID]['name'] = name
        chr_ID_Gband[chr][ID]['type'] = type
        chr_ID_Gband[chr][ID]['value'] = value
        ID +=1
    return chr_ID_Gband


def read_bin (fname, chr_choices=None):
    First = True
    ID_pos = {}
    name_ID_score = {}
    for line in open(fname):
        cols = line.strip().split()
        if First:
            names = cols[4:-1]
            First = False
            continue
        ID, chr, start, end = cols[:4]
        start, end = int(start), int(end)
        if chr_choices!=None and chr not in chr_choices:
            continue
        pos = int(float(start + end)/2)
        assert ID not in ID_pos
        ID_pos[ID] = pos
        control_count = float(cols[-2])
        for i in range(len(names)-1):
            name = names[i]
            #control_count = float(cols[-2])
            #if control_count <= 0:
            #    rcount = np.nan
            #else:
            #    #rcount = -np.log(float(cols[4+i]) / control_count)
            #    #rcount = float(cols[4+i]) / control_count
            #    rcount = np.log(control_count)
            #rcount = -np.log((float(cols[4+i])+1) / (control_count+1))
            #rcount = np.log(float(cols[4+i])+1)
            #rcount = float(cols[-1])
            score = -np.log((float(cols[4+i])+1) / (control_count+1))
            if name not in name_ID_score:
                name_ID_score[name] = {}
            assert ID not in name_ID_score[name]
            #name_ID_score[name][ID] = rcount
            name_ID_score[name][ID] = score

    return ID_pos, name_ID_score

# read files
chr_choices = ['chr1']
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'Ca', 'Mg']
agent_titrnum = {'sp':8, 'spd':6, 'CoH':5, 'PEG':6, 'Ca':5, 'Mg':5}
agent_ID_score = {}
for agent in agent_list:
    fname = 'H1_NCP_%s_10kb_bin.cn' % (agent)
    field_name = 'H1-NCP-%s-%d.bam' % (agent, agent_titrnum[agent])
    ID_pos, field_ID_value = read_bin(fname, chr_choices=chr_choices)
    agent_ID_score[agent] = field_ID_value[field_name]


# compare between condensing agents
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'Ca', 'Mg']
agent_charge = {'sp':'4+', 'spd':'3+', 'CoH':'3+', 'PEG':'', 'Ca':'2+', 'Mg':'2+'}
ID_list = ID_pos.keys()
pair_corr = {}
corr_matrix = [ [np.nan for i in range(len(agent_list))] for i in range(len(agent_list)) ]
for i in range(len(agent_list)-1):
    for j in range(i+1, len(agent_list)):
        agent1 = agent_list[i]
        agent2 = agent_list[j]

        X, Y = [], []
        for ID in ID_list:
            X.append(agent_ID_score[agent1][ID])
            Y.append(agent_ID_score[agent2][ID])

        corr = scipy.stats.spearmanr(X, Y)[0]
        #corr = scipy.stats.pearsonr(X, Y)[0]
        #print ("%s VS %s: %1.2f" % (agent1, agent2, corr))
            
        fig = plt.figure()
        plt.plot(X, Y, '.')
        plt.annotate("Spearman %1.2f" % (corr), xy=(0.2, 0.75), fontsize=12, xycoords='axes fraction')
        plt.title("%s VS %s" % (agent1, agent2))
        plt.xlabel("fold change (%s)" % (agent1))
        plt.ylabel("fold change (%s)" % (agent2))
        #plt.xscale('log', base=2)
        #plt.yscale('log', base=2)
        #plt.show()
        plt.close()

        pair_corr[(agent1, agent2)] = corr
        corr_matrix[i][j] = corr
        corr_matrix[j][i] = corr


fig, axes = plt.subplots(nrows=len(agent_list), ncols=len(agent_list))
for i in range(len(agent_list)):
    for j in range(len(agent_list)):
        idx = len(agent_list)*i + j
        agent1, agent2 = agent_list[i], agent_list[j]
        if i > j:
            X = [ agent_ID_score[agent1][ID] for ID in ID_list ]
            Y = [ agent_ID_score[agent2][ID] for ID in ID_list ]
            axes[i,j].plot(X, Y, 'k,', alpha=0.1)
            #axes[i,j].set_xscale('log', base=2)
            #axes[i,j].set_yscale('log', base=2)
            if j > 0 and i < len(agent_list) -1:
                axes[i,j].tick_params(axis='both', which='both', labelbottom=False, labelleft=False)
            if j == 0 and i < len(agent_list) -1:
                axes[i,j].tick_params(axis='x', which='both', labelbottom=False)
            if j > 0 and i == len(agent_list) - 1:
                axes[i,j].tick_params(axis='y', which='both', labelleft=False)
                
        elif i == j:
            matrix = np.zeros((len(agent_list), len(agent_list)))
            matrix[:] = np.nan
            axes[i,j].imshow(matrix, origin='lower')
            s = agent1 + '$^{%s}$' % (agent_charge[agent1])
            axes[i,j].text(len(agent_list)/2, len(agent_list)/2, s, ha="center", va="center", fontsize=10, weight='bold')
            axes[i,j].set_xlim([0, len(agent_list)-1])
            axes[i,j].set_ylim([0, len(agent_list)-1])
            axes[i,j].set_axis_off()
        else:
            assert i < j
            value = pair_corr[(agent1, agent2)]
            matrix = np.zeros((len(agent_list), len(agent_list)))
            matrix[:] = value
            img = axes[i,j].imshow(matrix, cmap="jet", vmin=0.1, vmax=0.8, origin='lower')
            if value < 0.3 or value > 0.7:
                color = "white"
            else:
                color = "black"
            axes[i,j].text(len(agent_list)/2, len(agent_list)/2, str(round(value,2)), ha="center", va="center", fontsize=10, color=color, weight='bold')
            axes[i,j].set_xlim([0, len(agent_list)-1])
            axes[i,j].set_ylim([0, len(agent_list)-1])
            axes[i,j].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labeltop=False, labelleft=False, labelright=False)

plt.subplots_adjust(wspace=0.1, hspace=0.1)
cbar=fig.colorbar(img, ax=axes, location='right', shrink=0.8)
cbar.ax.set_ylabel('Spearman correlation', rotation=-90, va="bottom")
plt.suptitle("Corrrelation betwen condensing agents (10 kb bin)")
#plt.show()
plt.close()

# make genome-wide plot
# set binning resolution
#i = 20
#bin_size = int(0.5*(10**6) / i) # binsize (unit of bp)
#blur_win = int(4*i + 1) # sliding window (unit of bin)

bin_size = 10**6

# set target names and feature names
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'Ca', 'Mg']
agent_minmax = {'sp':(-0.4, 0.6), 'spd':(-0.7, 0.5), 'CoH':(-0.5, 2), 'PEG':(-0.3, 0.5), 'Ca':(-0.15, 0.25), 'Mg':(0.1, 0.4)}
agent_color = {'sp':'tab:blue', 'spd':'tab:orange', 'CoH':'tab:green', 'PEG':'tab:red', 'Ca':'tab:purple', 'Mg':'tab:brown'}
target_names = agent_list
feature_names = ['Gene activity']
names = agent_list + feature_names

# chromosome cartoon aspect ratio
aspect = (0.05*bin_size) / (10**6)

#### Binning the data and draw the genome-wide plot by chromosome
for chr_choice in chr_choices:

    #### read files
    path = ""
    # read G-band information
    chr_gID_Gband = read_Gband(path+"Gband_information.txt", chr_choices=[chr_choice])

    # read annotation file
    #field_ID_value = load_file.read_tabular_file (path+"H1_NCP_sp_10kb_anot.cn", mode='col')
    #ID_pos = field_ID_value['PhysicalPosition']
    #ID_pos, field_ID_value = read_bin("H1_NCP_sp_10kb_bin.cn", chr_choices=[chr_choice])

    # read GTF file
    geneID_field_values, field_geneID_values = load_file.read_GTF (path+"ENCFF159KBI.gtf", mode="both", chr_list=[chr_choice])
    geneID_pos = {}
    for geneID in geneID_field_values:
        try:
            pos = geneID_field_values[geneID]['TSS']
            geneID_pos[geneID] = pos
        except:
            continue

    # read tsv file
    geneID_FPKM = load_file.read_tsv(path+"ENCFF174OMR.tsv")

    # read compartment score file
    #eigenbinID_value, eigenbinID_interval = read_eigenfile(path+"eigen_H1_100kb.txt", bin_size=100000)


    #### binning the data
    name_binID_mean = {}
    name_binID_count = {}

    # binning the condensability data
    for ID in ID_pos:
        binID = int(ID_pos[ID]) / int(bin_size)
        for name in agent_list:
            if name not in name_binID_mean:
                name_binID_mean[name] = {}
            if name not in name_binID_count:
                name_binID_count[name] = {}
            try:
                value = agent_ID_score[name][ID]
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

    # binning the RNA-seq data
    temp = {}
    for ID in ID_pos:
        binID = int(ID_pos[ID]) / int(bin_size)
        temp[binID] = 0

    name_binID_mean['Gene density'] = copy.deepcopy(temp)
    name_binID_count['Gene density'] = copy.deepcopy(temp)
    name_binID_mean['Gene activity'] = copy.deepcopy(temp)
    name_binID_count['Gene activity'] = copy.deepcopy(temp)

    min_FPKM = min(geneID_FPKM.values())
    for geneID in geneID_pos:
        binID = int(geneID_pos[geneID]) / int(bin_size)
        try:
            name_binID_mean['Gene density'][binID] += 1.0
        except:
            continue
        try:
            #name_binID_mean['Gene activity'][binID] += geneID_FPKM[geneID]
            #name_binID_mean['Gene activity'][binID] += np.log2(geneID_FPKM[geneID])
            name_binID_mean['Gene activity'][binID] += np.log2(geneID_FPKM[geneID] - min_FPKM + 1)
        except:
            #name_binID_mean['Gene density'][binID] += 0.0
            name_binID_mean['Gene activity'][binID] += np.nan
        name_binID_count['Gene activity'][binID] += 1

    for binID in name_binID_mean['Gene activity']:
        if name_binID_count['Gene activity'][binID] <= 0:
            continue
        name_binID_mean['Gene activity'][binID] /= name_binID_count['Gene activity'][binID] 

    #names.append('Gene density')
    names.append('Gene activity')

    ## binning the compartment score data
    #binID_eigens = {}
    #for value, interval in zip(eigenbinID_value, eigenbinID_interval):
    #    st, ed = interval
    #    st_binID, ed_binID = st / bin_size, ed / bin_size
    #    for binID in range(st_binID, ed_binID):
    #        if binID not in binID_eigens:
    #            binID_eigens[binID] = []
    #        binID_eigens[binID].append(value)

    #binID_eigen = {}
    #for binID in binID_eigens:
    #    try:
    #        binID_eigen[binID] = np.mean(binID_eigens[binID])
    #    except:
    #        binID_eigen[binID] = np.nan
    #name_binID_mean['eigen'] = binID_eigen


    # binning the G-banding data and make ideogram
    gID_Gband = chr_gID_Gband[chr_choice]
    binID_Gvalue, binID_Gtype = {}, {}
    gID_binterval = {}
    gband_img, gband_cenimg, gband_varimg = [], [], []
    gtick_locs, gtick_labels = [], []
    for gID in sorted(gID_Gband.keys()):
        st, ed = gID_Gband[gID]['interval']
        st_binID, ed_binID = st / bin_size, ed / bin_size
        gtype = gID_Gband[gID]['type']
        gvalue = gID_Gband[gID]['value']
        for binID in range(st_binID, ed_binID):
            binID_Gvalue[binID] = gID_Gband[gID]['value']
            binID_Gtype[binID] = gID_Gband[gID]['type']
        gID_binterval[gID] = (st_binID, ed_binID)
        bsize = ed_binID - st_binID
        gband_img += [[gvalue] for k in range(bsize)] 
        if gtype == 'acen':
            gband_cenimg += [[10] for k in range(bsize)]
        else:
            gband_cenimg += [[np.nan] for k in range(bsize)]
        if gtype == 'var':
            gband_varimg += [[10] for k in range(bsize)]
        else:
            gband_varimg += [[np.nan] for k in range(bsize)]
        mid = (st_binID + ed_binID)/2
        gname = gID_Gband[gID]['name']
        gtick_locs.append(mid)
        gtick_labels.append(gname)


    # set xtick labels along chromosome
    xtick_locs, xtick_labels = [], []
    for binID in sorted(binID_Gvalue.keys()):
        pos = bin_size*binID + bin_size/2
        Mb_pos = int(round(float(pos)/(10**6)))
        label = str(Mb_pos)
        if label not in xtick_labels:
            xtick_locs.append(binID)
            xtick_labels.append(label)


    # correation score VS feature
    binIDs = sorted(binID_Gvalue.keys())
    fig, axes = plt.subplots(nrows=2*len(feature_names), ncols=len(agent_list), gridspec_kw = {'hspace':0.1, 'wspace':0.12})
    for i in range(len(agent_list)):
        agent = agent_list[i]
        #X = name_sig[agent]
        X = [name_binID_mean[agent][binID] for binID in binIDs]
        for j in range(len(feature_names)):
            feature_name = feature_names[j]
            #Y = name_sig[feature_name]
            Y = [name_binID_mean[feature_name][binID] for binID in binIDs]
            corr = scipy.stats.spearmanr(X, Y)[0]

            axes[2*j, i].plot(X, Y, 'k.', markersize=2, alpha=0.5)
            
            if i > 0:
                axes[2*j, i].tick_params(axis='y', which='both', labelleft=False)
            if i == 0:
                axes[2*j, i].set_ylabel('Gene expression')
            axes[2*j, i].set_title(agent+'$^{%s}$' % (agent_charge[agent]), weight='bold')
            axes[2*j, i].set_aspect(np.diff(axes[2*j, i].get_xlim())/np.diff(axes[2*j, i].get_ylim()))

            matrix = np.zeros((len(agent_list), len(agent_list)))
            matrix[:] = corr
            img = axes[2*j+1, i].imshow(matrix, cmap="jet_r", vmin=-0.8, vmax=0, origin='lower')
            if abs(corr) < 0.3 or abs(corr) > 0.7:
                color = "white"
            else:
                color = "black"
            axes[2*j+1, i].text(len(agent_list)/2, len(agent_list)/2, str(round(corr,2)), ha="center", va="center", fontsize=10, color=color, weight='bold')
            axes[2*j+1, i].set_xlim([0, len(agent_list)-1])
            axes[2*j+1, i].set_ylim([0, len(agent_list)-1])
            axes[2*j+1, i].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    #plt.subplots_adjust(wspace=0.1, hspace=0.1)
    cbar=fig.colorbar(img, ax=axes, location='right', shrink=0.5)
    cbar.ax.set_ylabel('Spearman correlation', rotation=-90, va="bottom")
    plt.suptitle("Corrrelation with gene expression (1 Mb bin)")
    #plt.tight_layout()
    #fig.text(0.04, 0.5, 'Gene expression', va='center', rotation='vertical')
    plt.show()
    plt.close()

        
    # smoothing the binned data by sliding window average
    name_sig = {}
    for name in names:
        binID_mean = name_binID_mean[name]
        sig = []
        for binID in range(len(gband_img)):
            #size = binID_size[binID]
            size = 1
            try:
                sig += [binID_mean[binID] for k in range(size)]
            except:
                sig += [np.nan for k in range(size)]
        if name == 'eigen':
            sig = statis.slow_moving_average2(sig, int(4*i/5.0 + 1))
        else:
            sig = statis.slow_moving_average2(sig, blur_win)
        name_sig[name] = sig

    ## rescale the target signals
    #for name in agent_list:
    #    old_st, old_ed = agent_minmax[name]
    #    name_sig[name] = rescale(name_sig[name], old_st, old_ed, 0, 1)
        

    # plot genome cartoon and genome-wide data (target VS feature)
    if False:
        for i in range(len(feature_names)):
            fig = plt.figure(figsize=(15,5))
            ax1 = plt.subplot(211)
            ax2 = plt.subplot(212)

            for j in range(len(agent_list)):
                agent = agent_list[j]
                #ax1.plot(name_sig[agent], '#1f77b4', alpha=1)
                ax1.plot(name_sig[agent], agent_color[agent], alpha=0.5, label=agent)

            ax1.set_xticks(xtick_locs[::10])
            ax1.set_xticklabels(xtick_labels[::10])
            ax1.set_xlabel("Position (Mb)")
            ax1.set_ylabel("Condensability (Rescaled)", color='blue')
            ax1.tick_params('y', colors='blue')
            ax1.set_xlim([0, len(gband_img)+1])
            ax1.set_ylim([0, 1])
            ax1.legend()

            ax1p = ax1.twinx()
            feature_name = feature_names[i]
            ax1p.plot(name_sig[feature_name], '#d62728', alpha=0.2)
            #ax1p.plot(name_sig[feature_name], 'tab:orange', alpha=0.8)
            #ax1p.plot(np.log(-1*np.asarray(name_sig[feature_name])), '#d62728', alpha=0.5)
            ax1p.set_ylabel(feature_name, color='r')
            #ax1p.set_ylabel('Eigenvector', color='orangered')
            #ax1p.tick_params('y', colors='#d62728')
            #ax1p.tick_params('y', colors='orangered')
            ax1p.set_ylim([-0.1, 2.0])

            for gID in sorted(gID_binterval):
                st, ed = gID_binterval[gID]
                gtype = gID_Gband[gID]['type']
                if gtype == 'pos':
                    ax1.axvspan(st, ed-1, alpha=0.15, color='grey')
                if gtype == 'acen' or gtype == 'var':
                    ax1.axvspan(st, ed-1, alpha=1, color='white', zorder=10)
                    ax1p.axvspan(st, ed-1, alpha=1, color='white', zorder=10)

            ax2.imshow(np.transpose(gband_img), cmap='Greys', aspect=0.3/aspect)
            ax2.imshow(np.transpose(gband_cenimg), cmap ='Reds', vmin=0, vmax=20, aspect=0.3/aspect)
            ax2.imshow(np.transpose(gband_varimg), cmap ='Purples', vmin=0, vmax=20, aspect=0.3/aspect)
            ax2.set_yticks([])
            ax2.set_xticks(gtick_locs)
            ax2.set_xticklabels(gtick_labels)
            ax2.tick_params(axis="x", labelsize=5, rotation=90)
            ax2.set_xlim([0, len(gband_img)+1])

            plt.tight_layout()
            #plt.savefig("Gwide_" + chr_choice + '_' + target_name + '_' + feature_name + ".png", bbox_inches='tight', dpi=1000)
            plt.show()
            plt.close()

    # plot genome cartoon and genome-wide data (vertical)
    if False:
        figwidth = 1 + 2*len(agent_list)
        figheight = round(10*float(len(gband_img))/9970, 1)
        fig, axes = plt.subplots(nrows=1, ncols=1 + len(agent_list), figsize=(figwidth, figheight), sharey=True, gridspec_kw = {'hspace':0.02, 'wspace':0.12, 'width_ratios':[1] + [8]*len(agent_list)})

        for i in range(len(agent_list)):
            agent = agent_list[i]
            axes[i+1].plot(name_sig[agent], range(len(name_sig[agent])), '#1f77b4', alpha=1)

            ax1p = axes[i+1].twiny()
            for j in range(len(feature_names)):
                feature_name = feature_names[j]
                ax1p.plot(name_sig[feature_name], range(len(name_sig[feature_name])), '#d62728', alpha=0.25)

            for gID in sorted(gID_binterval):
                gtype = gID_Gband[gID]['type']
                if gtype == 'acen' or gtype == 'var':
                    st, ed = gID_binterval[gID]
                    axes[i+1].axhspan(st, ed-1, alpha=1, color='white', zorder=10)
                    ax1p.axhspan(st, ed-1, alpha=1, color='white', zorder=10)

            min,max = agent_minmax[agent]
            axes[i+1].set_xlim([min - 0.05*(max-min), max + 0.05*(max-min)])
            axes[i+1].spines['top'].set_visible(False)
            #axes[i+1].spines['bottom'].set_visible(False)
            axes[i+1].spines['left'].set_visible(False)
            axes[i+1].spines['right'].set_visible(False)
            axes[i+1].tick_params(top='off', bottom='on', left='off', right='off', labelbottom='on')
            axes[i+1].tick_params('x', colors='blue')
            axes[i+1].set_xticks([min, max])
            axes[i+1].set_xticklabels([str(min), str(max)], rotation=45, fontsize=7, color='blue')
            
            ax1p.set_xlabel(agent+'$^{%s}$' % (agent_charge[agent]),
                            rotation=40, fontsize=12, labelpad=30, ha='right')
            ax1p.set_xlim([-0.1, 2.0])
            ax1p.spines['top'].set_visible(False)
            ax1p.spines['bottom'].set_visible(False)
            ax1p.spines['left'].set_visible(False)
            ax1p.spines['right'].set_visible(False)
            ax1p.tick_params(top='off', bottom='off', left='off', right='off', labeltop='off')

        axes[0].imshow(gband_img, cmap='Greys', aspect='auto')
        axes[0].imshow(gband_cenimg, cmap ='Reds', vmin=0, vmax=20, aspect='auto')
        axes[0].imshow(gband_varimg, cmap ='Purples', vmin=0, vmax=20, aspect='auto')
        axes[0].set_xticks([])
        axes[0].set_yticks([])
        
        #plt.tight_layout()
        plt.savefig("Gwide_" + chr_choice + '_all_agents_' + feature_name + ".png", bbox_inches='tight')
        plt.show()
        plt.close()

    # plot genome cartoon and genome-wide data (horizontal)
    if True:
        figwidth = round(12*float(len(gband_img))/9970, 1)
        figheight = 1 + 2*len(agent_list)
        fig, axes = plt.subplots(nrows=1+len(agent_list), ncols=1, figsize=(figwidth, figheight), sharex=True, gridspec_kw = {'hspace':0.12, 'wspace':0.12, 'height_ratios':[7]*len(agent_list)+[1]})

        for i in range(len(agent_list)):
            agent = agent_list[i]
            axes[i].plot(name_sig[agent], '#1f77b4', alpha=1)

            ax1p = axes[i].twinx()
            for j in range(len(feature_names)):
                feature_name = feature_names[j]
                ax1p.plot(name_sig[feature_name], '#d62728', alpha=0.3)

            for gID in sorted(gID_binterval):
                gtype = gID_Gband[gID]['type']
                if gtype == 'acen' or gtype == 'var':
                    st, ed = gID_binterval[gID]
                    axes[i].axvspan(st, ed-1, alpha=1, color='white', zorder=10)
                    ax1p.axvspan(st, ed-1, alpha=1, color='white', zorder=10)

            min,max = agent_minmax[agent]
            #axes[i].set_ylim([min - 0.05*(max-min), max + 0.05*(max-min)])
            axes[i].set_ylim([min, max])
            axes[i].spines['top'].set_visible(False)
            axes[i].spines['bottom'].set_visible(False)
            axes[i].spines['left'].set_visible(False)
            axes[i].spines['right'].set_visible(False)
            axes[i].tick_params(top='off', bottom='off', left='on', right='off', labelbottom='off')
            axes[i].tick_params('y', colors='blue')
            axes[i].set_yticks([min + 0.2*(max-min), max - 0.2*(max-min)])
            axes[i].set_yticklabels([str(min + 0.2*(max-min)), str(max - 0.2*(max-min))], fontsize=7, color='blue')
            axes[i].set_ylabel(agent+'$^{%s}$\n(t#%d)' % (agent_charge[agent], agent_titrnum[agent]),
                               rotation=0, fontsize=13, labelpad=20, ha='center', va='center',
                               weight='bold')
            
            ax1p.set_ylim([-0.1, 2.0])
            ax1p.spines['top'].set_visible(False)
            ax1p.spines['bottom'].set_visible(False)
            ax1p.spines['left'].set_visible(False)
            ax1p.spines['right'].set_visible(False)
            ax1p.tick_params(top='off', bottom='off', left='off', right='off', labelright='off')


        axes[len(agent_list)].imshow(np.transpose(gband_img), cmap='Greys', aspect='auto')
        axes[len(agent_list)].imshow(np.transpose(gband_cenimg), cmap ='Reds', vmin=0, vmax=20, aspect='auto')
        axes[len(agent_list)].imshow(np.transpose(gband_varimg), cmap ='Purples', vmin=0, vmax=20, aspect='auto')
        axes[len(agent_list)].set_xticks([])
        axes[len(agent_list)].set_yticks([])
        
        #plt.tight_layout()
        #plt.savefig("Gwide_" + chr_choice + '_all_agents_' + feature_name + ".png", bbox_inches='tight')
        #plt.show()
        plt.close()
