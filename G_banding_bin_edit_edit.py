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
import matplotlib.gridspec as gridspec

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
            score = -np.log((float(cols[4+i])+1) / (control_count+1))
            if name not in name_ID_score:
                name_ID_score[name] = {}
            assert ID not in name_ID_score[name]
            name_ID_score[name][ID] = score

    return ID_pos, name_ID_score


# read A/B compartment annotation
def read_eigenfile (fname, bin_size=1000000):
    eigenBinID_value = []
    eigenBinID_interval = []
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
        eigenBinID_value.append(value)
        eigenBinID_interval.append((st,ed))
        i +=1
    return eigenBinID_value, eigenBinID_interval


#### set parameters
# set binning resolution
i = 20
bin_size = int(0.5*(10**6) / i) # binsize (unit of bp)
#bin_size = int(0.5*(10**5) / i)
#blur_win = 1
#blur_win = int(2*i + 1)
blur_win = int(4*i + 1) # sliding window (unit of bin)
#blur_win = int(6*i + 1)
#blur_win = int(10*i + 1)

# set chromosomes
#chr_choice = ['chr1']
#chr_choices = ['chr10']
#chr_choices = ['chr%d' % (i) for i in range(1, 13)]
#chr_choices = ['chr1']
#chr_choices = ['chrX', 'chrY']
chr_choices = ['chr%d' % (i) for i in range(1, 23)] + ['chrX', 'chrY']
#chr_choices = ['chr%d' % (i) for i in range(1, 4)]

# set target names and feature names
target_names = ["H1-NCP-sp-8.bam"]
feature_names = ['Gene activity']
names = target_names + feature_names

# plot mode
#plot_mode = 'targetVSfeature'
#plot_mode = 'targetVSfeature_nospace'
#plot_mode = 'targetonly_nospace'

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
    ID_pos, field_ID_value = read_bin("H1_NCP_sp_10kb_bin.cn", chr_choices=[chr_choice])

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
        

    # plot genome cartoon and genome-wide data (target VS feature)
    if False:
        for i in range(len(feature_names)):
            fig = plt.figure(figsize=(15,5))
            ax1 = plt.subplot(211)
            ax2 = plt.subplot(212)

            for j in range(len(target_names)):
                target_name = target_names[j]
                ax1.plot(name_sig[target_name], '#1f77b4', alpha=1)

            for gID in sorted(gID_binterval):
                st, ed = gID_binterval[gID]
                if gID_Gband[gID]['type'] == 'pos':
                    ax1.axvspan(st, ed-1, alpha=0.15, color='grey')

            ax1.set_xticks(xtick_locs[::10])
            ax1.set_xticklabels(xtick_labels[::10])
            ax1.set_xlabel("Position (Mb)")
            ax1.set_ylabel(target_name, color='blue')
            ax1.tick_params('y', colors='blue')
            ax1.set_xlim([0, len(gband_img)+1])
            ax1.set_ylim([-0.4, 0.6])

            ax1p = ax1.twinx()
            feature_name = feature_names[i]
            ax1p.plot(name_sig[feature_name], '#d62728', alpha=0.5)
            #ax1p.plot(name_sig[feature_name], 'tab:orange', alpha=0.8)
            #ax1p.plot(np.log(-1*np.asarray(name_sig[feature_name])), '#d62728', alpha=0.5)
            ax1p.set_ylabel(feature_name, color='r')
            #ax1p.set_ylabel('Eigenvector', color='orangered')
            #ax1p.tick_params('y', colors='#d62728')
            #ax1p.tick_params('y', colors='orangered')
            ax1p.set_ylim([-0.1, 2.0])

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
            plt.savefig("Gwide_" + chr_choice + '_' + target_name + '_' + feature_name + ".svg", format='svg', bbox_inches='tight')
            plt.show()
            plt.close()


    # plot genome cartoon and genome-wide data (no space, horizontal)
    if False:
        figwidth = round(10*float(len(gband_img))/9970, 1)
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(figwidth,1.5), sharex=True, gridspec_kw = {'hspace':0.02, 'wspace':None, 'height_ratios':[6,1]})

        for j in range(len(target_names)):
            target_name = target_names[j]
            axes[0].plot(name_sig[target_name], '#1f77b4', alpha=1)

        ax0p = axes[0].twinx()
        for j in range(len(feature_names)):
            feature_name = feature_names[j]
            ax0p.plot(name_sig[feature_name], '#d62728', alpha=0.25)

        for gID in sorted(gID_binterval):
            gtype = gID_Gband[gID]['type']
            if gtype == 'acen' or gtype == 'var':
                st, ed = gID_binterval[gID]
                axes[0].axvspan(st, ed-1, alpha=1, color='white', zorder=10)
                ax0p.axvspan(st, ed-1, alpha=1, color='white', zorder=10)


        axes[0].set_ylabel(chr_choice, rotation=90+40, fontsize=12, labelpad=10)
        axes[0].set_ylim([-0.5, 0.6])
        axes[0].spines['top'].set_visible(False)
        axes[0].spines['bottom'].set_visible(False)
        axes[0].spines['left'].set_visible(False)
        axes[0].spines['right'].set_visible(False)
        axes[0].tick_params(top='off', bottom='off', left='off', right='off', labelleft='off', labelbottom='off')

        ax0p.set_ylim([-0.1, 2.0])
        ax0p.spines['top'].set_visible(False)
        ax0p.spines['bottom'].set_visible(False)
        ax0p.spines['left'].set_visible(False)
        ax0p.spines['right'].set_visible(False)
        ax0p.tick_params(top='off', bottom='off', left='off', right='off', labelright='off', labelbottom='off')

        axes[1].imshow(np.transpose(gband_img), cmap='Greys', aspect='auto')
        axes[1].imshow(np.transpose(gband_cenimg), cmap ='Reds', vmin=0, vmax=20, aspect='auto')
        axes[1].imshow(np.transpose(gband_varimg), cmap ='Purples', vmin=0, vmax=20, aspect='auto')
        axes[1].set_xticks([])
        axes[1].set_yticks([])
        
        #plt.tight_layout()
        #plt.savefig("Gwide_" + chr_choice + '_' + target_name + '_' + feature_name + ".png", bbox_inches='tight', dpi=500)
        #plt.xticks(rotation=90)
        plt.show()
        plt.close()

    # plot genome cartoon and genome-wide data (no space, vertical)
    if True:
        figheight = round(8*float(len(gband_img))/9970, 1)
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(2, figheight), sharey=True, gridspec_kw = {'hspace':0.02, 'wspace':0.02, 'width_ratios':[1,8]})

        for j in range(len(target_names)):
            target_name = target_names[j]
            axes[1].plot(name_sig[target_name], range(len(name_sig[target_name])), '#1f77b4', alpha=1)

        ax1p = axes[1].twiny()
        for j in range(len(feature_names)):
            feature_name = feature_names[j]
            ax1p.plot(name_sig[feature_name], range(len(name_sig[feature_name])), '#d62728', alpha=0.25)

        for gID in sorted(gID_binterval):
            gtype = gID_Gband[gID]['type']
            if gtype == 'acen' or gtype == 'var':
                st, ed = gID_binterval[gID]
                axes[1].axhspan(st, ed-1, alpha=1, color='white', zorder=10)
                ax1p.axhspan(st, ed-1, alpha=1, color='white', zorder=10)

        axes[1].set_xlim([-0.5, 0.6])
        axes[1].spines['top'].set_visible(False)
        axes[1].spines['bottom'].set_visible(False)
        axes[1].spines['left'].set_visible(False)
        axes[1].spines['right'].set_visible(False)
        axes[1].tick_params(top='off', bottom='off', left='off', right='off', labelbottom='off')

        ax1p.set_xlabel(chr_choice, rotation=40, fontsize=15, labelpad=23, ha='right')
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
        #plt.savefig("Gwide_" + chr_choice + '_' + target_name + '_' + feature_name + ".png", bbox_inches='tight')
        plt.savefig("Gwide_" + chr_choice + '_' + target_name + '_' + feature_name + ".svg", format='svg', bbox_inches='tight')

        #plt.show()
        plt.close()
