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

#chr_choices = ['chr' + str(i+1) for i in range(17, 23)]
chr_choices = ['chr1']

for chr_choice in chr_choices:

    # read G-band information and build interval-dictionary
    path = ""
    #path = "./data/"
    gID_Gband = read_Gband(path+"Gband_information.txt", chr_choices=[chr_choice])[chr_choice]

    # draw the continuous singal in the 1Mb resolution
    path = ""
    #path = "./data/"
    #field_ID_value = load_file.read_tabular_file (path + "hg19_" + chr_choice + "_167win25step_anot.cn", mode='col', jump=10)
    field_ID_value = load_file.read_tabular_file (path + "H1_NCP_sp_" + chr_choice + "_167win25step_anot.cn", mode='col', jump=10)

    ID_pos = field_ID_value['PhysicalPosition']
    ID_AT = field_ID_value['ATcontent']
    ID_k27ac = field_ID_value['k27ac']
    ID_score1 = field_ID_value["data/sp_spd_tests_detail/sp7"]
    ID_score2 = field_ID_value["data/sp_spd_tests_detail/sp8"]

    #geneID_field_values, field_geneID_values = load_file.read_GTF (path + "Homo_sapiens.GRCh37.87.gtf", chr_list=[chr_choice], mode="both")
    #geneID_RPKM = load_file.read_RPKM (path+"GSE63124_all_gene_raw_readcounts.txt", path+"Homo_sapiens.GRCh37.87.gtf", chr_choice)

    geneID_field_values, field_geneID_values = load_file.read_GTF ("ENCFF159KBI.gtf", chr_list=[chr_choice], mode="both")
    geneID_RPKM = read_tsv("ENCFF174OMR.tsv")

    #geneID_pos = {}
    #for geneID in geneID_field_values:
    #    pos = geneID_field_values[geneID]['TSS']
    #    geneID_pos[geneID] = pos

    geneID_pos = {}
    for geneID in geneID_field_values:
        try:
            pos = geneID_field_values[geneID]['TSS']
            geneID_pos[geneID] = pos
        except:
            continue


    #names = ['data/sp_spd_tests_detail/sp7', 'ATcontent', 'CpGNumber', 'k27ac', 'k9ac', 'k4me3', 'k36me3_2', 'k9me2_2', 'k9me3_2', 'k27me3a_2']
    #names = ['data/sp_spd_tests_detail/sp7', 'ATcontent', 'CpGNumber', 'k27ac', 'k9ac', 'k4me3', 'k36me3', 'k9me2', 'k9me3', 'k27me3']
    #cmap_list = ['rainbow', 'hot_r', 'viridis_r', 'YlOrRd', 'Purples', 'Oranges', 'Greens', 'Blues', 'YlGnBu', 'pink_r']
    names = ['data/sp_spd_tests_detail/sp7', 'k9ac']
    cmap_list = ['rainbow']

    i = 20
    bin_size = int(0.5*(10**6) / i)
    #blur_win = int(2*i + 1)
    blur_win = int(4*i + 1)
    #blur_win = int(6*i + 1)
    #blur_win = int(10*i + 1)

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

    temp = {}
    for ID in ID_pos:
        binID = int(ID_pos[ID]) / int(bin_size)
        temp[binID] = 0

    name_binID_mean['Gene density'] = copy.deepcopy(temp)
    name_binID_count['Gene density'] = copy.deepcopy(temp)
    name_binID_mean['Gene activity'] = copy.deepcopy(temp)
    name_binID_count['Gene activity'] = copy.deepcopy(temp)

    min_RPKM = min(geneID_RPKM.values())
    for geneID in geneID_pos:
        binID = int(geneID_pos[geneID]) / int(bin_size)
        name_binID_mean['Gene density'][binID] += 1.0
        try:
            #name_binID_mean['Gene activity'][binID] += geneID_RPKM[geneID]
            name_binID_mean['Gene activity'][binID] += np.log2(geneID_RPKM[geneID] - min_RPKM + 1)
            #name_binID_mean['Gene activity'][binID] += np.log2(geneID_RPKM[geneID])
        except:
            #name_binID_mean['Gene density'][binID] += 0.0
            name_binID_mean['Gene activity'][binID] += np.nan
        name_binID_count['Gene activity'][binID] += 1

    for binID in name_binID_mean['Gene activity']:
        if name_binID_count['Gene activity'][binID] <= 0:
            continue
        name_binID_mean['Gene activity'][binID] /= name_binID_count['Gene activity'][binID] 

    names.append('Gene density')
    cmap_list.append('jet')
    names.append('Gene activity')
    cmap_list.append('jet')

    for binID in name_binID_mean['k9ac']:
        try:
            name_binID_mean['k9ac'][binID] /= name_binID_mean['Gene density'][binID]
        except:
            name_binID_mean['k9ac'][binID] = np.nan

    # G-banding ideogram
    gID_size = {}
    gband_img = []
    gband_cenimg, gband_varimg = [], []
    ytick_locs, ytick_labels = [], []
    prev_size = 0
    for gID in sorted(gID_Gband.keys()):
        st, ed = gID_Gband[gID]['interval']
        size = int(round(float(ed - st)/(bin_size)))
        gID_size[gID] = size
        value = gID_Gband[gID]['value']
        gband_img += [[value] for k in range(size)]
        if gID_Gband[gID]['type'] == 'var':
            gband_varimg += [[10] for k in range(size)]
        else:
            gband_varimg += [[np.nan] for k in range(size)]
        if gID_Gband[gID]['type'] == 'acen':
            gband_cenimg += [[10] for k in range(size)]
        else:
            gband_cenimg += [[np.nan] for k in range(size)]
        gname = gID_Gband[gID]['name']
        mid = prev_size + size/2
        ytick_locs.append(mid)
        ytick_labels.append(gname)
        prev_size += size

    name_sig = {}
    name_img = {}
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
        sig = statis.slow_moving_average2(sig, blur_win)
        name_sig[name] = sig
        name_img[name] = [[value] for value in sig]

    aspect = (0.05*bin_size) / (10**6)
    """
    fig = plt.figure()
    aspect = (0.05*bin_size) / (10**6)
    plot_num = len(names) + 1
    plt.subplot(1, plot_num, 1)
    plt.imshow(gband_cenimg, cmap ='Reds', vmin=0, vmax=20, aspect=aspect)
    plt.imshow(gband_varimg, cmap ='Purples', vmin=0, vmax=20, aspect=aspect)
    #cmap = plt.cm.Greys
    #cmap.set_bad((1, 0, 0, 1))
    plt.imshow(gband_img, cmap='Greys', aspect=aspect)
    plt.xticks([])
    plt.yticks(ytick_locs, ytick_labels)
    ax = plt.gca()
    ax.tick_params(axis="y", labelsize=5)
    plt.xlabel('G-banding', rotation=45)
    for i in range(len(names)):
        name = names[i]
        img = name_img[name]
        plt.subplot(1, plot_num, i+2)
        #if name in ['k9me2', 'k9me3', 'k27me3a', 'k9me3_2']:
        #    img = np.log2(np.asarray(img))
        im = plt.imshow(img, cmap=cmap_list[i], aspect=aspect)
        plt.colorbar(im, ax=[plt.gca()], shrink=0.2, anchor=(0,0.95))
        plt.xticks([])
        plt.yticks([])
        if name == 'data/sp_spd_tests_detail/sp7':
            name = "Condensability"
        name = name.split('_')[0]
        if name.startswith('k27me3'):
            name = 'k27me3'
        plt.xlabel(name, rotation=45)
    #plt.show()
    plt.close()
    """
    
    cond_sig = name_sig['data/sp_spd_tests_detail/sp7']
    xtick_locs = []
    xtick_labels = []
    for i in range(len(cond_sig)):
        pos = bin_size*i + bin_size/2
        Mb_pos = int(round(float(pos)/(10**6)))
        label = str(Mb_pos)
        if label not in xtick_labels:
            xtick_locs.append(i)
            xtick_labels.append(label)

    """
    fig = plt.figure(figsize=(15,5))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
    ax1.plot(cond_sig, alpha=1)
    #color_list = np.linspace(0, 1, num=len(keys))
    #cmap = cm.get_cmap('Spectral')
    #cond_sig = np.asarray(cond_sig)
    #ax1.scatter(range(len(cond_sig)), cond_sig, c=cond_sig, s=3, cmap='RdYlBu_r')
    pt = 0
    for gID in sorted(gID_size):
        size = gID_size[gID]
        st, ed = pt, pt+size
        #ax1.plot(range(pt, pt+size), cond_sig[pt:pt+size])
        if gID_Gband[gID]['type'] == 'pos':
            ax1.axvspan(st, ed-1, alpha=0.15, color='grey')
        pt += size
    assert pt == len(cond_sig)
    ax1.set_xticks(xtick_locs[::10] + [xtick_locs[-1]])
    ax1.set_xticklabels(xtick_labels[::10] + [xtick_labels[-1]])
    ax1.set_xlabel("Position (Mb)")
    ax1.set_ylabel("Condensability (A.U.)")
    ax1.set_xlim([0, len(gband_img)+1])
    #ax1p = ax1.twinx()
    #ax1p.plot(name_sig['ATcontent'], 'r', alpha=0.5)
    #ax1p.set_ylabel('AT content', color='r')
    #ax1p.tick_params('y', colors='r')
    ax2.imshow(np.transpose(gband_img), cmap='Greys', aspect=0.3/aspect)
    ax2.imshow(np.transpose(gband_cenimg), cmap ='Reds', vmin=0, vmax=20, aspect=0.3/aspect)
    ax2.imshow(np.transpose(gband_varimg), cmap ='Purples', vmin=0, vmax=20, aspect=0.3/aspect)
    ax2.set_yticks([])
    ax2.set_xticks(ytick_locs)
    ax2.set_xticklabels(ytick_labels)
    ax2.tick_params(axis="x", labelsize=5, rotation=90)
    ax2.set_xlim([0, len(gband_img)+1])
    plt.tight_layout()
    #plt.savefig("Gwide_" + 'conden' + ".png", bbox_inches='tight', dpi=1000)
    #plt.show()
    plt.close()
    """
    
    for name in name_sig:
        if name == 'data/sp_spd_tests_detail/sp7':
            continue
        fig = plt.figure(figsize=(15,5))
        ax1 = plt.subplot(211)
        ax2 = plt.subplot(212)
        ax1.plot(name_sig[name], '#1f77b4', alpha=1)
        #ax1.plot(cond_sig, '#1f77b4', alpha=1)
        pt = 0
        for gID in sorted(gID_size):
            size = gID_size[gID]
            st, ed = pt, pt+size
            if gID_Gband[gID]['type'] == 'pos':
                ax1.axvspan(st, ed-1, alpha=0.15, color='grey')
            pt += size
        assert pt == len(cond_sig)
        ax1.set_xticks(xtick_locs[::10])
        ax1.set_xticklabels(xtick_labels[::10])
        #ax1.set_xticks(xtick_locs[::10] + [xtick_locs[-1]])
        #ax1.set_xticklabels(xtick_labels[::10] + [xtick_labels[-1]])
        ax1.set_xlabel("Position (Mb)")
        #ax1.set_ylabel("Condensability (A.U.)", color='#1f77b4')
        #ax1.set_ylabel("Condensability (A.U.)", color='blue')
        ax1.set_ylabel(name, color='blue')
        ax1.tick_params('y', colors='blue')
        ax1.set_xlim([0, len(gband_img)+1])
        ax1p = ax1.twinx()
        ax1p.plot(name_sig['Gene activity'], '#d62728', alpha=0.5)
        ax1p.set_ylabel('Expression ($log$ RPKM)', color='r')

        """
        if name == 'Gene activity':
            #ax1p.scatter(range(len(name_sig[name])), name_sig[name], c='#d62728', s=1.5, alpha=0.1)
            ax1p.plot(name_sig[name], '#d62728', alpha=0.5)
            ax1p.set_ylabel('Expression ($log$ RPKM)', color='r')
        else:
            ax1p.plot(name_sig[name], '#d62728', alpha=0.5)
            ax1p.set_ylabel(name, color='#d62728')
        """
        #ax1p.plot(name_sig[name], '#d62728', alpha=0.5)
        #ax1p.set_ylabel(name, color='#d62728')
        ax1p.tick_params('y', colors='#d62728')
        ax2.imshow(np.transpose(gband_img), cmap='Greys', aspect=0.3/aspect)
        ax2.imshow(np.transpose(gband_cenimg), cmap ='Reds', vmin=0, vmax=20, aspect=0.3/aspect)
        ax2.imshow(np.transpose(gband_varimg), cmap ='Purples', vmin=0, vmax=20, aspect=0.3/aspect)
        ax2.set_yticks([])
        ax2.set_xticks(ytick_locs)
        ax2.set_xticklabels(ytick_labels)
        ax2.tick_params(axis="x", labelsize=5, rotation=90)
        ax2.set_xlim([0, len(gband_img)+1])
        plt.tight_layout()
        plt.savefig("Gwide_" + chr_choice + '_' + name + ".png", bbox_inches='tight', dpi=1000)
        #plt.show()
        plt.close()


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

    eigenBinID_value, eigenBinID_interval = read_eigenfile(path+"eigen_WT_50kb.txt", bin_size=50000)

    BinID_eigens = {}
    for value, interval in zip(eigenBinID_value, eigenBinID_interval):
        st, ed = interval
        st_BinID, ed_BinID = st / bin_size, ed / bin_size
        for BinID in range(st_BinID, ed_BinID):
            if BinID not in BinID_eigens:
                BinID_eigens[BinID] = []
            BinID_eigens[BinID].append(value)

    BinID_eigen = {}
    for BinID in range(len(gband_img)):
        try:
            BinID_eigen[BinID] = np.mean(BinID_eigens[BinID])
        except:
            BinID_eigen[BinID] = np.nan

    eigen_sig = [BinID_eigen[BinID] for BinID in sorted(BinID_eigen.keys())]
    A_sig, B_sig = [], []
    for value in eigen_sig:
        if value >= 0:
            A_sig.append(value)
            B_sig.append(0.0)
        else:
            A_sig.append(0.0)
            B_sig.append(value)

    fig = plt.figure(figsize=(15,5))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
    #ax1.plot(cond_sig, '#1f77b4', alpha=1)
    ax1.plot(cond_sig, 'k', alpha=1)
    pt = 0
    for gID in sorted(gID_size):
        size = gID_size[gID]
        st, ed = pt, pt+size
        if gID_Gband[gID]['type'] == 'pos':
            ax1.axvspan(st, ed-1, alpha=0.15, color='grey')
        pt += size
    assert pt == len(cond_sig)
    ax1.set_xticks(xtick_locs[::10] + [xtick_locs[-1]])
    ax1.set_xticklabels(xtick_labels[::10] + [xtick_labels[-1]])
    ax1.set_xlabel("Position (Mb)")
    #ax1.set_ylabel("Condensability (A.U.)", color='#1f77b4')
    ax1.set_ylabel("Condensability (A.U.)", color='k')
    ax1.tick_params('y', colors='k')
    ax1.set_xlim([0, len(gband_img)+1])
    ax1p = ax1.twinx()
    #ax1p.plot(eigen_sig, '#d62728', alpha=0.5)
    #ax1p.plot(eigen_sig, 'k', alpha=0.15)
    ax1p.fill_between(range(len(A_sig)), [0.0]*len(A_sig), A_sig, alpha=0.2, color='C1')
    ax1p.fill_between(range(len(B_sig)), B_sig, [0.0]*len(B_sig), alpha=0.2, color='C0')
    #ax1p.plot(A_sig, 'red', alpha=0.5)
    #ax1p.plot(B_sig, 'blue', alpha=0.5)
    ax1p.set_ylabel('Eigenvector', color='k')
    ax1p.tick_params('k', colors='#d62728')
    ax2.imshow(np.transpose(gband_img), cmap='Greys', aspect=0.3/aspect)
    ax2.imshow(np.transpose(gband_cenimg), cmap ='Reds', vmin=0, vmax=20, aspect=0.3/aspect)
    ax2.imshow(np.transpose(gband_varimg), cmap ='Purples', vmin=0, vmax=20, aspect=0.3/aspect)
    ax2.set_yticks([])
    ax2.set_xticks(ytick_locs)
    ax2.set_xticklabels(ytick_labels)
    ax2.tick_params(axis="x", labelsize=5, rotation=90)
    ax2.set_xlim([0, len(gband_img)+1])
    plt.tight_layout()
    #plt.savefig("Gwide_eigen" + ".png", bbox_inches='tight', dpi=1000)
    #plt.show()
    plt.close()

# Check TADs



    
    


"""
# divide genome into G-band bins

gID_interval = {}
for gID in gID_Gband:
    interval = gID_Gband[gID]['interval']
    gID_interval[gID] = interval

ginterval_dict = Interval_dict.double_hash(gID_interval, 100000, 250000000)

# read annotation file
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+"hg19_chr1_NCP_anot.cn")

ID_score1 = name_ID_value['data/sp_spd_tests_detail/sp7']
ID_AT = name_ID_value['ATcontent']
ID_CpG = name_ID_value['CpGNumber']
ID_me = name_ID_value['meGCNumber']

ID_mefrac = {}
for ID in ID_CpG:
    CG = ID_CpG[ID]
    if CG <= 0:
        continue
    me = ID_me[ID]
    mefrac = float(me) / (2*CG)
    ID_mefrac[ID] = mefrac
name_ID_value['meCpG density'] = ID_mefrac

name_gID_values = {}
for ID in ID_pos:
    pos = ID_pos[ID]
    gIDs = ginterval_dict.find(pos)
    if len(gIDs) <= 0:
        continue
    for name in name_ID_value:
        if name not in name_gID_values:
            name_gID_values[name] = {}
        try:
            value = name_ID_value[name][ID]
        except:
            continue
        if np.isnan(value):
            continue
        for gID in gIDs:
            if gID not in name_gID_values[name]:
                name_gID_values[name][gID] = []
            name_gID_values[name][gID].append(value)

print "sorting is done"

name_gID_mean = {}
for name in name_gID_values:
    if name not in name_gID_mean:
        name_gID_mean[name] = {}
    for gID in name_gID_values[name]:
        assert gID not in name_gID_mean[name]
        name_gID_mean[name][gID] = np.mean(name_gID_values[name][gID])

# G-banding ideogram
gID_size = {}
gband_img = []
gband_cenimg = []
gband_varimg = []
ytick_locs, ytick_labels = [], []
prev_size = 0
for gID in sorted(gID_Gband.keys()):
    st, ed = gID_Gband[gID]['interval']
    size = int(round(float(ed - st)/(10**6)))
    gID_size[gID] = size
    value = gID_Gband[gID]['value']
    gband_img += [[value] for k in range(size)]
    if gID_Gband[gID]['type'] == 'var':
        gband_varimg += [[10] for k in range(size)]
    else:
        gband_varimg += [[np.nan] for k in range(size)]
    if gID_Gband[gID]['type'] == 'acen':
        gband_cenimg += [[10] for k in range(size)]
    else:
        gband_cenimg += [[np.nan] for k in range(size)]
    gname = gID_Gband[gID]['name']
    mid = prev_size + size/2
    ytick_locs.append(mid)
    ytick_labels.append(gname)
    prev_size += size

#names = ['data/sp_spd_tests_detail/sp7', 'data/sp_spd_tests_detail/sp8', 'ATcontent', 'CpGNumber', 'meGCNumber', 'k27ac', 'k27me3a', 'k27me3a_2', 'k27me3b', 'k27me3b_2', 'k36me3', 'k36me3_2', 'k4me3', 'k9ac', 'k9me2', 'k9me2_2', 'k9me3', 'k9me3_2']

names = ['data/sp_spd_tests_detail/sp7', 'ATcontent', 'CpGNumber', 'k27ac', 'k9ac', 'k4me3', 'k36me3_2', 'k9me2_2', 'k9me3_2', 'k27me3a_2']
cmap_list = ['rainbow', 'hot_r', 'viridis_r', 'YlOrRd', 'Purples', 'Oranges', 'Greens', 'Blues', 'YlGnBu', 'pink_r']


name_img = {}
for name in names:
    gID_mean = name_gID_mean[name]
    img = []
    for gID in sorted(gID_Gband.keys()):
        size = gID_size[gID]
        try:
            img += [[gID_mean[gID]] for k in range(size)]
        except:
            img += [[np.nan] for k in range(size)]
    name_img[name] = img

fig = plt.figure()
plot_num = len(names) + 1
plt.subplot(1, plot_num, 1)
plt.imshow(gband_cenimg, cmap ='Reds', vmin=0, vmax=20, aspect=0.05)
plt.imshow(gband_varimg, cmap ='Purples', vmin=0, vmax=20, aspect=0.05)
#cmap = plt.cm.Greys
#cmap.set_bad((1, 0, 0, 1))
plt.imshow(gband_img, cmap='Greys', aspect=0.05)
plt.xticks([])
plt.yticks(ytick_locs, ytick_labels)
ax = plt.gca()
ax.tick_params(axis="y", labelsize=5)
plt.xlabel('G-banding', rotation=45)
for i in range(len(names)):
    name = names[i]
    img = name_img[name]
    plt.subplot(1, plot_num, i+2)
    if name in ['k9me2', 'k9me3', 'k27me3a', 'k9me3_2']:
        img = np.log2(np.asarray(img))
    im = plt.imshow(img, cmap=cmap_list[i], aspect=0.05)
    plt.colorbar(im, ax=[plt.gca()], shrink=0.2, anchor=(0,0.95))
    plt.xticks([])
    plt.yticks([])
    if name == 'data/sp_spd_tests_detail/sp7':
        name = "Condensability"
    name = name.split('_')[0]
    if name.startswith('k27me3'):
        name = 'k27me3'
    plt.xlabel(name, rotation=45)
#plt.show()
plt.close()
"""
