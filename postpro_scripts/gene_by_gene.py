import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import load_file
import statis

def single_plot (profiles, offset=-1000, xtick_loc_name=None, xlabel='Distance from TSS (bp)', ylabel='Nucleosome Occupancy', names=None, title="", note=""):
    for profile in profiles:
        pad_len = int(len(profile)*0.1)
        profile[:pad_len] = [np.NaN]*pad_len
        profile[len(profile)-pad_len:] = [np.NaN]*pad_len
    fig, ax = plt.subplots(figsize=(10,5))
    #fig, ax1 = plt.subplots()
    X = [ i + offset for i in range(len(profiles[0]))]
    for k in range(len(profiles)):
        profile = profiles[k]
        if names !=None:
            ax.plot(X, profile, label=names[k])
        else:
            ax.plot(X. profile)
    if xtick_loc_name:
        xtick_locs, xtick_names = xtick_loc_name
        ax.set_xticks(xtick_locs)
        ax.set_xticklabels(xtick_names)
    fig.tight_layout()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.title(title)
    if names != None:
        plt.legend()
    plt.savefig("single_" + note + ".png",bbox_inches='tight')
    #plt.show()
    plt.close()


def pair_plot (profile1, profile2, offset=-1000, xtick_loc_name=None, xlabel='Distance from TSS (bp)', ylabel1='Condensability (A.U.)', ylabel2="", title="", note=""):
    assert len(profile1) == len(profile2)
    pad_len = int(len(profile1)*0.1)
    profile1[:pad_len] = [np.NaN]*pad_len
    profile1[len(profile1)-pad_len:] = [np.NaN]*pad_len
    profile2[:pad_len] = [np.NaN]*pad_len
    profile2[len(profile2)-pad_len:] = [np.NaN]*pad_len
    fig, ax1 = plt.subplots(figsize=(10,5))
    #fig, ax1 = plt.subplots()
    X = [ i + offset for i in range(len(profile1))]
    ax1.plot(X, profile1, 'b')
    #ax1.scatter(X, profile1, s=3, color='b')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel1, color='b')
    ax1.tick_params('y', colors='b')
    ax2 = ax1.twinx()
    ax2.plot(X, profile2, 'r')
    #ax2.scatter(X, profile2, s=3, color='r')
    ax2.set_ylabel(ylabel2, color='r')
    ax2.tick_params('y', colors='r')
    if xtick_loc_name:
        xtick_locs, xtick_names = xtick_loc_name
        ax2.set_xticks(xtick_locs)
        ax2.set_xticklabels(xtick_names)
    fig.tight_layout()
    plt.title(title)
    plt.savefig("pair_" + note + ".png",bbox_inches='tight')
    #plt.show()
    plt.close()


def pair_plot_mini (profile1, profile2, offset=-1000, xtick_loc_name=None, vline_locs=[], xlabel='Distance from TSS (bp)', ylabel1='Condensability (A.U.)', ylabel2="", title=""):
    assert len(profile1) == len(profile2)
    pad_len = int(len(profile1)*0.1)
    profile1[:pad_len] = [np.NaN]*pad_len
    profile1[len(profile1)-pad_len:] = [np.NaN]*pad_len
    profile2[:pad_len] = [np.NaN]*pad_len
    profile2[len(profile2)-pad_len:] = [np.NaN]*pad_len
    #fig, ax1 = plt.subplots(figsize=(10,5))
    ax1 = plt.gca()
    #fig, ax1 = plt.subplots()
    X = [ i + offset for i in range(len(profile1))]
    ax1.plot(X, profile1, 'b')
    #ax1.scatter(X, profile1, s=3, color='b')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel1, color='b')
    ax1.tick_params('y', colors='b')
    ax2 = ax1.twinx()
    ax2.plot(X, profile2, 'r', alpha=0.5)
    #ax2.scatter(X, profile2, s=3, color='r')
    ax2.set_ylabel(ylabel2, color='r')
    ax2.tick_params('y', colors='r')
    if xtick_loc_name:
        xtick_locs, xtick_names = xtick_loc_name
        ax2.set_xticks(xtick_locs)
        ax2.set_xticklabels(xtick_names)
    for vline_loc in vline_locs:
        plt.axvline(x=vline_loc, color='k', linestyle='--')
    #fig.tight_layout()
    plt.title(title)
    #plt.savefig("pair_" + note + ".png",bbox_inches='tight')
    #plt.show()
    #plt.close()


#path = "/home/spark159/../../media/spark159/sw/sp_spd_tests_detail/"
#path = "./data/"
path = ""

cell = 'H1'
sample = 'NCP'
agent = 'sp'
chr = 'chr1'

# set TSS parameters
feature = 'TSS'
profile_fname = path + '_'.join([cell, sample, agent, chr, feature]) + "_profile.txt"
#profile_fname = path + 'H1_chr1_gtf_' + feature + "_profile.txt"
#profile_fname = path + 'hg19_chr1_' + feature + "_check_profile.txt"
occprofile_fname = path + '_'.join([cell, sample, agent, chr, feature]) + "_occ_profile.txt"
#occprofile_fname = path + 'hg19_chr1_gtf_' + feature + "_occ_profile.txt"
#occprofile_fname = path + 'hg19_chr1_' + feature + "_occ_check_profile.txt"
moving_average_win = 100
offset = -1000
xtick_locs = [-500, 0, 500, 1000, 1500]
xtick_names = ['-500', 'TSS', '500', '1000', '1500']
xtick_loc_name = [xtick_locs, xtick_names]
vline_locs = [0]

## set TSS-TTS parameters
#feature = 'TSS_TTS'
#profile_fname = path + '_'.join([cell, sample, agent, chr, feature]) + "_profile.txt"
#occprofile_fname = path + '_'.join([cell, sample, agent, chr, feature]) + "_occ_profile.txt"
##profile_fname = path + 'hg19_chr1_gtf_' + feature + "_profile.txt"
##occprofile_fname = path + 'hg19_chr1_gtf_' + feature + "_occ_profile.txt"
#moving_average_win = 20
#offset = -200
#xtick_locs = [-100, 0, 600, 700]
#xtick_names = ["-2.5kb", "TSS", "TTS", "2.5kb"]
#xtick_loc_name = [xtick_locs, xtick_names]
#vline_locs = [0, 600]

# load GTF file
gID_field_values = load_file.read_GTF ("ENCFF159KBI.gtf", chr_list=[chr], mode="gene")

# load profiles
name_mean_occprofile, name_ID_occprofile = load_file.read_profile(occprofile_fname)
name_mean_profile, name_ID_profile = load_file.read_profile(profile_fname)

ID_meCpGfrac_profile = {}
for ID in (set(name_ID_profile['CNumber(CpG)'].keys()) & set(name_ID_profile['meCNumber(CpG)'].keys())):
    meGC_profile = np.asarray(statis.NN_interpolate(name_ID_profile['meCNumber(CpG)'][ID]))
    CpG_profile = np.asarray(statis.NN_interpolate(name_ID_profile['CNumber(CpG)'][ID]))
    meCpGfrac_profile =  (meGC_profile+1) / (CpG_profile+1)
    ID_meCpGfrac_profile[ID] = meCpGfrac_profile
name_ID_profile['meCpGfrac'] = ID_meCpGfrac_profile

# load RPKM
gID_RPKM = load_file.read_tsv("ENCFF174OMR.tsv")

# sort genes in the order of expression level
target1 = "work/2021_06_07_H1_sp_detail/H1-NCP-sp-8"
target2 = "ATcontent"

common_IDs = set(name_ID_profile[target1]) & set(gID_RPKM)
if target2 == "occupancy":
    common_IDs &= set(name_ID_occprofile["work/2021_06_07_H1_sp_detail/H1-NCP-sp-0"])
else:
    common_IDs &= set(name_ID_profile[target2])
common_IDs = list(common_IDs)

RPKMgID = []
for ID in common_IDs:
    RPKM = gID_RPKM[ID]
    RPKMgID.append([RPKM, ID])

RPKMgID = sorted(RPKMgID, reverse=True)

# plot gene-by-gene score pdf
note = "score2"
row = 8
col = 3
page_nums = int(math.ceil(len(common_IDs)/float(row*col))) # number of pages
page_nums = 5

pdf = matplotlib.backends.backend_pdf.PdfPages('All_genes_' + note + ".pdf")
for i in range(page_nums):
    fig = plt.figure(figsize=(15,18))
    j = 0
    while j < min(row*col, len(common_IDs)-row*col*i):
        RPKM, gID = RPKMgID[row*col*i + j]

        target_profile1 = name_ID_profile[target1][gID]
        if target2 == 'occupancy':
            target_profile2 = name_ID_occprofile["work/2021_06_07_H1_sp_detail/H1-NCP-sp-0"][gID]
        else:
            target_profile2 = name_ID_profile[target2][gID]

        # smoothing
        target_profile1 = statis.moving_average(target_profile1, moving_average_win)
        target_profile2 = statis.moving_average(target_profile2, moving_average_win)

        # plotting
        plt.subplot(row, col, j+1)
        pair_plot_mini(target_profile1, target_profile2, offset=offset, xlabel='', ylabel1='Condensability', ylabel2=target2, xtick_loc_name = xtick_loc_name, vline_locs=vline_locs, title=gID_field_values[gID]['geneName'])

        j +=1
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close()

pdf.close()

