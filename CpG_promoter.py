import sys
import Interval_dict
import numpy as np
import matplotlib.pyplot as plt
import load_file
import statis

def read_CpGTable (fname, chr_choice=None):
    id = 0
    id_field_values = {}
    for line in open(fname):
        cols = line.strip().split('\t')
        _, chr, st, ed, name, length, cpgNum, gcNum, perCpg, perGC, obsExp = cols
        if chr_choice and chr not in chr_choice:
            continue
        #num = int(num.strip().split(":_")[1])
        interval = (int(st), int(ed))
        id_field_values[id] = {}
        id_field_values[id]["num"] = int(cpgNum)
        id_field_values[id]["interval"] = interval
        id +=1
    return id_field_values

path = './data/'

# find CpG promoter
gID_field_values, field_gID_values = load_file.read_GTF (path+"Homo_sapiens.GRCh37.87.gtf", chr_list=["chr1"], mode="both")

gID_ginterval = {}
for gID in gID_field_values:
    TSS = gID_field_values[gID]['TSS']
    TTS = gID_field_values[gID]['TTS']
    strand = gID_field_values[gID]['strand']
    #interval = (TSS-500, TSS+500)
    if strand == '+':
        #interval = (TSS, TTS)
        interval = (TSS-1000, TSS+500)
        #interval = (TSS, TSS+2500)
    else:
        #interval = (TTS, TSS)
        interval = (TSS-500, TSS+1000)
        #interval = (TSS-2500, TSS)
    gID_ginterval[gID] = interval

ginterval_dict = Interval_dict.double_hash(gID_ginterval, 10000, 250000000)
CpGID_field_values = read_CpGTable(path+"cpgIslandExt.txt", chr_choice=["chr1"])

CpG_genes = []
for CpGID in CpGID_field_values:
    st, ed = CpGID_field_values[CpGID]["interval"]
    IDs = ginterval_dict.find_range(st, ed)
    CpG_genes += IDs

CpG_genes = set(CpG_genes)
NonCpG_genes = set(gID_ginterval.keys()) - CpG_genes
CpG_genes, NonCpG_genes = list(CpG_genes), list(NonCpG_genes)

print "completed to sort CpG genes"
print "CpG genes: ", len(CpG_genes)
print "Non CpG genes: ", len(NonCpG_genes)

#sys.exit(1)

# plot profile
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
    plt.show()
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

# mean TSS plot
feature = 'TSS'
ID_choice = NonCpG_genes
profile_fname = path + 'hg19_chr1_gtf_' + feature + "_profile.txt"
#profile_fname = path + 'hg19_chr1_' + feature + "_check_profile.txt"
occprofile_fname = path + 'hg19_chr1_gtf_' + feature + "_occ_profile.txt"
#occprofile_fname = path + 'hg19_chr1_' + feature + "_occ_check_profile.txt"
moving_average_win = 100
offset = -1000
xtick_loc_name = None

name_mean_profile, name_ID_profile = load_file.read_profile(profile_fname, ID_choice=ID_choice)
for name in name_mean_profile:
    name_mean_profile[name] = statis.moving_average(name_mean_profile[name], moving_average_win)

name_mean_occprofile, name_ID_occprofile = load_file.read_profile(occprofile_fname, ID_choice=ID_choice)
mean_occprofile = statis.moving_average(name_mean_occprofile["data/sp_spd_tests_detail/sp1"], moving_average_win)
mean_occprofile2 = statis.moving_average(name_mean_occprofile["data/sp_spd_tests_detail/sp7"], moving_average_win)
mean_occprofile3 = statis.moving_average(name_mean_occprofile["data/sp_spd_tests_detail/sp8"], moving_average_win)
#mean_occprofile = name_mean_occprofile["data/sp_spd_tests_detail/sp1"]
ID_occprofile = name_ID_occprofile["data/sp_spd_tests_detail/sp1"]

occprofiles = [mean_occprofile/sum(mean_occprofile),mean_occprofile2/sum(mean_occprofile2), mean_occprofile3/sum(mean_occprofile3)]
single_plot(occprofiles, names=['Input', 'Sp6', 'Sp7'])

# average plot of all genes
mean_score1_profile = name_mean_profile["data/sp_spd_tests_detail/sp7"]
pair_plot(mean_score1_profile, mean_occprofile, offset=offset, xlabel='Distance from ' + feature +' (bp)', ylabel1='Condensability (A.U.)', ylabel2="Occupancy", xtick_loc_name = xtick_loc_name, note="Occ" + '_' + feature)

for name in sorted(name_mean_profile):
    if name in  ["data/sp_spd_tests_detail/sp7", "data/sp_spd_tests_detail/sp8"]:
        continue
    mean_profile = name_mean_profile[name]
    pair_plot(mean_score1_profile, mean_profile, offset=offset, xlabel='Distance from ' + feature +' (bp)', ylabel1='Condensability (A.U.)', ylabel2=name,  xtick_loc_name = xtick_loc_name, note=name.split('/')[-1] + '_' + feature)

mean_meCpGfrac_profile = name_mean_profile["meGCNumber"] / (2*name_mean_profile["CpGNumber"])
pair_plot(mean_score1_profile, mean_meCpGfrac_profile, offset=offset, xlabel='Distance from ' + feature +' (bp)', ylabel1='Condensability (A.U.)', ylabel2="meCpG frac",  xtick_loc_name = xtick_loc_name, note="meCpGfrac" + '_' + feature)









"""
def pair_plot (profile1, profile2, offset=-1000, xtick_loc_name=None, xlabel='Distance from TSS (bp)', ylabel1='Condensability (A.U.)', ylabel2="", note=""):
    assert len(profile1) == len(profile2)
    profile1[:100] = [np.NaN]*100
    profile1[len(profile1)-100:] = [np.NaN]*100
    profile2[:100] = [np.NaN]*100
    profile2[len(profile2)-100:] = [np.NaN]*100
    fig, ax1 = plt.subplots(figsize=(10,5))
    #fig, ax1 = plt.subplots()
    X = [ i + offset for i in range(len(profile1))]
    ax1.plot(X, profile1, 'b')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel1, color='b')
    ax1.tick_params('y', colors='b')
    ax2 = ax1.twinx()
    ax2.plot(X, profile2, 'r')
    ax2.set_ylabel(ylabel2, color='r')
    ax2.tick_params('y', colors='r')
    if xtick_loc_name:
        xtick_locs, xtick_names = xtick_loc_name
        ax2.set_xticks(xtick_locs)
        ax2.set_xticklabels(xtick_names)
    fig.tight_layout()
    plt.savefig("pair_" + note + ".png",bbox_inches='tight')
    #plt.show()
    plt.close()

# plot profile
feature = 'TSS'
#profile_fname = feature + "_profile.txt"
#profile_fname = feature + "_gtf_profile.txt"
profile_fname = "TSS_GTF_newmetric_profile.txt"
moving_average_win = 100
offset = -1000
#xtick_locs = [-100, 0, 600, 700]
#xtick_names = ["-2.5kb", "TSS", "TTS", "2.5kb"]
#xtick_loc_name = [xtick_locs, xtick_names]
xtick_loc_name = None
ID_choice = CpG_genes
#ID_choice = NonCpG_genes
note = "_CpGgenes_newmetric"
#note = "_NonCpGgenes_newmetric"
#note = "_CpGgenes"
#note = "_NonCpGgenes"


name_mean_profile, name_ID_profile = load_file.read_profile(profile_fname, ID_choice=ID_choice)
#name_mean_occprofile, name_ID_occprofile = load_file.read_profile("data/hg19_chr1_" + feature + "_profile.txt")
name_mean_occprofile, name_ID_occprofile = load_file.read_profile(feature + "_occ_gtf_profile.txt", ID_choice=ID_choice)
mean_occprofile = statis.moving_average(name_mean_occprofile["work/condense_seq/sp1_hg19"], moving_average_win)

ID_occprofile = name_ID_occprofile["work/condense_seq/sp1_hg19"]

for name in name_mean_profile:
    name_mean_profile[name] = statis.moving_average(name_mean_profile[name], moving_average_win)

# average plot of all genes
mean_score2_profile = name_mean_profile["work/condense_seq/sp10_hg19_chr1"]
pair_plot(mean_score2_profile, mean_occprofile, offset=offset, xlabel='Distance from ' + feature +' (bp)', ylabel1='Condensability (A.U.)', ylabel2="Occupancy", xtick_loc_name = xtick_loc_name, note="Occ" + '_' + feature + note)

for name in sorted(name_mean_profile):
    if name == "work/condense_seq/sp10_hg19_chr1":
        continue
    #if name == "meGCNumber" or name == "CpGNumber":
    #    continue
    mean_profile = name_mean_profile[name]
    pair_plot(mean_score2_profile, mean_profile, offset=offset, xlabel='Distance from ' + feature +' (bp)', ylabel1='Condensability (A.U.)', ylabel2=name,  xtick_loc_name = xtick_loc_name, note=name.split('/')[-1] + '_' + feature + note)

mean_meCpGfrac_profile = name_mean_profile["meGCNumber"] / (2*name_mean_profile["CpGNumber"])
pair_plot(mean_score2_profile, mean_meCpGfrac_profile, offset=offset, xlabel='Distance from ' + feature +' (bp)', ylabel1='Condensability (A.U.)', ylabel2="meCpG frac",  xtick_loc_name = xtick_loc_name, note="meCpGfrac" + '_' + feature + note)



# paritions according to gene expression
gID_RPKM = load_file.read_RPKM ("data/GSE63124_all_gene_raw_readcounts.txt", "data/Homo_sapiens.GRCh37.87.gtf", "chr1")
ID_profile = name_ID_profile["work/condense_seq/sp10_hg19_chr1"]
name_ID_profile['Occ'] = name_ID_occprofile["work/condense_seq/sp1_hg19"]

for name in name_ID_profile:
    ID_profile = name_ID_profile[name]
    IDs = list(set(ID_profile) & set(gID_RPKM))

    gID_score = {}
    scoregID = []
    for ID in IDs:
        score = gID_RPKM[ID]
        gID_score[ID] = score
        scoregID.append([score, ID])

    frac = [1 for i in range(5)]
    groups = statis.quantile(gID_score, 5, frac=frac)

    fig = plt.figure()
    for i in range(len(groups)):
        IDs = groups[i]
        #profile = [ ID_occprofile[ID] for ID in IDs ]
        profile = [ ID_profile[ID] for ID in IDs ]
        profile = np.nanmean(profile, axis=0)
        profile = statis.moving_average(profile, moving_average_win)
        profile[:100] = [np.NaN]*100
        profile[len(profile)-100:] = [np.NaN]*100
        profile[:100] = [np.NaN]*100
        profile[len(profile)-100:] = [np.NaN]*100
        X = [ k + offset for k in range(len(profile))]
        plt.plot(X, profile, label='RPKM quantile '+str(i+1), alpha=0.8)
    if xtick_loc_name:
        xtick_locs, xtick_names = xtick_loc_name
        plt.xticks(xtick_locs, xtick_names)
    plt.xlabel('Distance from ' + feature +' (bp)')
    plt.title(name)
    plt.legend()
    plt.savefig("RPKM_quantile_profile_" + feature + "_" + name.split('/')[-1] + note + ".png",bbox_inches='tight')
    #plt.show()
    plt.close()
"""
