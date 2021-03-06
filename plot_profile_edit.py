import sys
import numpy as np
import matplotlib.pyplot as plt
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

#path = "/home/spark159/../../media/spark159/sw/sp_spd_tests_detail/"
path = "./data/"

"""
# mean TSS plot
feature = 'TSS'
profile_fname = path + 'hg19_chr1_gtf_' + feature + "_profile.txt"
#profile_fname = path + 'hg19_chr1_' + feature + "_check_profile.txt"
occprofile_fname = path + 'hg19_chr1_gtf_' + feature + "_occ_profile.txt"
#occprofile_fname = path + 'hg19_chr1_' + feature + "_occ_check_profile.txt"
moving_average_win = 100
offset = -1000
xtick_loc_name = None

name_mean_profile, name_ID_profile = load_file.read_profile(profile_fname)
for name in name_mean_profile:
    name_mean_profile[name] = statis.moving_average(name_mean_profile[name], moving_average_win)

name_mean_occprofile, name_ID_occprofile = load_file.read_profile(occprofile_fname)
mean_occprofile = statis.moving_average(name_mean_occprofile["data/sp_spd_tests_detail/sp1"], moving_average_win)
mean_occprofile2 = statis.moving_average(name_mean_occprofile["data/sp_spd_tests_detail/sp7"], moving_average_win)
mean_occprofile3 = statis.moving_average(name_mean_occprofile["data/sp_spd_tests_detail/sp8"], moving_average_win)
#mean_occprofile = name_mean_occprofile["data/sp_spd_tests_detail/sp1"]
ID_occprofile = name_ID_occprofile["data/sp_spd_tests_detail/sp1"]

occprofiles = [mean_occprofile/sum(mean_occprofile),mean_occprofile2/sum(mean_occprofile2), mean_occprofile3/sum(mean_occprofile3)]
#single_plot(occprofiles, names=['Input', 'Sp6', 'Sp7'])

# average plot of all genes
mean_score1_profile = name_mean_profile["data/sp_spd_tests_detail/sp7"]
#pair_plot(mean_score1_profile, mean_occprofile, offset=offset, xlabel='Distance from ' + feature +' (bp)', ylabel1='Condensability (A.U.)', ylabel2="Occupancy", xtick_loc_name = xtick_loc_name, note="Occ" + '_' + feature)

for name in sorted(name_mean_profile):
    if name in  ["data/sp_spd_tests_detail/sp7", "data/sp_spd_tests_detail/sp8"]:
        continue
    mean_profile = name_mean_profile[name]
    pair_plot(mean_score1_profile, mean_profile, offset=offset, xlabel='Distance from ' + feature +' (bp)', ylabel1='Condensability (A.U.)', ylabel2=name,  xtick_loc_name = xtick_loc_name, note=name.split('/')[-1] + '_' + feature)

mean_meCpGfrac_profile = name_mean_profile["meGCNumber"] / (2*name_mean_profile["CpGNumber"])
#pair_plot(mean_score1_profile, mean_meCpGfrac_profile, offset=offset, xlabel='Distance from ' + feature +' (bp)', ylabel1='Condensability (A.U.)', ylabel2="meCpG frac",  xtick_loc_name = xtick_loc_name, note="meCpGfrac" + '_' + feature)


# paritions according to gene expression
gID_RPKM = load_file.read_RPKM (path+"GSE63124_all_gene_raw_readcounts.txt", path+"Homo_sapiens.GRCh37.87.gtf", "chr1")
ID_profile = name_ID_profile["data/sp_spd_tests_detail/sp7"]
name_ID_profile['Occ'] = name_ID_occprofile["data/sp_spd_tests_detail/sp1"]


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
        X = [ k + offset for k in range(len(profile))]
        plt.plot(X, profile, label='RPKM quantile '+str(i+1), alpha=0.8)
    if xtick_loc_name:
        xtick_locs, xtick_names = xtick_loc_name
        plt.xticks(xtick_locs, xtick_names)
    plt.xlabel('Distance from ' + feature +' (bp)')
    plt.title(name)
    plt.legend()
    #plt.savefig("RPKM_quantile_profile_" + feature + "_" + name.split('/')[-1] + ".png",bbox_inches='tight')
    #plt.show()
    plt.close()


# heatmap plot by gene expression level
name_list = ['data/sp_spd_tests_detail/sp7', 'k27ac', 'k4me3', 'k9ac', 'Occ']
#name_list = ['data/sp_spd_tests_detail/sp7']

common_IDs = set(gID_RPKM)
for name in name_list:
    ID_profile = name_ID_profile[name]
    common_IDs = common_IDs & set(ID_profile)

RPKMgID = []
for ID in common_IDs:
    RPKM = gID_RPKM[ID]
    RPKMgID.append([RPKM, ID])

RPKMgID = sorted(RPKMgID)

name_img = {}
for name in name_list:
    ID_profile = name_ID_profile[name]
    img = []
    # smoothing
    for RPKM, ID in RPKMgID:
        profile = ID_profile[ID]
        profile = statis.moving_average(profile, moving_average_win)
        profile[:100] = [profile[100]]*100
        profile[len(profile)-100:] = [profile[len(profile)-101]]*100
        img.append(profile)
    img = np.asarray(img)
    # rescaling
    if not name.split('/')[-1].startswith('sp'):
        mini = np.min(img)
        img = np.log2(img - mini + 1)
    else:
        img = 2*img
        
    name_img[name] = img

cmap_list = ['rainbow', 'YlOrRd', 'YlGn', 'Purples', 'Greys', 'Blues', 'Oranges']
vlim_list = [ [-4, 4], [None, None], [None, None], [None, None], [None, None], [None, None]]
xtick_labels = [ str(k + offset) for k in range(0, len(profile), 1000)]
for i in range(len(name_list)):
    name = name_list[i]
    img = name_img[name]
    #fig = plt.figure(figsize=(5,10))
    fig = plt.figure(figsize=(6.5,10))
    plt.imshow(img, aspect='auto', cmap=cmap_list[i], vmin=vlim_list[i][0], vmax=vlim_list[i][1])
    #plt.imshow(img, aspect='auto', cmap='jet')
    #plt.title(name)
    plt.xticks(range(0, len(profile), 1000), xtick_labels)
    plt.savefig("heatmap_profile_byRPKM_" + feature + "_" + name.split('/')[-1] + ".png",bbox_inches='tight')
    #plt.show()
    plt.close()



#scoregID = sorted(scoregID, reverse=True)
#img = []
#for score, ID in scoregID:
#    #profile = ID_occprofile[ID]
#    profile = ID_profile[ID]
#    img.append(statis.moving_average(profile, 100))

#fig = plt.figure()
#plt.imshow(img, aspect='auto')
#plt.show()
#plt.close()

#order = []
#for ID in ID_occprofile:
#    profile = statis.moving_average(ID_occprofile[ID], 100)
#    value = sum(profile)
#    order.append([value, ID])

#order = sorted(order, reverse=True)

#img = []
#for value, ID in order:
#    profile = statis.moving_average(ID_occprofile[ID], 100)
#    fig = plt.figure()
#    plt.plot(profile)
#    plt.show()
#    plt.close()
#    #img.append(profile[100:len(profile)-100])



# mean TTS plot
feature = 'TTS'
profile_fname = path + 'hg19_chr1_gtf_' + feature + "_profile.txt"
occprofile_fname = path + 'hg19_chr1_gtf_' + feature + "_occ_profile.txt"
moving_average_win = 100
offset = -2000
xtick_loc_name = None

name_mean_profile, name_ID_profile = load_file.read_profile(profile_fname)
for name in name_mean_profile:
    name_mean_profile[name] = statis.moving_average(name_mean_profile[name], moving_average_win)

name_mean_occprofile, name_ID_occprofile = load_file.read_profile(occprofile_fname)
mean_occprofile = statis.moving_average(name_mean_occprofile["data/sp_spd_tests_detail/sp1"], moving_average_win)
#mean_occprofile = name_mean_occprofile["data/sp_spd_tests_detail/sp1"]
ID_occprofile = name_ID_occprofile["data/sp_spd_tests_detail/sp1"]

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

# paritions according to gene expression
gID_RPKM = load_file.read_RPKM ("/home/spark159/../../media/spark159/sw/dataforcondense/GSE63124_all_gene_raw_readcounts.txt", "/home/spark159/../../media/spark159/sw/dataforcondense/Homo_sapiens.GRCh37.87.gtf", "chr1")
ID_profile = name_ID_profile["data/sp_spd_tests_detail/sp7"]
name_ID_profile['Occ'] = name_ID_occprofile["data/sp_spd_tests_detail/sp1"]

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
        X = [ k + offset for k in range(len(profile))]
        plt.plot(X, profile, label='RPKM quantile '+str(i+1), alpha=0.8)
    if xtick_loc_name:
        xtick_locs, xtick_names = xtick_loc_name
        plt.xticks(xtick_locs, xtick_names)
    plt.xlabel('Distance from ' + feature +' (bp)')
    plt.title(name)
    plt.legend()
    plt.savefig("RPKM_quantile_profile_" + feature + "_" + name.split('/')[-1] + ".png",bbox_inches='tight')
    #plt.show()
    plt.close()
"""

# mean TSS-TTS plot
feature = 'TSS_TTS'
profile_fname = path + 'hg19_chr1_gtf_' + feature + "_profile.txt"
occprofile_fname = path + 'hg19_chr1_gtf_' + feature + "_occ_profile.txt"
moving_average_win = 20
offset = -200
xtick_locs = [-100, 0, 600, 700]
xtick_names = ["-2.5kb", "TSS", "TTS", "2.5kb"]
xtick_loc_name = [xtick_locs, xtick_names]

name_mean_profile, name_ID_profile = load_file.read_profile(profile_fname)
for name in name_mean_profile:
    name_mean_profile[name] = statis.moving_average(name_mean_profile[name], moving_average_win)

name_mean_occprofile, name_ID_occprofile = load_file.read_profile(occprofile_fname)
#mean_occprofile = name_mean_occprofile["data/sp_spd_tests_detail/sp1"]
mean_occprofile = statis.moving_average(name_mean_occprofile["data/sp_spd_tests_detail/sp1"], moving_average_win)
ID_occprofile = name_ID_occprofile["data/sp_spd_tests_detail/sp1"]

# average plot of all genes
mean_score1_profile = name_mean_profile["data/sp_spd_tests_detail/sp7"]
#pair_plot(mean_score1_profile, mean_occprofile, offset=offset, xlabel='Distance from ' + feature +' (bp)', ylabel1='Condensability (A.U.)', ylabel2="Occupancy", xtick_loc_name = xtick_loc_name, note="Occ" + '_' + feature)

for name in sorted(name_mean_profile):
    if name in  ["data/sp_spd_tests_detail/sp7", "data/sp_spd_tests_detail/sp8"]:
        continue
    mean_profile = name_mean_profile[name]
    #pair_plot(mean_score1_profile, mean_profile, offset=offset, xlabel='Distance from ' + feature +' (bp)', ylabel1='Condensability (A.U.)', ylabel2=name,  xtick_loc_name = xtick_loc_name, note=name.split('/')[-1] + '_' + feature)

mean_meCpGfrac_profile = name_mean_profile["meGCNumber"] / (2*name_mean_profile["CpGNumber"])
#pair_plot(mean_score1_profile, mean_meCpGfrac_profile, offset=offset, xlabel='Distance from ' + feature +' (bp)', ylabel1='Condensability (A.U.)', ylabel2="meCpG frac",  xtick_loc_name = xtick_loc_name, note="meCpGfrac" + '_' + feature)

# paritions according to gene expression
gID_RPKM = load_file.read_RPKM (path+"GSE63124_all_gene_raw_readcounts.txt", path+"Homo_sapiens.GRCh37.87.gtf", "chr1")
ID_profile = name_ID_profile["data/sp_spd_tests_detail/sp7"]
name_ID_profile['Occ'] = name_ID_occprofile["data/sp_spd_tests_detail/sp1"]

ID_meCpGfrac_profile = {}
for ID in (set(name_ID_profile["meGCNumber"].keys()) & set(name_ID_profile["CpGNumber"].keys())):
    #meGC_profile = statis.moving_average(name_ID_profile["meGCNumber"][ID], moving_average_win)
    #CpG_profile = statis.moving_average(name_ID_profile["CpGNumber"][ID], moving_average_win)
    meGC_profile = np.asarray(statis.NN_interpolate(name_ID_profile["meGCNumber"][ID]))
    CpG_profile = np.asarray(statis.NN_interpolate(name_ID_profile["CpGNumber"][ID]))
    meCpGfrac_profile =  (meGC_profile+1) / (2*CpG_profile+1)
    ID_meCpGfrac_profile[ID] = meCpGfrac_profile
name_ID_profile['meCpGfrac'] = ID_meCpGfrac_profile

"""
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
        X = [ k + offset for k in range(len(profile))]
        plt.plot(X, profile, label='RPKM quantile '+str(i+1), alpha=0.8)
    if xtick_loc_name:
        xtick_locs, xtick_names = xtick_loc_name
        plt.xticks(xtick_locs, xtick_names)
    #plt.xlabel('Distance from ' + feature +' (bp)')
    plt.title(name)
    plt.legend()
    plt.savefig("RPKM_quantile_profile_" + feature + "_" + name.split('/')[-1] + ".png",bbox_inches='tight')
    #plt.show()
    plt.close()
"""


# heatmap plot by gene expression level
#name_list = ['data/sp_spd_tests_detail/sp7', 'ATcontent', 'k36me3', 'k9me2', 'k27ac', 'k4me3', 'k9ac', 'Occ']
#name_list = ['data/sp_spd_tests_detail/sp7']
name_list = ['data/sp_spd_tests_detail/sp7', 'meCpGfrac', 'k36me3', 'ATcontent']
#name_list = ['data/sp_spd_tests_detail/sp7', 'ATcontent', 'k36me3', 'k9me2']

common_IDs = set(gID_RPKM)
for name in name_list:
    #if name == 'meCpGfrac':
    #    continue
    ID_profile = name_ID_profile[name]
    common_IDs = common_IDs & set(ID_profile)

#if 'meCpGfrac' in name_list:
#    common_IDs = common_IDs & set(name_ID_profile['CpGNumber'])
#    common_IDs = common_IDs & set(name_ID_profile['meGCNumber'])

RPKMgID = []
for ID in common_IDs:
    RPKM = gID_RPKM[ID]
    RPKMgID.append([RPKM, ID])

RPKMgID = sorted(RPKMgID)

name_img = {}
for name in name_list:
    ID_profile = name_ID_profile[name]
    img = []
    # smoothing
    for RPKM, ID in RPKMgID:
        profile = ID_profile[ID]
        profile = statis.moving_average(profile, moving_average_win)
        profile[:50] = [profile[50]]*50
        profile[len(profile)-50:] = [profile[len(profile)-51]]*50
        img.append(profile)
    img = np.asarray(img)
    # rescaling
    if not name.split('/')[-1].startswith('sp') or name in ['meCpGfrac', 'ATcontent']:
        mini = np.min(img)
        img = np.log2(img - mini + 1)
    else:
        img = 2*img
        
    name_img[name] = img

#if 'meCpGfrac' in name_list:
#    ID_meCpGfrac_profile = {}
#    for ID in common_IDs:
#        meGC_profile = statis.moving_average(name_ID_profile["meGCNumber"][ID], moving_average_win)
#        CpG_profile = statis.moving_average(name_ID_profile["CpGNumber"][ID], moving_average_win)
#        meCpGfrac_profile =  meGC_profile / (2*CpG_profile+1)
#        ID_meCpGfrac_profile[ID] = meCpGfrac_profile
#    img = []
#    for RPKM, ID in RPKMgID:
#        profile = ID_meCpGfrac_profile[ID]
#        profile[:100] = [profile[100]]*100
#        profile[len(profile)-100:] = [profile[len(profile)-101]]*100
#        img.append(profile)
#    img = np.asarray(img)
#    mini = np.min(img)
#    img = np.log2(img - mini + 1)
#    name_img['meCpGfrac'] = img

cmap_list = ['rainbow', 'cool', 'viridis', 'hot']
vlim_list = [ [-3, 3], [None, None], [None, None], [None, None], [None, None], [None, None]]
for i in range(len(name_list)):
    name = name_list[i]
    img = name_img[name]
    #fig = plt.figure(figsize=(5,10))
    fig = plt.figure(figsize=(6.5,10))
    plt.imshow(img, aspect='auto', cmap=cmap_list[i], vmin=vlim_list[i][0], vmax=vlim_list[i][1])
    #plt.imshow(img, aspect='auto', cmap='jet')
    #plt.title(name)
    xtick_locs, xtick_names = xtick_loc_name
    plt.xticks(np.asarray(xtick_locs) - offset, xtick_names)
    #plt.xticks(range(0, len(profile), 1000), xtick_labels)
    plt.savefig("heatmap_profile_byRPKM_" + feature + "_" + name.split('/')[-1] + ".png",bbox_inches='tight')
    #plt.show()
    plt.close()
