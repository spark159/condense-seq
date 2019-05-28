import numpy as np
import matplotlib.pyplot as plt
import load_file
import statis

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


feature = 'TTS'
#profile_fname = feature + "_profile.txt"
profile_fname = feature + "_gtf_profile.txt"
moving_average_win = 100
offset = -2000
#xtick_locs = [-100, 0, 600, 700]
#xtick_names = ["-2.5kb", "TSS", "TTS", "2.5kb"]
#xtick_loc_name = [xtick_locs, xtick_names]
xtick_loc_name = None

name_mean_profile, name_ID_profile = load_file.read_profile(profile_fname)
#name_mean_occprofile, name_ID_occprofile = load_file.read_profile("data/hg19_chr1_" + feature + "_profile.txt")
name_mean_occprofile, name_ID_occprofile = load_file.read_profile(feature + "_occ_gtf_profile.txt")
mean_occprofile = statis.moving_average(name_mean_occprofile["work/condense_seq/sp1_hg19"], moving_average_win)

ID_occprofile = name_ID_occprofile["work/condense_seq/sp1_hg19"]

for name in name_mean_profile:
    name_mean_profile[name] = statis.moving_average(name_mean_profile[name], moving_average_win)

# average plot of all genes
mean_score2_profile = name_mean_profile["work/condense_seq/sp10_hg19_chr1"]
pair_plot(mean_score2_profile, mean_occprofile, offset=offset, xlabel='Distance from ' + feature +' (bp)', ylabel1='Condensability (A.U.)', ylabel2="Occupancy", xtick_loc_name = xtick_loc_name, note="Occ" + '_' + feature)

for name in sorted(name_mean_profile):
    if name == "work/condense_seq/sp10_hg19_chr1":
        continue
    #if name == "meGCNumber" or name == "CpGNumber":
    #    continue
    mean_profile = name_mean_profile[name]
    pair_plot(mean_score2_profile, mean_profile, offset=offset, xlabel='Distance from ' + feature +' (bp)', ylabel1='Condensability (A.U.)', ylabel2=name,  xtick_loc_name = xtick_loc_name, note=name.split('/')[-1] + '_' + feature)

mean_meCpGfrac_profile = name_mean_profile["meGCNumber"] / (2*name_mean_profile["CpGNumber"])
pair_plot(mean_score2_profile, mean_meCpGfrac_profile, offset=offset, xlabel='Distance from ' + feature +' (bp)', ylabel1='Condensability (A.U.)', ylabel2="meCpG frac",  xtick_loc_name = xtick_loc_name, note="meCpGfrac" + '_' + feature)


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
    plt.savefig("RPKM_quantile_profile_" + feature + "_" + name.split('/')[-1] + ".png",bbox_inches='tight')
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
    

