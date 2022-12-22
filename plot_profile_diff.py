import sys
import copy
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
        #ax.set_xticklabels(xtick_names)
        ax.set_xticklabels(xtick_names, rotation=45)
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


def multi_plot (profile_list,
                offset=-1000,
                xtick_loc_name=None,
                xlabel='Distance from TSS (bp)',
                ylabel='Condensability (A.U.)',
                label_list=None, 
                title="",
                note=""):

    for i in range(len(profile_list)-1):
        for j in range(i+1, len(profile_list)):
            assert len(profile_list[i]) == len(profile_list[j])

    if label_list == None:
        label_list = [None]*len(profile_list)


    fig = plt.figure(figsize=(10,5))
    
    for i in range(len(profile_list)):
        profile = profile_list[i]
        newprofile = copy.deepcopy(profile)
        pad_len = int(len(profile)*0.1)
        newprofile[:pad_len] = [np.NaN]*pad_len
        newprofile[len(profile)-pad_len:] = [np.NaN]*pad_len

        X = [k + offset for k in range(len(newprofile))]
        plt.plot(X, newprofile, label=label_list[i])

    if xtick_loc_name:
        xtick_locs, xtick_names = xtick_loc_name
        plt.xticks(xtick_locs, xtick_names)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.savefig("multi_" + note + ".png",bbox_inches='tight')
    #plt.show()
    plt.close()



### parameters
#path = "/home/spark159/../../media/spark159/sw/sp_spd_tests_detail/"
#path = "./data/"
#path = ""
path = "/home/spark159/../../media/spark159/sw/"


### experiment information
#cell1, cell2 = "H1", "GM"
cell1, cell2, cell3 = "mCD8T", "mCD8T", "mCD8T"
#sample1, sample2 = "NCP", "NCP"
sample1, sample2, sample3 = "WT-NCP", "inht-NCP", "KO-NCP"
agent = "sp"

exp_list = [(cell1, sample1, agent),
            (cell2, sample2, agent),
            (cell3, sample3, agent)]

### select regions
chr_list = ['chr1']
#chr_list = ['chr' + str(i) for i in range(1, 20)] #mouse
#chr_list += ['chrX']

domain = 'TSS-TTS'
#domain = 'TSS'


# RNA-seq data fname
rnafname = "GSE136898_rawCounts.txt"
gtfname = "gencodeM21pri-UCSC-tRNAs-ERCC-phiX.gtf"
ID_field_values = load_file.read_GTF ("gencodeM21pri-UCSC-tRNAs-ERCC-phiX.gtf",
                                      mode="gene",
                                      chr_list=chr_list) #mouse    


# promoter classification file
pcfname = '41586_2007_BFnature06008_MOESM456_ESM.csv'
loc_field_value = load_file.read_tabular_file(pcfname)
gname_pclass = {}
for loc in loc_field_value:
    gname = loc_field_value[loc]['Gene(s)']
    CpG_class = loc_field_value[loc]['Class']
    ptm_class = loc_field_value[loc]['ESC']
    pclass = (CpG_class, ptm_class)
    gname_pclass[gname] = pclass

# rank file (dzscore around TSS 5kb window)
rnkfname1 = "inht-WT_gID.rnk"
rnkfname2 = "KO-WT_gID.rnk"
ID_rnkvalue1 = load_file.read_tabular_file(rnkfname1, header=False, mode='col')[0]
ID_rnkvalue2 = load_file.read_tabular_file(rnkfname2, header=False, mode='col')[0]


# plot parameters
if domain == 'TSS':
    moving_average_win = 20
    profile_len = 1000
    up_win = 2500
    down_win = 5000
    pad_len = int(profile_len*0.1)
    left_len = int(round(up_win*float(profile_len)/(up_win+down_win+1)))
    right_len = profile_len - left_len - 1
    offset = -left_len
    xtick_locs = [-left_len + pad_len,
                  0,
                  right_len - pad_len]
    xtick_names = ["-2.5kb", "TSS", "5kb"]
    xtick_loc_name = [xtick_locs, xtick_names]


elif domain == 'TSS-TTS':
    moving_average_win = 20
    #offset = -200
    #xtick_locs = [-100, 0, 600, 700]
    #xtick_names = ["-2.5kb", "TSS", "TTS", "2.5kb"]
    profile_len = 1000
    domain_frac = 0.5
    margin_frac = 1.0 - domain_frac
    up_win = 5000
    down_win = 2500
    left_len = int((margin_frac*profile_len)*(float(up_win)/(up_win + down_win)))
    right_len = int((margin_frac*profile_len)*(float(down_win)/(up_win + down_win)))
    pad_len = int(profile_len*0.1)
    offset = -left_len
    xtick_locs = [-left_len+pad_len,
                  0,
                  int(domain_frac*profile_len),
                  int(domain_frac*profile_len) + right_len - pad_len ]
    xtick_names = ["-5kb", "TSS", "TTS", "2.5kb"]
    xtick_loc_name = [xtick_locs, xtick_names]
    
elif domain == 'TTS':
    moving_average_win = 100
    offset = -2000
    xtick_loc_name = None
    name_list = ["mCD8T-KO-NCP-sp-8", 'H3K27ac', 'H3K4me3', 'H3K4me1', 'Occ']
    cmap_list = ['rainbow', 'YlOrRd', 'YlGn', 'Purples', 'Greys', 'Blues', 'Oranges']
    vlim_list = [[None, None], [None, None], [None, None], [None, None], [None, None], [None, None]]


# field names
target_names = ['mCD8T-WT-NCP-sp-8', 'mCD8T-inht-NCP-sp-8', 'mCD8T-KO-NCP-sp-8']
feature_names = ['H3K27me3', 'H3K4me3', 'H3K27ac', 'H3K36me3', 'H3K4me1', 'H3K9me3', 'ATcontent']
names = target_names + feature_names


# labels for replacing field name
#labels = ['WT', '+inht', 'ODC KO']
#labels = ['WT', '+inht', 'ODC KO', 'H3K27me3', 'H3K4me3']
labels = ['WT', '+inht', 'ODC KO'] + feature_names

# load occupancy choice
load_occ = False
    
name_ID_profile = {}
# read profile files
for cell, sample, agent in exp_list:
    print "loading %s-%s-%s" % (cell, sample, agent)

    for chr_name in chr_list:
        print "\t reading %s" % (chr_name)
        
        # load profile files
        #profile_fname = path + '_'.join([cell, sample, agent, chr_name, domain]) + "_profile.txt"
        profile_fname = path + '_'.join([cell, sample, agent, chr_name, 'zscore', domain]) + "_profile.txt"
        name_ID_data = load_file.read_profile(profile_fname,
                                              name_choice=names,
                                              average = False)

        for name in name_ID_data:
            newname = name.rsplit('/', 1)[-1].rsplit('.', 1)[0]
            if newname not in name_ID_profile:
                name_ID_profile[newname] = {}
            name_ID_profile[newname].update(name_ID_data[name])

        # load occ files
        if load_occ:
            occ_fname = path + '_'.join([cell, sample, agent, chr_name, 'occ', domain]) + "_profile.txt"
            name_ID_data = load_file.read_profile(occ_fname,
                                                  name_choice=names,
                                                  average=False)

            for name in name_ID_data:
                newname = name.rsplit('/', 1)[-1].rsplit('.', 1)[0] + '_occ'
                if newname not in name_ID_profile:
                    name_ID_profile[newname] = {}
                name_ID_profile[newname].update(name_ID_data[name])

        del name_ID_data
    

# get mean profile for all genes
name_mprofile = {}
name_smprofile = {}
if False:
    print "data average and smoothing"
    for name in names:
        name_mprofile[name] = np.nanmean(name_ID_profile[name].values(), axis=0)
        name_smprofile[name] = statis.moving_average(name_mprofile[name], moving_average_win)


# plot mean profiles
if False:
    profile_list, label_list = [], []
    for name in names:
        profile_list.append(name_smprofile[name])
        label_list.append(name)

    multi_plot (profile_list,
                offset=offset,
                xtick_loc_name=xtick_loc_name,
                xlabel='Distance from ' + domain +' (bp)',
                ylabel='Condensability (A.U.)',
                label_list=label_list, 
                title="",
                note=domain)

# partition quantiles according to gene expression level
q_IDs = None
if False:
    ID_RPKM = load_file.read_RPKM_new (rnafname, gtfname, chr_list=chr_list)

    common_IDs = set(ID_RPKM.keys())
    for name in names:
        ID_profile = name_ID_profile[name]
        common_IDs &= set(ID_profile)
    
    common_IDs = list(common_IDs)
    q_IDs = statis.quantile({ID:ID_RPKM[ID] for ID in common_IDs}, 5)


# plot quantiles for each features
if False:
    for name in names:
        ID_profile = name_ID_profile[name]

        fig = plt.figure(figsize=(5,3.5))
        for i in range(len(q_IDs)):
            profile = np.nanmean([ID_profile[ID] for ID in q_IDs[i]], axis=0)
            profile = statis.moving_average(profile, moving_average_win)
            profile[:100] = [np.NaN]*100
            profile[len(profile)-100:] = [np.NaN]*100
            X = [ k + offset for k in range(len(profile))]
            #plt.plot(X, profile, label='RPKM quantile '+str(i+1), alpha=0.8)
            plt.plot(X, profile, label='FPKM quantile '+str(i+1), alpha=0.8, lw=4)
        if xtick_loc_name:
            xtick_locs, xtick_names = xtick_loc_name
            plt.xticks(xtick_locs, xtick_names, fontsize=15, rotation=45)
        plt.xlabel('Distance from ' + domain +' (bp)', fontsize=15)
        if name in target_names:
            plt.title('Condensability', fontsize=20)
        else:
            plt.title(name, fontsize=20)
        plt.gca().tick_params(axis='both', which='major', labelsize=15)
        plt.gca().tick_params(axis='both', which='minor', labelsize=15)
        plt.legend(fontsize=14, frameon=False)
        plt.savefig("RPKM_quantile_profile_" + domain + "_" + name + ".svg",
                    format='svg',
                    bbox_inches='tight')

        #plt.show()
        plt.close()

# plot features for each quantiles
if False:
    for i in range(len(q_IDs)):
        fig = plt.figure(figsize=(5,3.5))
        for name, label in zip(names, labels):
            ID_profile = name_ID_profile[name]
            profile = np.nanmean([ID_profile[ID] for ID in q_IDs[i]], axis=0)
            profile = statis.moving_average(profile, moving_average_win)
            profile[:100] = [np.NaN]*100
            profile[len(profile)-100:] = [np.NaN]*100
            X = [ k + offset for k in range(len(profile))]
            #plt.plot(X, profile, label='RPKM quantile '+str(i+1), alpha=0.8)
            plt.plot(X, profile, label=label, alpha=0.8, lw=4)
        if xtick_loc_name:
            xtick_locs, xtick_names = xtick_loc_name
            plt.xticks(xtick_locs, xtick_names, fontsize=10, rotation=45)
        plt.xlabel('Distance from ' + domain +' (bp)', fontsize=15)
        if name in target_names:
            plt.ylabel('Condensability', fontsize=15)
        else:
            plt.ylabel(name, fontsize=15)
        plt.gca().tick_params(axis='both', which='major', labelsize=15)
        plt.gca().tick_params(axis='both', which='minor', labelsize=15)
        plt.legend(fontsize=14, frameon=False, loc='lower right')
        plt.savefig("RPKM_quantile" + str(i) + "_" + domain + "_" + name + ".svg",
                    format='svg',
                    bbox_inches='tight')

        #plt.show()
        plt.close()

# plot profile quantile and heatmap ored by rank score
if True:
    common_IDs = set(ID_rnkvalue1.keys()) & set(ID_rnkvalue2.keys())
    for name in names:
        ID_profile = name_ID_profile[name]
        common_IDs &= set(ID_profile)

    name_ID_dprofile = {}
    for name, label in zip(names, labels):
        if name not in target_names:
            continue
        ID_profile = name_ID_profile[name]
        if label == 'WT':
            control_ID_profile = ID_profile
            continue
        for ID in common_IDs:
            dprofile = ID_profile[ID] - control_ID_profile[ID]
            if name not in name_ID_dprofile:
                name_ID_dprofile[name] = {}
            name_ID_dprofile[name][ID] = dprofile

    # plot quantiles for each features
    q_IDs = statis.quantile({ID:ID_rnkvalue2[ID] for ID in common_IDs}, 5)
    for name, label in zip(names, labels):
        if label == 'WT':
            continue
        if name in target_names:
            ID_dprofile = name_ID_dprofile[name]
        else:
            ID_dprofile = name_ID_profile[name]
        #fig = plt.figure(figsize=(3.5, 2.5))
        #fig = plt.figure(figsize=(4, 3))
        #fig = plt.figure(figsize=(5,3.5))
        fig = plt.figure(figsize=(3,2))
        for i in range(len(q_IDs)):
            dprofile = np.nanmean([ID_dprofile[ID] for ID in q_IDs[i]], axis=0)
            dprofile = statis.moving_average(dprofile, moving_average_win)
            dprofile[:pad_len] = [np.NaN]*pad_len
            dprofile[len(dprofile)-pad_len:] = [np.NaN]*pad_len
            X = [ k + offset for k in range(len(dprofile))]
            plt.plot(X, dprofile, label='Q '+str(i+1), alpha=0.8, lw=3.5)
        if xtick_loc_name:
            xtick_locs, xtick_names = xtick_loc_name
            plt.xticks(xtick_locs, xtick_names, fontsize=10, rotation=45)
            plt.xlim([xtick_locs[0]-pad_len, xtick_locs[-1]+pad_len])
        plt.xlabel('Distance from ' + domain +' (bp)', fontsize=10)
        if name in target_names:
            plt.title('Condensability', fontsize=10)
        else:
            plt.title(name, fontsize=10)
        plt.gca().tick_params(axis='both', which='major', labelsize=10)
        plt.gca().tick_params(axis='both', which='minor', labelsize=10)
        plt.legend(fontsize=11, frameon=False)
        plt.savefig("dzscore_quantile_profile_" + domain + "_" + label + ".svg",
                    format='svg',
                    bbox_inches='tight')

        #plt.show()
        plt.close()

    # heatmap ordred by rnk score
    rnkvalue_ID1 = sorted([(ID_rnkvalue1[ID], ID) for ID in common_IDs])
    rnkvalue_ID2 = sorted([(ID_rnkvalue2[ID], ID) for ID in common_IDs])
    ordered_IDs1 = [ID for rnkvalue, ID in rnkvalue_ID1]
    ordered_IDs2 = [ID for rnkvalue, ID in rnkvalue_ID2]

    label_cmap = {'+inht':'binary_r', 'ODC KO':'binary_r', 'H3K4me3':'Oranges', 'H3K27me3':'Blues'}
    #label_vlims = {'+inht':(-0.5, 0.0), 'ODC KO':(-0.5, 1.0), 'H3K27me3':(0, 100), 'H3K27me3':(0, 200)}
    name_img = {}
    for name, label in zip(names, labels):
        if label == 'WT':
            continue
        if name in target_names:
            ID_dprofile = name_ID_dprofile[name]
        else:
            ID_dprofile = name_ID_profile[name]
        img = []
        for ID in ordered_IDs2:
            dprofile = ID_dprofile[ID]
            dprofile = statis.moving_average(dprofile, moving_average_win) # smoothing
            #dprofile[:pad_len] = [np.NaN]*pad_len
            #dprofile[len(dprofile)-pad_len:] = [np.NaN]*pad_len
            dprofile[:pad_len] = [dprofile[pad_len]]*pad_len
            dprofile[len(dprofile)-pad_len:] = [dprofile[len(dprofile)-pad_len-1]]*pad_len
            img.append(dprofile)
        img = np.asarray(img)
        # rescaling
        if name not in target_names:
            mini = np.min(img)
            img = np.log2(img - mini + 1)
        name_img[name] = img

    for name, label in zip(names[1:], labels[1:]):
        img = name_img[name]
        try:
            cmap = label_cmap[label]
        except:
            cmap = 'viridis'
        try:
            vmin, vmax = label_vlims[label]
        except:
            median = np.nanmedian(img)
            std = np.nanstd(img)
            if name in target_names:
                vmin, vmax = median-0.5*std, median+0.5*std
            else:
                vmin, vmax = None, None
        fig = plt.figure(figsize=(3,5))
        #fig = plt.figure(figsize=(3.5,10))
        #fig = plt.figure(figsize=(4,10))
        #fig = plt.figure(figsize=(5,8))
        #fig = plt.figure(figsize=(7,10))
        #fig = plt.figure(figsize=(3,8))
        plt.imshow(img, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
        #plt.imshow(img, aspect='auto', cmap='jet')
        #plt.title(name)
        #plt.xticks(range(0, len(profile), 1000), xtick_labels, fontsize=15)
        if xtick_loc_name:
            xtick_locs, xtick_names = xtick_loc_name
            new_xtick_locs = [loc-offset for loc in xtick_locs]
            plt.xticks(new_xtick_locs, xtick_names, fontsize=10, rotation=45)
            #plt.xlim([0, int(domain_frac*profile_len)+right_len])
        plt.gca().tick_params(axis='both', which='major', labelsize=10)
        plt.gca().tick_params(axis='both', which='minor', labelsize=10)
        plt.yticks([], [])
        plt.savefig("heatmap_dzscore_profile_" + domain + "_" + label + ".svg",
                    format='svg',
                    bbox_inches='tight')
        #plt.show()
        plt.close()

    


# plot mean profile according to promoter classification
if False:
    common_IDs = set([])
    for name in names:
        if len(common_IDs) <= 0:
            common_IDs = set(name_ID_profile[name].keys())
        common_IDs &= set(name_ID_profile[name].keys())

    pclass_name_profiles = {}
    for ID in common_IDs:
        gname = ID_field_values[ID]['geneName']
        try:
            #pclass = gname_pclass[gname]
            pclass = (gname_pclass[gname][1])
        except:
            continue
        for name in names:
            profile = name_ID_profile[name][ID]
            if pclass not in pclass_name_profiles:
                pclass_name_profiles[pclass] = {}
            if name not in pclass_name_profiles[pclass]:
                pclass_name_profiles[pclass][name] = []
            pclass_name_profiles[pclass][name].append(profile)

    for pclass in pclass_name_profiles:
        fig = plt.figure(figsize=(5,3.5))
        for name, label in zip(names, labels):
            profile = np.nanmean(pclass_name_profiles[pclass][name], axis=0)
            profile = statis.moving_average(profile, moving_average_win)
            profile[:100] = [np.NaN]*100
            profile[len(profile)-100:] = [np.NaN]*100
            X = [ k + offset for k in range(len(profile))]
            #plt.plot(X, profile, label='RPKM quantile '+str(i+1), alpha=0.8)
            plt.plot(X, profile, label=label, alpha=0.8, lw=4)
        if xtick_loc_name:
            xtick_locs, xtick_names = xtick_loc_name
            plt.xticks(xtick_locs, xtick_names, fontsize=10, rotation=45)
        plt.xlabel('Distance from ' + domain +' (bp)', fontsize=15)
        if name in target_names:
            plt.ylabel('Condensability', fontsize=15)
        else:
            plt.ylabel(name, fontsize=15)
        plt.gca().tick_params(axis='both', which='major', labelsize=15)
        plt.gca().tick_params(axis='both', which='minor', labelsize=15)
        #plt.title("%s %s" % pclass)
        plt.title("%s" % pclass)
        #plt.legend(fontsize=14, frameon=False, loc='lower right')
        #plt.savefig("%s-%s_promoter_class.svg" % pclass,
        #            format='svg',
        #            bbox_inches='tight')
        plt.savefig("%s_promoter_class.svg" % pclass,
                    format='svg',
                    bbox_inches='tight')
        #plt.show()
        plt.close()

    name_pclass_dprofile = {}
    for pclass in pclass_name_profiles:
        for name, label in zip(names, labels):
            profile = np.nanmean(pclass_name_profiles[pclass][name], axis=0)
            if label == 'WT':
                control = statis.moving_average(profile, moving_average_win)
                continue
            profile = statis.moving_average(profile, moving_average_win)
            dprofile = profile - control
            if name not in name_pclass_dprofile:
                name_pclass_dprofile[name] = {}
            name_pclass_dprofile[name][pclass] = dprofile

    pclass_list = ['None', 'K4', 'K27', 'K4+K27']
    color_list = ['tab:gray', 'orangered', 'tab:cyan', 'purple']
    for name, label in zip(names[1:], labels[1:]):
        fig = plt.figure(figsize=(5,3.5))
        for pclass, color in zip(pclass_list, color_list):
            dprofile = name_pclass_dprofile[name][pclass]
            dprofile[:100] = [np.NaN]*100
            dprofile[len(dprofile)-100:] = [np.NaN]*100
            X = [ k + offset for k in range(len(dprofile))]
            #plt.plot(X, profile, label='RPKM quantile '+str(i+1), alpha=0.8)
            plt.plot(X, dprofile, color=color, label=pclass, alpha=0.8, lw=4)
        if xtick_loc_name:
            xtick_locs, xtick_names = xtick_loc_name
            plt.xticks(xtick_locs, xtick_names, fontsize=10, rotation=45)
        plt.xlabel('Distance from ' + domain +' (bp)', fontsize=15)
        if name in target_names:
            plt.ylabel('$\Delta$zscore', fontsize=15)
        else:
            plt.ylabel(name, fontsize=15)
        plt.gca().tick_params(axis='both', which='major', labelsize=15)
        plt.gca().tick_params(axis='both', which='minor', labelsize=15)
        #plt.title("%s %s" % pclass)
        plt.title("%s" % (label))
        plt.legend(fontsize=14, frameon=False, loc='lower right')
        plt.savefig("%s_dprofile.svg" % (label),
                    format='svg',
                    bbox_inches='tight')
        #plt.show()
        plt.close()

            

# clustering based on profile shape
if False:
    print
    
        
"""
# plot experiment by experiment
if False:
    # plot mean nucleosome occupancy
    names, occs = [], []
    for name, occ in sorted(name_mocc.items()):
        names.append(name.split('/')[-1])
        occs.append(occ)
    single_plot(occs, names=names, title="Nucleosome occupancy near TSS")


    # mean profile plot (score VS occupancy)
    for score_name in score_names:
        pair_plot(name_mprofile[score_name],
                  name_mocc['-'.join([cell, sample, agent, 0])],
                  offset=offset,
                  xlabel='Distance from ' + domain +' (bp)',
                  ylabel1='Condensability (A.U.)',
                  ylabel2="Occupancy",
                  xtick_loc_name = xtick_loc_name,
                  note="Occ" + "VS" +score_name.split('/')[-1] + domain)


    # mean profile plot (score VS others)
    for score_name in score_names:
        for feature_name in feature_names:
            pair_plot(name_mprofile[score_name],
                      name_mprofile[feature_name],
                      offset=offset,
                      xlabel='Distance from ' + domain +' (bp)',
                      ylabel1='Condensability (A.U.)',
                      ylabel2=feature_name.split('/')[-1],
                      xtick_loc_name = xtick_loc_name,
                      note=score_name.split('/')[-1] + 'VS' + feature_name.split('/')[-1] + domain)


    # paritions according to gene expression
    ID_RPKM = load_file.read_RPKM_new (rna_fname, gtfname, chr_list=[chr_choice]) 
    name_ID_profile['Occ'] = name_ID_occ["mCD8T-KO-NCP-sp-0"]

    for name in name_ID_profile:
        ID_profile = name_ID_profile[name]
        IDs = list(set(ID_profile) & set(ID_RPKM))

        groups = statis.quantile(ID_RPKM, 5)

        fig = plt.figure(figsize=(5,3.5))
        for i in range(len(groups)):
            profile = np.nanmean([ID_profile[ID] for ID in groups[i]], axis=0)
            profile = statis.moving_average(profile, moving_average_win)
            profile[:100] = [np.NaN]*100
            profile[len(profile)-100:] = [np.NaN]*100
            X = [ k + offset for k in range(len(profile))]
            #plt.plot(X, profile, label='RPKM quantile '+str(i+1), alpha=0.8)
            plt.plot(X, profile, label='FPKM quantile '+str(i+1), alpha=0.8, lw=4)
        if xtick_loc_name:
            xtick_locs, xtick_names = xtick_loc_name
            plt.xticks(xtick_locs, xtick_names, fontsize=15, rotation=45)
        plt.xlabel('Distance from ' + feature +' (bp)', fontsize=15)
        if name.split('/')[-1].startswith('-'.join([cell, sample, agent])):
            plt.title('Condensability', fontsize=20)
        else:
            plt.title(name.split('/')[-1], fontsize=20)
        plt.gca().tick_params(axis='both', which='major', labelsize=15)
        plt.gca().tick_params(axis='both', which='minor', labelsize=15)
        plt.legend(fontsize=14, frameon=False)
        plt.savefig("RPKM_quantile_profile_" + domain + "_" + name.split('/')[-1] + ".svg",
                    format='svg',
                    bbox_inches='tight')

        #plt.show()
        plt.close()

    # heatmap plot by gene expression level
    common_IDs = set(ID_RPKM)
    for name in name_list:
        ID_profile = name_ID_profile[name]
        common_IDs &= set(ID_profile)

    RPKM_ID = sorted([(ID_RPKM[ID], ID) for ID in common_IDs])

    name_img = {}
    for name in name_list:
        ID_profile = name_ID_profile[name]
        img = []
        # smoothing
        for RPKM, ID in RPKM_ID:
            profile = ID_profile[ID]
            profile = statis.moving_average(profile, moving_average_win)
            profile[:100] = [profile[100]]*100
            profile[len(profile)-100:] = [profile[len(profile)-101]]*100
            img.append(profile)
        img = np.asarray(img)
        # rescaling
        if name.split('/')[-1].startswith('-'.join([cell, sample, agent])):
            mini = np.min(img)
            img = np.log2(img - mini + 1)
        else:
            img = 2*img

        name_img[name] = img
    
    for i in range(len(name_list)):
        name = name_list[i]
        img = name_img[name]
        #fig = plt.figure(figsize=(5,10))
        #fig = plt.figure(figsize=(6.5,10))
        fig = plt.figure(figsize=(5,8))
        plt.imshow(img, aspect='auto', cmap=cmap_list[i], vmin=vlim_list[i][0], vmax=vlim_list[i][1])
        #plt.imshow(img, aspect='auto', cmap='jet')
        #plt.title(name)
        if xtick_loc_name:
            xtick_locs, xtick_names = xtick_loc_name
            plt.xticks(np.asarray(xtick_locs) - offset, xtick_names, fontsize=15, rotation=45)
        else:
            xtick_labels = [ str(k + offset) for k in range(0, len(profile), 1000)]
            plt.xticks(range(0, len(profile), 1000), xtick_labels, fontsize=15)
        plt.gca().tick_params(axis='both', which='major', labelsize=15)
        plt.gca().tick_params(axis='both', which='minor', labelsize=15)
        plt.savefig("heatmap_profile_byRPKM_" + domain + "_" + name.split('/')[-1] + ".svg",
                    format='svg',
                    bbox_inches='tight')
        #plt.show()
        plt.close()
        
"""
