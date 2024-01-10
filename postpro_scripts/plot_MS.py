import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

def tuple_cmp (a, b):
    if a[0] < b[0]:
        return -1
    elif a[0] > b[0]:
        return 1
    else:
        if a[1] < b[1]:
            return -1
        elif a[1] > b[1]:
            return 1
        else:
            return 0


def sort_dict (dic, reverse=False):
    value_key = sorted([(value, key) for key, value in dic.items()], cmp=tuple_cmp, reverse=reverse)
    sorted_key = [key for value, key in value_key]
    return sorted_key

def read_MS (fname):
    sample_peptide_ptm_ratio = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        if First:
            fields = line.split('\t')
            First = False
            continue
        if line.startswith('Peptide'):
            col_choices = []
            cols = line.split('\t')
            assert len(cols) == len(fields) + 1
            for i in range(len(cols)):
                if cols[i].strip() == 'Ratio':
                    col_choices.append(i)
            continue
        cols = line.split('\t')
        if not cols[0].startswith('H'):
            peptide = cols[0]
            continue
        else:
            try:
                float(cols[1])
            except:
                peptide = cols[0]
                continue

        ptm = cols[0]
        for idx in col_choices:
            ratio = cols[idx]
            #print peptide, ptm, ratio
            ratio = float(ratio)
            sample = fields[idx-1]
            if sample not in sample_peptide_ptm_ratio:
                sample_peptide_ptm_ratio[sample] = {}
            if peptide not in sample_peptide_ptm_ratio[sample]:
                sample_peptide_ptm_ratio[sample][peptide] = {}
            #print ptm
            assert ptm not in sample_peptide_ptm_ratio[sample][peptide]
            sample_peptide_ptm_ratio[sample][peptide][ptm] = ratio

    return sample_peptide_ptm_ratio

def parse_peptide (peptide):
    seq, histone_info = peptide[:-1].split('(')
    histone, st, ed = histone_info.split('_')
    st, ed = int(st), int(ed)
    return histone, st, ed, seq

def parse_ptm (ptm_string):
    pattern = '([A-Z][0-9]+(?:ac|me|ph|ub)[1-3]?)'
    ptm_list = re.findall(pattern, ptm_string)
    return ptm_list

sample_peptide_ptm_ratio = read_MS("Sangwoo_histone_ratios_092122.csv")


# aggregate data into single ptm ratio
histone_singleptm_ratios = {}
histones = ['H3', 'H4']
for peptide in sample_peptide_ptm_ratio['input'].keys():
    histone, st, ed, seq = parse_peptide(peptide)

    if histone not in histones:
        continue

    label_score1, label_score2 = {}, {}
    for ptm in sample_peptide_ptm_ratio['input'][peptide]:
        _, ptm_string = ptm.split(' ')

        singleptm_list = parse_ptm (ptm_string)

        input_ratio = sample_peptide_ptm_ratio['input'][peptide][ptm]
        sup_ratio = sample_peptide_ptm_ratio['supernatant'][peptide][ptm]
        pel_ratio = sample_peptide_ptm_ratio['pellet'][peptide][ptm]

        for singleptm in singleptm_list:

            if histone not in histone_singleptm_ratios:
                histone_singleptm_ratios[histone] = {}
            if singleptm not in histone_singleptm_ratios[histone]:
                histone_singleptm_ratios[histone][singleptm] = [0.0, 0.0, 0.0]

            histone_singleptm_ratios[histone][singleptm][0] += input_ratio
            histone_singleptm_ratios[histone][singleptm][1] += pel_ratio
            histone_singleptm_ratios[histone][singleptm][2] += sup_ratio
            

ptmname_folds = {}
ptmname_state = {}
for histone in histone_singleptm_ratios:
    for singleptm in histone_singleptm_ratios[histone]:
        ptmname = histone + singleptm
        ratios = histone_singleptm_ratios[histone][singleptm]

        if sum(ratios) <= 0:
            continue

        mean = np.mean(ratios)
        std = np.std(ratios)
        state = [float(ratio-mean)/std for ratio in ratios]
        assert ptmname not in ptmname_state
        ptmname_state[ptmname] = state

        if ratios[0] <= 0:
            continue
        if ratios[1] <= 0:
            continue
        if ratios[2] <= 0:
            continue

        pel_fold = np.log2(float(ratios[1])/ratios[0])
        sup_fold = np.log2(float(ratios[2])/ratios[0])
        assert ptmname not in ptmname_folds
        ptmname_folds[ptmname] = [pel_fold, sup_fold]


# plot fold change heatamp

fold_ptmname = sorted([(folds[1], ptmname) for ptmname, folds in ptmname_folds.items()], reverse=True)
ptmname_list = [ptmname for fold, ptmname in fold_ptmname]
data = []
for ptmname in ptmname_list:
    folds = ptmname_folds[ptmname]
    data.append(folds)

fig = plt.figure()
plt.imshow(data, cmap='coolwarm', vmin=-1.5, vmax=1.5)
plt.yticks(range(len(ptmname_list)), ptmname_list)
plt.colorbar()
#plt.show()
plt.close()


        
#hierarchical clustering heatmap

ptmname_list = ptmname_state.keys()
data = []
for ptmname in ptmname_list:
    state = ptmname_state[ptmname]
    data.append(state)

pastel_jet = LinearSegmentedColormap.from_list('white_viridis',
                                             [(0, 'tab:blue'),
                                              (0.2, 'tab:cyan'),
                                              (0.5, 'white'),
                                              (0.8, 'tomato'),
                                              (1, 'tab:red'),
                                             ], N=256)


hmap = sns.clustermap(data,
                      method='average',
                      metric='euclidean',
                      figsize=(3,10),
                      cbar_kws=None,
                      row_cluster=True,
                      col_cluster=False,
                      dendrogram_ratio=0.2,
                      colors_ratio=0.03,
                      tree_kws=None,
                      cmap='Spectral_r',
                      center=0,
                      xticklabels=['Input', 'Pellet', 'Supnt'],
                      annot=True,
                      annot_kws={"size": 8})
                      #cbar_pos=None)

new_labels = []
for tick_label in hmap.ax_heatmap.axes.get_yticklabels():
    idx = int(tick_label.get_text())
    ptmname = ptmname_list[idx]
    tick_label.set_text(ptmname)
    new_labels.append(tick_label)
    #tick_label.set_color(lut[species_name])

hmap.ax_heatmap.axes.set_yticklabels(new_labels)

plt.gca().xaxis.tick_top()
plt.xticks(rotation=70)
plt.savefig("hmap.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()


sys.exit(1)

sns.clustermap(data,
               method='average',
               metric='euclidean',
               figsize=(10, 10),
               cbar_kws=None,
               row_cluster=True,
               col_cluster=False,
               row_linkage=None,
               col_linkage=False,
               row_colors=None,
               col_colors=False,
               dendrogram_ratio=0.2,
               colors_ratio=0.03,
               cbar_pos=(0.02, 0.8, 0.05, 0.18),
               tree_kws=None)
            



histones = ['H2A', 'H2B', 'H1']

small = 10**(-10)
for peptide in sample_peptide_ptm_ratio['input'].keys():
    histone, st, ed, seq = parse_peptide(peptide)

    if histone not in histones:
        continue

    label_score1, label_score2 = {}, {}
    for ptm in sample_peptide_ptm_ratio['input'][peptide]:
        _, only_ptm = ptm.split(' ')

        input_ratio = sample_peptide_ptm_ratio['input'][peptide][ptm]
        sup_ratio = sample_peptide_ptm_ratio['supernatant'][peptide][ptm]
        pel_ratio = sample_peptide_ptm_ratio['pellet'][peptide][ptm]

        if input_ratio + sup_ratio + pel_ratio <= 0:
            continue

        input_ratio +=small
        sup_ratio +=small
        pel_ratio +=small

        sup_score = np.log2(sup_ratio/input_ratio)
        pel_score = np.log2(pel_ratio/input_ratio)

        label = '%s %s' % (histone, only_ptm)

        label_score1[label] = sup_score
        label_score2[label] = pel_score

    if len(label_score1) <= 0:
        continue

    sorted_labels = sort_dict(label_score1, reverse=True)

    img = []
    for label in sorted_labels:
        img.append([label_score1[label], label_score2[label]])


    row_count = len(img)
    col_count = len(img[0])

    cell_width, cell_height = 0.3, 0.2 #inch
    marg_top = 0.5
    marg_bott = 0.5
    marg_left = 2
    
    marg_right = 0.2

    fig_width = cell_width*col_count + marg_left + marg_right
    fig_height = cell_height*row_count + marg_top + marg_bott

    fig = plt.figure(figsize=(fig_width, fig_height))
    
    # adjust margins (relative numbers) according to absolute values
    fig.subplots_adjust(bottom =marg_bott/fig_height ,top=1.-marg_top/fig_height,
                        left=marg_left/fig_width, right=1.-marg_right/fig_width)

    im = plt.imshow(img, cmap='seismic', aspect='auto', vmin=-4, vmax=4)
    plt.yticks(range(len(sorted_labels)), sorted_labels)
    plt.xticks([0,1], ['sup', 'pel'])
    title = '%s %d-%d' % (histone, st, ed)
    plt.title(title)
    #plt.colorbar()
    #plt.tight_layout()
    plt.savefig(title + '.svg', format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

# save colorbar
fig = plt.figure(figsize=(1.1,1.2))
plt.subplot(1,2,1)
cbar = plt.colorbar(im, cax=plt.gca(), ticks=[-4, 4])
cbar.ax.set_yticklabels(['-4', '4'], fontsize=8)
cbar.ax.set_ylabel('log2 Fold change', rotation=-90, va="bottom", fontsize=8)
plt.tight_layout()
plt.savefig('MS_cbar.svg', format='svg', bbox_inches='tight')
plt.close()



"""
for peptide in sample_peptide_ptm_ratio['input'].keys():
    histone, st, ed, seq = parse_peptide(peptide)

    fig = plt.figure()

    input_bott, sup_bott, pel_bott = 0.0, 0.0, 0.0    
    for ptm in sample_peptide_ptm_ratio['input'][peptide]:
        _, only_ptm = ptm.split(' ')

        input_ratio = sample_peptide_ptm_ratio['input'][peptide][ptm]        
        sup_ratio = sample_peptide_ptm_ratio['supernatant'][peptide][ptm]
        pel_ratio = sample_peptide_ptm_ratio['pellet'][peptide][ptm]

        plt.bar(['Input', 'Supernatant', 'Pellet'],
                [input_ratio, sup_ratio, pel_ratio],
                bottom=[input_bott, sup_bott, pel_bott],
                label='%s %s' % (histone, only_ptm), width=0.4)

        input_bott +=input_ratio
        sup_bott +=sup_ratio
        pel_bott +=pel_ratio

    plt.title('%s %d-%d' % (histone, st, ed))
    plt.ylim([0.0, 1.1])
    plt.show()
    plt.close()

histone_ptm_sample_score = {}
histone_frag_ptm = {}

histones = ['H3', 'H4']
addon = 10**(-10)

for peptide in sample_peptide_ptm_ratio['input'].keys():
    histone, st, ed, seq = parse_peptide(peptide)
    if histone not in histones:
        continue
    frag = (st, ed)
    for ptm in sample_peptide_ptm_ratio['input'][peptide]:
        _, only_ptm = ptm.split(' ')

        if only_ptm == 'unmod':
            continue

        input_ratio = sample_peptide_ptm_ratio['input'][peptide][ptm]

        if input_ratio <=0:
            continue
        
        sup_ratio = sample_peptide_ptm_ratio['supernatant'][peptide][ptm]
        pel_ratio = sample_peptide_ptm_ratio['pellet'][peptide][ptm]

        sup_score = np.log2((addon+sup_ratio) / (addon+input_ratio))
        pel_score = np.log2((addon+pel_ratio) / (addon+input_ratio))

        #sup_score = np.log2(1.0 + float(sup_ratio)/input_ratio)
        #pel_score = np.log2(1.0 + float(pel_ratio)/input_ratio)

        #sup_score = float(sup_ratio) / input_ratio
        #pel_score = float(pel_ratio) / input_ratio

        if histone not in histone_ptm_sample_score:
            histone_ptm_sample_score[histone] = {}

        assert only_ptm not in histone_ptm_sample_score[histone]
        histone_ptm_sample_score[histone][only_ptm] = {}
        histone_ptm_sample_score[histone][only_ptm]['sup'] = sup_score
        histone_ptm_sample_score[histone][only_ptm]['pel'] = pel_score

        if histone not in histone_frag_ptm:
            histone_frag_ptm[histone] = {}
        if frag not in histone_frag_ptm[histone]:
            histone_frag_ptm[histone][frag] = []
        histone_frag_ptm[histone][frag].append(only_ptm)


labels = []
img = []
for histone in histones:
    for frag in sorted(histone_frag_ptm[histone]):
        for ptm in histone_frag_ptm[histone][frag]:
           labels.append(' '.join([histone, ptm]))
           sup_score = histone_ptm_sample_score[histone][ptm]['sup']
           pel_score = histone_ptm_sample_score[histone][ptm]['pel']
           img.append([sup_score, pel_score])

fig = plt.figure()
plt.imshow(img, cmap='seismic', aspect='auto', vmin=-5, vmax=5)
plt.yticks(range(len(labels)), labels)
plt.colorbar()
plt.tight_layout()
plt.show()



#GM_peptide_ptm_ratio = read_MS("Sangwoo_human_062222.csv")
#HeLa_peptide_ptm_ratio = read_MS("HeLa_histone_ratios_Jan2021.csv")

for HeLa in HeLa_peptide_ptm_ratio:
    for peptide in sorted(HeLa_peptide_ptm_ratio[HeLa]):
        fig = plt.figure()
        ptms = sorted(HeLa_peptide_ptm_ratio[HeLa][peptide].keys())
        HeLa_bottom, GM_AE_bottom, GM_Heat_bottom = 0, 0, 0
        for ptm in ptms:
            HeLa_ratio = HeLa_peptide_ptm_ratio[HeLa][peptide][ptm]
            GM_AE_ratio = GM_peptide_ptm_ratio['1,20220622_Bhanu_HaSangwoo_GM_AE'][peptide][ptm]
            GM_Heat_ratio = GM_peptide_ptm_ratio['2,20220622_Bhanu_HaSangwoo_GM_Heat'][peptide][ptm]
            plt.bar(range(3),
                    [HeLa_ratio, GM_AE_ratio, GM_Heat_ratio],
                    bottom=[HeLa_bottom, GM_AE_bottom, GM_Heat_bottom],
                    label=ptm)
            #print HeLa_ratio
            HeLa_bottom += HeLa_ratio
            GM_AE_bottom += GM_AE_ratio
            GM_Heat_bottom += GM_Heat_ratio
        plt.xticks(range(3), ['HeLa', 'GM AE', 'GM Heat'], rotation=45)
        #plt.xlabel("Sample")
        plt.ylabel("Ratios")
        plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
        #plt.legend(loc='upper right')
        plt.title(peptide)
        plt.savefig("%s.png" % (peptide), bbox_inches='tight')
        #plt.show()
        plt.close()



# profile plot
for HeLa in HeLa_peptide_ptm_ratio:
    for peptide in sorted(HeLa_peptide_ptm_ratio[HeLa]):
        fig = plt.figure()
        ptms = sorted(HeLa_peptide_ptm_ratio[HeLa][peptide].keys())
        HeLa_profile, GM_AE_profile, GM_Heat_profile = [], [], []
        for ptm in ptms:
            HeLa_ratio = HeLa_peptide_ptm_ratio[HeLa][peptide][ptm]
            GM_AE_ratio = GM_peptide_ptm_ratio['1,20220622_Bhanu_HaSangwoo_GM_AE'][peptide][ptm]
            GM_Heat_ratio = GM_peptide_ptm_ratio['2,20220622_Bhanu_HaSangwoo_GM_Heat'][peptide][ptm]
            HeLa_profile.append(HeLa_ratio)
            GM_AE_profile.append(GM_AE_ratio)
            GM_Heat_profile.append(GM_Heat_ratio)
        plt.plot(range(len(ptms)), HeLa_profile, 'ko--', alpha=0.25, label='HeLa')
        plt.plot(range(len(ptms)), GM_AE_profile, 'bo--', alpha=0.5, label='GM AE')
        plt.plot(range(len(ptms)), GM_Heat_profile, 'ro--', alpha=0.5, label='GM Heat')
        plt.xticks(range(len(ptms)), ptms, rotation=45)
        plt.xlabel("PTMs")
        plt.ylabel("Ratios")
        plt.legend()
        plt.title(peptide)
        plt.show()
        plt.close()

"""
