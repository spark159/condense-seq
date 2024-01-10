import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import copy
import scipy

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

# reorganize the data by replicates and get average value
def reorganize_data (sample_peptide_ptm_ratio):
    # organize the data by replicate
    sample_peptide_ptm_rep_ratio = {}
    for sample in sample_peptide_ptm_ratio:
        peptide_ptm_ratio = sample_peptide_ptm_ratio[sample]
        cols = sample.strip().split('_')
        new_sample, rep = cols[-2], cols[-1]
        rep = int(rep)
        for peptide in peptide_ptm_ratio:
            for ptm  in peptide_ptm_ratio[peptide]:
                ratio = peptide_ptm_ratio[peptide][ptm]
                if new_sample not in sample_peptide_ptm_rep_ratio:
                    sample_peptide_ptm_rep_ratio[new_sample] = {}
                if peptide not in sample_peptide_ptm_rep_ratio[new_sample]:
                    sample_peptide_ptm_rep_ratio[new_sample][peptide] = {}
                if ptm not in sample_peptide_ptm_rep_ratio[new_sample][peptide]:
                    sample_peptide_ptm_rep_ratio[new_sample][peptide][ptm] = {}
                assert rep not in sample_peptide_ptm_rep_ratio[new_sample][peptide][ptm]
                sample_peptide_ptm_rep_ratio[new_sample][peptide][ptm][rep] = ratio

    # average the replicates
    sample_peptide_ptm_mratio = {}
    for sample in sample_peptide_ptm_rep_ratio:
        for peptide in sample_peptide_ptm_rep_ratio[sample]:
            for ptm in sample_peptide_ptm_rep_ratio[sample][peptide]:
                ratios = sample_peptide_ptm_rep_ratio[sample][peptide][ptm].values()
                if sample not in sample_peptide_ptm_mratio:
                    sample_peptide_ptm_mratio[sample] = {}
                if peptide not in sample_peptide_ptm_mratio[sample]:
                    sample_peptide_ptm_mratio[sample][peptide] = {}
                sample_peptide_ptm_mratio[sample][peptide][ptm] = np.mean(ratios)

    return sample_peptide_ptm_rep_ratio, sample_peptide_ptm_mratio


def get_score (sample_peptide_ptm_ratio):
    #histones = ['H2A', 'H2B', 'H1']
    histones = ['H3', 'H4']
    small = 10**(-10)
    label_score1, label_score2 = {}, {}
    for peptide in sample_peptide_ptm_ratio['input'].keys():
        histone, st, ed, seq = parse_peptide(peptide)

        if histone not in histones:
            continue

        for ptm in sample_peptide_ptm_ratio['input'][peptide]:
            _, only_ptm = ptm.split(' ')

            input_ratio = sample_peptide_ptm_ratio['input'][peptide][ptm]
            #sup_ratio = sample_peptide_ptm_ratio['supernatant'][peptide][ptm]
            sup_ratio = sample_peptide_ptm_ratio['sup'][peptide][ptm]
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

    return label_score1, label_score2


def draw_data (sample_peptide_ptm_ratio, histones=['H3', 'H4']):
    small = 10**(-10)
    for peptide in sample_peptide_ptm_ratio['input'].keys():
        histone, st, ed, seq = parse_peptide(peptide)

        if histone not in histones:
            continue

        label_score1, label_score2 = {}, {}
        for ptm in sample_peptide_ptm_ratio['input'][peptide]:
            _, only_ptm = ptm.split(' ')

            input_ratio = sample_peptide_ptm_ratio['input'][peptide][ptm]
            #sup_ratio = sample_peptide_ptm_ratio['supernatant'][peptide][ptm]
            sup_ratio = sample_peptide_ptm_ratio['sup'][peptide][ptm]
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


# read data
sample_peptide_ptm_ratio0 = read_MS("Sangwoo_histone_ratios_092122.csv")
sample_peptide_ptm_ratio1 = read_MS("Sangwoo_histone_ratios_old_full.csv")
sample_peptide_ptm_ratio2 = read_MS("Sangwoo_histone_ratios_new_full.csv")

# reorganizing data
data_list = [sample_peptide_ptm_ratio0,
             sample_peptide_ptm_ratio1,
             sample_peptide_ptm_ratio2]

sample_ptm_value = {}
for i in range(len(data_list)):
    sample_peptide_ptm_ratio = data_list[i]
    for sample in sample_peptide_ptm_ratio:
        cols = sample.strip().split('_')
        if len(cols) <= 1:
            condition  = cols[0]
            rep = 0
        else:
            condition, rep = cols[-2], cols[-1]
            rep = int(rep)

        if i <= 1:
            date = 'Old'
        else:
            date = 'New'

        newsample = '%s-%s-%d' % (date, condition, rep)
        peptide_ptm_ratio = sample_peptide_ptm_ratio[sample]
        for peptide in peptide_ptm_ratio:
            for ptm in peptide_ptm_ratio[peptide]:
                if newsample not in sample_ptm_value:
                    sample_ptm_value[newsample] = {}
                assert ptm not in sample_ptm_value[newsample]
                ratio = peptide_ptm_ratio[peptide][ptm]
                sample_ptm_value[newsample][ptm] = ratio

# get fold change
def get_fold_change(sample_ptm_value, small=10**(-10)):
    sample_ptm_fold = {}
    for date in ['Old', 'New']:
        if date == 'Old':
            rep_list = range(4)
        else:
            rep_list = range(1,4)
        for rep in rep_list:
            control_sample = '%s-%s-%d' % (date, 'input', rep)
            pellet_sample = '%s-%s-%d' % (date, 'pellet', rep)
            sup_sample = '%s-%s-%d' % (date, 'sup', rep)
            for ptm in sample_ptm_value[control_sample]:
                control_value = sample_ptm_value[control_sample][ptm]
                pellet_value = sample_ptm_value[pellet_sample][ptm]
                sup_value = sample_ptm_value[sup_sample][ptm]

                if control_value*pellet_value*sup_value <=0:
                    continue
                
                #if control_value + pellet_value + sup_value <= 0:
                #    continue
                #control_value += small
                #pellet_value += small
                #sup_value += small

                pellet_fold = np.log2(float(pellet_value)/control_value)
                sup_fold = np.log2(float(sup_value)/control_value)

                if pellet_sample not in sample_ptm_fold:
                    sample_ptm_fold[pellet_sample] = {}
                sample_ptm_fold[pellet_sample][ptm] = pellet_fold

                if sup_sample not in sample_ptm_fold:
                    sample_ptm_fold[sup_sample] = {}
                sample_ptm_fold[sup_sample][ptm] = sup_fold

    return sample_ptm_fold
sample_ptm_fold = get_fold_change(sample_ptm_value)

#sys.exit(1)                
            
# QC check correlation between replicates
def plot_corr (sample_ptm_value, sample_list=None, note=''):
    if sample_list == None:
        sample_list = sorted(sample_ptm_value.keys())
    
    pair_corr = {}
    fig, axes = plt.subplots(figsize=(7,6), nrows=len(sample_list), ncols=len(sample_list))
    for i in range(len(sample_list)):
        for j in range(len(sample_list)):
            sample1, sample2 = sample_list[i], sample_list[j]
            ptm_value1, ptm_value2 = sample_ptm_value[sample1], sample_ptm_value[sample2]
            ptm_list = list(set(ptm_value1.keys()) & set(ptm_value2.keys()))
            X = [ptm_value1[ptm] for ptm in ptm_list]
            Y = [ptm_value2[ptm] for ptm in ptm_list]
            if i > j:
                axes[i,j].plot(X, Y, 'k.', markersize=0.5)
                wspace = 0.1*(max(X) - min(X))
                hspace = 0.1*(max(Y) - min(Y))
                axes[i,j].set_xticks([min(X), max(X)])
                axes[i,j].set_xticklabels([str(round(min(X),1)), str(round(max(X),1))], rotation=45)
                axes[i,j].set_yticks([min(Y), max(Y)])
                axes[i,j].set_yticklabels([str(round(min(Y),1)), str(round(max(Y),1))])
                axes[i,j].set_xlim(min(X)-wspace, max(X)+wspace)
                axes[i,j].set_ylim(min(Y)-hspace, max(Y)+hspace)
                axes[i,j].tick_params(axis='both', which='major', labelsize=5)
                axes[i,j].tick_params(axis='both', which='minor', labelsize=5)

                if j > 0 and i < len(sample_list) -1:
                    axes[i,j].tick_params(axis='both', which='both', labelbottom=False, labelleft=False)
                if j == 0 and i < len(sample_list) -1:
                    axes[i,j].tick_params(axis='x', which='both', labelbottom=False)
                if j > 0 and i == len(sample_list) - 1:
                    axes[i,j].tick_params(axis='y', which='both', labelleft=False)

            elif i == j:
                matrix = np.zeros((len(sample_list), len(sample_list)))
                matrix[:] = np.nan
                axes[i,j].imshow(matrix, origin='lower')
                s = sample1
                axes[i,j].text(len(sample_list)/2, len(sample_list)/2, s, ha="center", va="center", fontsize=5, weight='bold')
                axes[i,j].set_xlim([0, len(sample_list)-1])
                axes[i,j].set_ylim([0, len(sample_list)-1])
                axes[i,j].set_axis_off()
            else:
                assert i < j
                corr = scipy.stats.spearmanr(X, Y)[0]
                #corr = scipy.stats.pearsonr(X, Y)[0]
                assert (sample1, sample2) not in pair_corr
                pair_corr[(sample1, sample2)] = corr
                matrix = np.zeros((len(sample_list), len(sample_list)))
                matrix[:] = corr
                img = axes[i,j].imshow(matrix, cmap="jet", vmin=0.0, vmax=1.0, origin='lower')
                if corr < 0.3 or corr > 0.7:
                    color = "white"
                else:
                    color = "black"
                axes[i,j].text(len(sample_list)/2, len(sample_list)/2, str(round(corr,2)), ha="center", va="center", fontsize=6, color=color, weight='bold')
                axes[i,j].set_xlim([0, len(sample_list)-1])
                axes[i,j].set_ylim([0, len(sample_list)-1])
                axes[i,j].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    plt.subplots_adjust(wspace=0.2, hspace=0.2)
    cbar=fig.colorbar(img, ax=axes, location='right', shrink=0.6, aspect=30, ticks=[0, 1.0])
    cbar.ax.set_yticklabels([str(0), str(1)], fontsize=5)
    cbar.ax.set_ylabel('Spearman correlation', rotation=-90, va="bottom", fontsize=5, labelpad=-5)
    #plt.suptitle("Corrrelation betwen condensing samples")
    #plt.show()
    #plt.savefig("MS_corr_btw_sample.svg", format='svg', bbox_inches='tight')
    #plt.savefig("MS_corr_btw_sample.png", dpi=300, bbox_inches='tight')
    plt.savefig("MS_corr_btw_sample_" + note + ".png", dpi=300, bbox_inches='tight')
    plt.close()
    return pair_corr

sample_list = sorted(sample_ptm_value.keys())
old_sample_list = []
new_sample_list = []
for sample in sample_list:
    date, condition, rep = sample.split('-')
    rep = int(rep)
    if condition == 'input':
        continue
    if sample.startswith('New'):
        new_sample_list.append(sample)
    else:
        old_sample_list.append(sample)

plot_corr(sample_ptm_fold, old_sample_list, note='old')
plot_corr(sample_ptm_fold, new_sample_list, note='new')
#sys.exit(1)




# aggregate data into single ptms
def aggregate_data(sample_ptm_value, histone_list = ['H3', 'H4']):
    sample_sptm_value = {}
    for sample in sample_ptm_value:
        for ptm in sample_ptm_value[sample]:
            peptide_string, ptm_string = ptm.split(' ')
            histone = peptide_string.split('_')[0]
            if histone not in histone_list:
                continue
            value = sample_ptm_value[sample][ptm]
            for sptm in parse_ptm (ptm_string):
                label = ' '.join([histone, sptm])
                if sample not in sample_sptm_value:
                    sample_sptm_value[sample] = {}
                if label not in sample_sptm_value[sample]:
                    sample_sptm_value[sample][label] = 0.0
                sample_sptm_value[sample][label] += value
    return sample_sptm_value         
sample_sptm_value = aggregate_data(sample_ptm_value)
sample_sptm_fold = get_fold_change(sample_sptm_value)
plot_corr(sample_ptm_fold, old_sample_list, note='old_single')
plot_corr(sample_ptm_fold, new_sample_list, note='new_single')


#hierarchical clustering heatmap
def plot_hcluster (name_state,
                  figsize = (3,10),
                  xticklabels=None,
                  cmap='Spectral_r',
                  note=""):

    name_list = name_state.keys()
    data = []
    for name in name_list:
        state = name_state[name]
        data.append(state)

    if xticklabels:
        assert len(data[0]) == len(xticklabels)
    else:
        xticklabels = [str(i) for i in range(len(data[0]))]

    hmap = sns.clustermap(data,
                          method='average',
                          metric='euclidean',
                          figsize=figsize,
                          cbar_kws=None,
                          row_cluster=True,
                          col_cluster=False,
                          dendrogram_ratio=0.2,
                          colors_ratio=0.03,
                          tree_kws=None,
                          cmap=cmap,
                          center=0,
                          vmin=-0.5,
                          vmax=0.5,
                          xticklabels=xticklabels,
                          yticklabels=1,
                          annot=True,
                          annot_kws={"size": 8},
                          cbar_pos=None)

    new_labels = []
    for tick_label in hmap.ax_heatmap.axes.get_yticklabels():
        idx = int(tick_label.get_text())
        name = name_list[idx]
        tick_label.set_text(name)
        new_labels.append(tick_label)
        #tick_label.set_color(lut[species_name])

    hmap.ax_heatmap.axes.set_yticklabels(new_labels)

    plt.gca().xaxis.tick_top()
    plt.xticks(rotation=70)
    plt.yticks(rotation=0)
    plt.savefig("hmap_" + note + ".svg", format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

#date = 'Old'
#old_sample_list = []
#for condition in ['input', 'pellet', 'sup']:
#    for rep in range(0,4):
#        sample = '%s-%s-%s' % (date, condition, rep)
#        old_sample_list.append(sample)

    
old_sptm_state = {}
for sptm in sample_sptm_fold['Old-pellet-0']:
    state = []
    for sample in old_sample_list:
        try:
            fold = sample_sptm_fold[sample][sptm]
            state.append(fold)
        except:
            break
    #if sum(state) <= 0:
    #x    continue
    if len(state) == len(old_sample_list):
        #mean = np.mean(state)
        #std = np.std(state)
        #state = [float(value-mean)/std for value in state]
        old_sptm_state[sptm] = state

plot_hcluster(old_sptm_state, xticklabels=old_sample_list, cmap='bwr')
        

# taking weighted average over replicates
def taking_average (sample_ptm_value, date_choice, rep_weight=None):
    condition_ptm_values = {}
    for sample in sample_ptm_value:
        for ptm in sample_ptm_value[sample]:
            date, condition, rep = sample.split('-')
            rep = int(rep)
            if date != date_choice:
                continue
            if rep_weight and rep not in rep_weight:
                continue
            if condition not in condition_ptm_values:
                condition_ptm_values[condition] = {}
            if ptm not in condition_ptm_values:
                condition_ptm_values[condition][ptm] = []
            value = sample_ptm_value[sample][ptm]
            condition_ptm_values[condition][ptm].append(rep_weight[rep]*value)
    condition_ptm_mvalue = {}
    for condition in condition_ptm_values:
        for ptm in condition_ptm_values[condition]:
            if condition not in condition_ptm_mvalue:
                condition_ptm_mvalue[condition] = {}
            condition_ptm_mvalue[condition][ptm] = np.mean(condition_ptm_values[condition][ptm])
    return condition_ptm_mvalue
cdt_ptm_value0 = taking_average (sample_ptm_value, 'Old', {0:1})
cdt_ptm_value1 = taking_average (sample_ptm_value, 'Old', {1:1, 2:1, 3:1})

cdt_sptm_value0 = aggregate_data (cdt_ptm_value0)
cdt_sptm_value1 = aggregate_data (cdt_ptm_value1)

conditions = ['input', 'pellet', 'sup']
sptm_state0 = {}
for sptm in cdt_sptm_value0['input']:
    state = []
    for cdt in conditions:
        state.append(cdt_sptm_value0[cdt][sptm])
    if sum(state) <= 0:
        continue
    mean = np.mean(state)
    std = np.std(state)
    state = [float(value-mean)/std for value in state]
    sptm_state0[sptm] = state

sptm_state1 = {}
for sptm in cdt_sptm_value1['input']:
    state = []
    for cdt in conditions:
        state.append(cdt_sptm_value1[cdt][sptm])
    if sum(state) <= 0:
        continue
    mean = np.mean(state)
    std = np.std(state)
    state = [float(value-mean)/std for value in state]    
    sptm_state1[sptm] = state

plot_hcluster(sptm_state0, xticklabels = conditions, note='first')
plot_hcluster(sptm_state1, xticklabels = conditions, note='old_average')

# get fold change and p-value for 3 replica old data







sys.exit(1)




#plot_hcluster(cdt_sptm_value0)
#plot_hcluster(cdt_sptm_value1)








# reorganizing data by replicates
def reorganize (sample_peptide_ptm_ratio):
    rep_data = {}
    for sample in sample_peptide_ptm_ratio:
        cols = sample.strip().split('_')
        newsample, rep = cols[-2], cols[-1]
        rep = int(rep)
        if rep not in rep_data:
            rep_data[rep] = {}
        assert newsample not in rep_data[rep]
        peptide_ptm_ratio = sample_peptide_ptm_ratio[sample]
        rep_data[rep][newsample] = copy.deepcopy(peptide_ptm_ratio)
    return rep_data
rep_data1 = reorganize(sample_peptide_ptm_ratio1)
rep_data1[0] = sample_peptide_ptm_ratio0
rep_data2 = reorganize(sample_peptide_ptm_ratio2)


# QC:correlation between replicates
def check_replicates (rep_data, title=None):
    rep_list = sorted(rep_data.keys())
    for i in range(len(rep_list)-1):
        for j in range(i+1, len(rep_list)):
            rep1 = rep_list[i]
            rep2 = rep_list[j]
            sample_peptide_ptm_ratio1 = rep_data[rep1]
            sample_peptide_ptm_ratio2 = rep_data[rep2]

            X, Y = [], []
            for sample in sample_peptide_ptm_ratio1:
                for peptide in sample_peptide_ptm_ratio1[sample]:
                    for ptm in sample_peptide_ptm_ratio1[sample][peptide]:
                        ratio1 = sample_peptide_ptm_ratio1[sample][peptide][ptm]
                        ratio2 = sample_peptide_ptm_ratio2[sample][peptide][ptm]
                        X.append(ratio1)
                        Y.append(ratio2)

            fig = plt.figure()
            plt.plot(X, Y, '.')
            plt.title(title)
            plt.xlabel("replicate %d" % (rep1))
            plt.ylabel("replicate %d" % (rep2))
            plt.show()
            plt.close()

check_replicates (rep_data1, title='Old sample')
check_replicates (rep_data2, title='New sample')

# get average value across the replicates

            

sys.exit(1)


    
    
# read data
sample_peptide_ptm_mratio1 = read_MS("Sangwoo_histone_ratios_092122.csv")
sample_peptide_ptm_ratio2 = read_MS("Sangwoo_histone_ratios_old_full.csv")
#sample_peptide_ptm_ratio2 = read_MS("Sangwoo_histone_ratios_old.csv")
#sample_peptide_ptm_ratio2 = read_MS("Sangwoo_histone_ratios_new.csv")
sample_peptide_ptm_mratio1['sup'] = copy.deepcopy(sample_peptide_ptm_mratio1['supernatant'])
del sample_peptide_ptm_mratio1['supernatant']

# reorganize data by replicates and get the average value
#sample_peptide_ptm_rep_ratio1, sample_peptide_ptm_mratio1 = reorganize_data(sample_peptide_ptm_ratio1)
sample_peptide_ptm_rep_ratio2, sample_peptide_ptm_mratio2 = reorganize_data(sample_peptide_ptm_ratio2)


# draw data
#draw_data (sample_peptide_ptm_mratio2)

#sys.exit(1)

# compare two data set
X, Y = [], []
for sample in sample_peptide_ptm_mratio1:
    for peptide in sample_peptide_ptm_mratio1[sample]:
        for ptm in sample_peptide_ptm_mratio1[sample][peptide]:
            try:
                mratio1 = sample_peptide_ptm_mratio1[sample][peptide][ptm]
                mratio2 = sample_peptide_ptm_mratio2[sample][peptide][ptm]
                X.append(mratio1)
                Y.append(mratio2)
            except:
                pass


fig =plt.figure()
plt.plot(X, Y, '.')
#plt.show()
plt.close()



histone_singleptm_ratios1 = aggregate_data(sample_peptide_ptm_mratio1)
histone_singleptm_ratios2 = aggregate_data(sample_peptide_ptm_mratio2)

X, Y = [], []
for histone in histone_singleptm_ratios1:
    for singleptm in histone_singleptm_ratios1[histone]:
        try:
            X += histone_singleptm_ratios1[histone][singleptm]
            Y += histone_singleptm_ratios2[histone][singleptm]
        except:
            continue

fig = plt.figure()
plt.plot(X, Y, '.')
#plt.show()
plt.close()

plot_hmap(histone_singleptm_ratios2)

sys.exit(1)            

label_sup1, label_pellet1 = get_score(sample_peptide_ptm_mratio1)
label_sup2, label_pellet2 = get_score(sample_peptide_ptm_mratio2)

X, Y = [], []
for label in label_sup1:
    try:
        sup1 = label_sup1[label]
        sup2 = label_sup2[label]
        X.append(sup1)
        Y.append(sup2)
    except:
        continue

fig = plt.figure()
plt.plot(X, Y, '.')
#plt.show()
plt.close()

                                                         




sys.exit(1)
            

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
                      yticklabels=1,
                      annot=True,
                      annot_kws={"size": 8},
                      cbar_pos=None)

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
plt.yticks(rotation=0)
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


# organize the data by replicate
sample_peptide_ptm_ratios = {}
rep_sample_peptide_ptm_ratio = {}
for sample in sample_peptide_ptm_ratio:
    peptide_ptm_ratio = sample_peptide_ptm_ratio[sample]
    cols = sample.strip().split('_')
    new_sample, rep = cols[-2], cols[-1]
    rep = int(rep)

    if rep not in rep_sample_peptide_ptm_ratio:
        rep_sample_peptide_ptm_ratio[rep] = {}
    rep_sample_peptide_ptm_ratio[rep][new_sample] = peptide_ptm_ratio

    if new_sample not in sample_peptide_ptm_ratios:
        sample_peptide_ptm_ratios[new_sample] = {}

"""
# compare replicates (old)
data_list1 = [[] for i in range(3)]
for sample in sample_peptide_ptm_rep_ratio1:
    for peptide in sample_peptide_ptm_rep_ratio1[sample]:
        for ptm in sample_peptide_ptm_rep_ratio1[sample][peptide]:
            ratio1 = sample_peptide_ptm_rep_ratio1[sample][peptide][ptm][1]
            ratio2 = sample_peptide_ptm_rep_ratio1[sample][peptide][ptm][2]
            ratio3 = sample_peptide_ptm_rep_ratio1[sample][peptide][ptm][3]
            data_list1[0].append(ratio1)
            data_list1[1].append(ratio2)
            data_list1[2].append(ratio3)

            
for i in range(len(data_list1)-1):
    for j in range(i+1, len(data_list1)):
        fig = plt.figure()
        plt.plot(data_list1[i], data_list1[j], '.')
        plt.xlabel("replicate %d" % (i+1))
        plt.ylabel("replicate %d" % (j+1))
        plt.title("Old data")
        #plt.show()
        plt.close()

# compare replicates
data_list2 = [[] for i in range(3)]
for sample in sample_peptide_ptm_rep_ratio2:
    for peptide in sample_peptide_ptm_rep_ratio2[sample]:
        for ptm in sample_peptide_ptm_rep_ratio2[sample][peptide]:
            ratio1 = sample_peptide_ptm_rep_ratio2[sample][peptide][ptm][1]
            ratio2 = sample_peptide_ptm_rep_ratio2[sample][peptide][ptm][2]
            ratio3 = sample_peptide_ptm_rep_ratio2[sample][peptide][ptm][3]
            data_list2[0].append(ratio1)
            data_list2[1].append(ratio2)
            data_list2[2].append(ratio3)

            
for i in range(len(data_list2)-1):
    for j in range(i+1, len(data_list2)):
        fig = plt.figure()
        plt.plot(data_list2[i], data_list2[j], '.')
        plt.xlabel("replicate %d" % (i+1))
        plt.ylabel("replicate %d" % (j+1))
        plt.title("New data")
        #plt.show()
        plt.close()
