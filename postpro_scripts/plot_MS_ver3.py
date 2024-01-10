import sys
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import copy
import scipy

def weighted_mean (values, weights):
    assert len(values) == len(weights)
    wmean = 0.0
    for value, weight in zip(values, weights):
        wmean += value*weight
    wmean /= float(sum(weights))
    return wmean

def rescale (value, old_st, old_ed, new_st, new_ed):
    assert value >= old_st and value <= old_ed
    new_value = new_st + (new_ed - new_st)*float(value-old_st)/(old_ed-old_st)
    return new_value

# two independent sample t-test
# equal sample size and variance
def get_two_sample_t (data1, data2):
    assert len(data1) == len(data2)
    assert len(data1) > 1
    mean1 = np.mean(data1)
    mean2 = np.mean(data2)
    std1 = np.std(data1)
    std2 = np.std(data2)
    t = float(mean1 - mean2)/np.sqrt((std1**2 + std2**2)/len(data1))
    return t

# two independent sample t-test
# unequal sample size and similar variance
def get_two_sample_t2 (data1, data2):
    size1 = len(data1)
    size2 = len(data2)
    assert size1 > 1
    assert size2 > 1
    mean1 = np.mean(data1)
    mean2 = np.mean(data2)
    std1 = np.std(data1)
    std2 = np.std(data2)
    std = np.sqrt(float((size1-1)*(std1**2) + (size2-1)*(std2**2)) / (size1+size2 - 2))
    t = float(mean1 - mean2)/std
    return t

# get pvalue (t-distribution and null hypothesis mean = 0)
def get_pvalue (t, deg):
    pval = scipy.stats.t.sf(np.abs(t), deg)
    return pval


def read_MS_table (fname):
    sample_peptide_ptm_signal = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        if First:
            fields = line.split('\t')
            fields = [field.strip('"') for field in fields]
            First = False
            continue
        if line.startswith('Peptide'):
            col_choices = []
            cols = line.split('\t')
            assert len(cols) == len(fields) + 1
            for i in range(len(cols)):
                if cols[i].strip() == 'Area':
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
            signal = cols[idx]
            #print peptide, ptm, ratio
            signal = float(signal)
            sample = fields[idx-1]
            if sample not in sample_peptide_ptm_signal:
                sample_peptide_ptm_signal[sample] = {}
            if peptide not in sample_peptide_ptm_signal[sample]:
                sample_peptide_ptm_signal[sample][peptide] = {}
            #print ptm
            assert ptm not in sample_peptide_ptm_signal[sample][peptide]
            sample_peptide_ptm_signal[sample][peptide][ptm] = signal

    return sample_peptide_ptm_signal

def parse_peptide (peptide):
    seq, histone_info = peptide[:-1].split('(')
    histone, st, ed = histone_info.split('_')
    st, ed = int(st), int(ed)
    return histone, st, ed, seq

def parse_ptm (ptm_string):
    pattern = '([A-Z][0-9]+(?:ac|me|ph|ub)[1-3]?)'
    ptm_list = re.findall(pattern, ptm_string)
    return ptm_list


# load MS data table
path = './'
fname = 'Garcialab_PTM_result_summary.txt'

sample_peptide_ptm_signal = read_MS_table(path+fname)




# pre-processing #1:
# filter out data in the case of input + sup <=0
# convert the signal to relative fraction for each peptide
# add small dummy number to fractions to eliminate zero

cdts = ['input', 'pellet', 'sup']
reps = range(4)
histones = ['H3', 'H4']

small = 10**(-10)

peptide_ptm_signal = sample_peptide_ptm_signal['input-0']

sample_peptide_ptm_frac = {}
for rep in reps:
    for peptide in peptide_ptm_signal:

        histone, st, ed, seq = parse_peptide(peptide)

        if histone not in histones:
            continue

        input_ptm_sig = sample_peptide_ptm_signal['input-%d' % (rep)][peptide]
        pellet_ptm_sig = sample_peptide_ptm_signal['pellet-%d' % (rep)][peptide]
        sup_ptm_sig = sample_peptide_ptm_signal['sup-%d' % (rep)][peptide]
        
        input_total = 0.0
        pellet_total = 0.0
        sup_total = 0.0

        input_ptm_frac = {}
        pellet_ptm_frac = {}
        sup_ptm_frac = {}

        for ptm in input_ptm_sig:
            input_sig = input_ptm_sig[ptm]
            pellet_sig = pellet_ptm_sig[ptm]
            sup_sig = sup_ptm_sig[ptm]

            #if input_sig <= 0:
            #    continue
            #if sup_sig <= 0:
            #    continue
            #if pellet_sig <= 0:
            #    continue
            
            if input_sig + sup_sig + pellet_sig <=0:
                continue

            #input_sig +=10000
            #pellet_sig +=10000
            #sup_sig +=10000
            
            input_total += input_sig
            pellet_total += pellet_sig
            sup_total += sup_sig

            input_ptm_frac[ptm] = input_sig
            pellet_ptm_frac[ptm] = pellet_sig
            sup_ptm_frac[ptm] = sup_sig

        for ptm in input_ptm_frac:
            input_ptm_frac[ptm] /= float(input_total)
            pellet_ptm_frac[ptm] /= float(pellet_total)
            sup_ptm_frac[ptm] /= float(input_total)

            input_ptm_frac[ptm] += small
            pellet_ptm_frac[ptm] += small
            sup_ptm_frac[ptm] += small

        if 'input-%d' % (rep) not in sample_peptide_ptm_frac:
            sample_peptide_ptm_frac['input-%d' % (rep)] = {}
        sample_peptide_ptm_frac['input-%d' % (rep)][peptide] = input_ptm_frac
    
        if 'pellet-%d' % (rep) not in sample_peptide_ptm_frac:
            sample_peptide_ptm_frac['pellet-%d' % (rep)] = {}
        sample_peptide_ptm_frac['pellet-%d' % (rep)][peptide] = pellet_ptm_frac
    
        if 'sup-%d' % (rep) not in sample_peptide_ptm_frac:
            sample_peptide_ptm_frac['sup-%d' % (rep)] = {}
        sample_peptide_ptm_frac['sup-%d' % (rep)][peptide] = sup_ptm_frac
        
    
# pre-processing #2:
# get fold-change compared to input

cdts = ['pellet', 'sup']
reps = range(4)
histones = ['H3', 'H4']

sample_peptide_ptm_fold = {}

peptides = sample_peptide_ptm_frac['input-0'].keys()

for rep in reps:
    for peptide in peptides:
        histone, st, ed, seq = parse_peptide(peptide)

        if histone not in histones:
            continue

        input_ptm_frac = sample_peptide_ptm_frac['input-%d' % (rep)][peptide]
        pellet_ptm_frac = sample_peptide_ptm_frac['pellet-%d' % (rep)][peptide]
        sup_ptm_frac = sample_peptide_ptm_frac['sup-%d' % (rep)][peptide]

        pellet_ptm_fold = {}
        sup_ptm_fold = {}
        for ptm in input_ptm_frac:
            input_frac = input_ptm_frac[ptm]
            pellet_frac = pellet_ptm_frac[ptm]
            sup_frac = sup_ptm_frac[ptm]

            pellet_fold = np.log2(pellet_frac/input_frac)
            sup_fold = np.log2(sup_frac/input_frac)

            pellet_ptm_fold[ptm] = pellet_fold
            sup_ptm_fold[ptm] = sup_fold

        if 'pellet-%d' % (rep) not in sample_peptide_ptm_fold:
            sample_peptide_ptm_fold['pellet-%d' % (rep)] = {}
        sample_peptide_ptm_fold['pellet-%d' % (rep)][peptide] = pellet_ptm_fold
    
        if 'sup-%d' % (rep) not in sample_peptide_ptm_fold:
            sample_peptide_ptm_fold['sup-%d' % (rep)] = {}
        sample_peptide_ptm_fold['sup-%d' % (rep)][peptide] = sup_ptm_fold


# post-processing #1:
# get mean fold-change difference compared to reference
# get p-value as significant difference compared to reference

cdts = ['sup']
reps = range(4)
histones = ['H3', 'H4']

cdt_peptide_ptm_delta = {}
cdt_peptide_ptm_pvalue = {}

peptides = sample_peptide_ptm_fold['sup-0'].keys()

for cdt in cdts:
    for peptide in peptides:
        histone, st, ed, seq = parse_peptide(peptide)

        if histone not in histones:
            continue

        # combine all replicate data
        ptm_fracs, ptm_folds = {}, {}
        for rep in reps:
            input_sample = 'input-%d' % (rep)
            cdt_sample = '%s-%d' % (cdt, rep)
            input_ptm_frac = sample_peptide_ptm_frac[input_sample][peptide]
            cdt_ptm_fold = sample_peptide_ptm_fold[cdt_sample][peptide]

            for ptm in input_ptm_frac:
                frac = input_ptm_frac[ptm]
                fold = cdt_ptm_fold[ptm]

                if ptm not in ptm_fracs:
                    ptm_fracs[ptm] = []
                ptm_fracs[ptm].append(frac)

                if ptm not in ptm_folds:
                    ptm_folds[ptm] = []
                ptm_folds[ptm].append(fold)

        # take the reference value first
        ref_ptm = histone_aarange + ' unmod'
        ref_fracs = ptm_fracs[ref_ptm]
        ref_folds = ptm_folds[ref_ptm]

        assert len(ref_fracs) == len(ref_folds)

        if len(ref_fracs) <= 1:
            continue

        mean_ref_fold = np.mean(ref_folds) # uniform mean
        wmean_ref_fold = weighted_mean(ref_folds, ref_fracs) # weighted mean
        
        # compare each ptm to the reference data
        ptms = list(set(ptm_folds.keys()) - set([ref_ptm]))        
        for ptm in ptms:
            target_fracs = ptm_fracs[ptm]
            target_folds = ptm_folds[ptm]

            assert len(target_fracs) == len(target_folds)

            if len(target_fracs) <= 1:
                continue
            
            mean_target_fold = np.mean(target_folds) # uniform mean
            wmean_target_fold = weighted_mean(target_folds, target_fracs) # weighted mean
            delta = mean_target_fold - mean_ref_fold
            #delta = wmean_target_fold - wmean_ref_fold
            
            t = get_two_sample_t2 (target_folds, ref_folds)
            pvalue = get_pvalue(t, len(target_folds) + len(ref_folds) - 2)

            if cdt not in cdt_peptide_ptm_delta:
                cdt_peptide_ptm_delta[cdt] = {}
            if peptide not in cdt_peptide_ptm_delta[cdt]:
                cdt_peptide_ptm_delta[cdt][peptide] = {}
            cdt_peptide_ptm_delta[cdt][peptide][ptm] = delta

            if cdt not in cdt_peptide_ptm_pvalue:
                cdt_peptide_ptm_pvalue[cdt] = {}
            if peptide not in cdt_peptide_ptm_pvalue[cdt]:
                cdt_peptide_ptm_pvalue[cdt][peptide] = {}
            cdt_peptide_ptm_pvalue[cdt][peptide][ptm] = pvalue



            
# draw delta and p-value in circle plot
# color: up (red) or down (blue) fold change compared to control
# size: significance proportional to -log10(p-value)
for cdt in cdt_peptide_ptm_delta:
    mark_delta = {}
    mark_pvalue = {}

    peptide_ptm_delta = cdt_peptide_ptm_delta[cdt]
    peptide_ptm_pvalue = cdt_peptide_ptm_pvalue[cdt]

    # find max and min values
    deltas, sigs = [], []
    for peptide in peptide_ptm_delta:
        for ptm in peptide_ptm_delta[peptide]:
            delta = peptide_ptm_delta[peptide][ptm]
            sig = -np.log10(peptide_ptm_pvalue[peptide][ptm])
            deltas.append(delta)
            sigs.append(sig)

    min_delta, max_delta = min(deltas), max(deltas)
    min_sig, max_sig = min(sigs), max(sigs)

    up_cmap = mpl.cm.get_cmap("Reds")
    down_cmap = mpl.cm.get_cmap("Blues_r")
    
    for peptide in cdt_peptide_ptm_delta[cdt]:
        ptm_delta = cdt_peptide_ptm_delta[cdt][peptide]
        
        # order by delta
        delta_ptm = sorted([(delta, ptm) for ptm, delta in ptm_delta.items()])
        ptms = [ptm for _, ptm in delta_ptm]
        colors, sizes = [], []
        labels = []
        for ptm in ptms:
            delta = cdt_peptide_ptm_delta[cdt][peptide][ptm]
            sig = -np.log10(cdt_peptide_ptm_pvalue[cdt][peptide][ptm])
            if delta >=0:
                color = up_cmap(rescale(delta, 0.0, max_delta, 0.0, 1.0))
            else:
                color = down_cmap(rescale(delta, min_delta, 0.0, 0.0, 1.0))
            size = rescale(sig, min_sig, max_sig, 5, 30)
            colors.append(color)
            sizes.append(size)

            his_range, marks = ptm.split(' ')
            histone = his_range.split('_')[0]
            label = ' '.join([histone, marks])
            labels.append(label)


        row_count = len(ptms)
        col_count = 1
        cell_width, cell_height = 1, 1 #inch
        marg_top = 0.1
        marg_bott = 0.1
        marg_left = 2
        marg_right = 0.1
        
        fig_width = cell_width*col_count + marg_left + marg_right
        fig_height = cell_height*row_count + marg_top + marg_bott

        fig = plt.figure(figsize=(fig_width, fig_height))
        fig.subplots_adjust(bottom=marg_bott/fig_height,
                            top=1.-marg_top/fig_height,
                            left=marg_left/fig_width,
                            right=1.-marg_right/fig_width)

        #fig = plt.figure()
        
        for i in range(len(ptms)):
            ptm = ptms[i]
            color = colors[i]
            size = sizes[i]

            plt.plot([0],[i],
                     marker='o',
                     markersize=size,
                     mfc=color,
                     markeredgewidth=0.5,
                     mec='k',
                     clip_on=False)

        #plt.xlim([-1, 1])
        #plt.ylim([-50, 50*len(ptms)])
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(False)
        plt.gca().spines['left'].set_visible(True)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().tick_params(left=False)
        plt.xticks([], [])
        plt.yticks(range(len(labels)), labels)
        plt.savefig('%s_delta.png' % (peptide),
                    bbox_inches='tight')
        plt.close()
