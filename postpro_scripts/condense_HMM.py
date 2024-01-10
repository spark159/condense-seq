import matplotlib.pyplot as plt
import numpy as np
import math
import random
import matplotlib.cm as cm
import scipy
import seaborn as sns
import matplotlib as mpl
from hmmlearn import hmm
import scipy.stats as stats
import pickle

def rescale (value_list, old_st, old_ed, new_st, new_ed):
    output = []
    for value in value_list:
        assert value >= old_st and value <= old_ed
        new_value = new_st + (new_ed - new_st)*float(value-old_st)/(old_ed-old_st)
        output.append(new_value)
    return output


def read_bin_file (fname, chr_choices=None):
    First = True
    ID_pos = {}
    name_ID_value = {}
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        if First:
            names = []
            for name in cols[4:]:
                #agent = name.rsplit('.', 1)[0].split('-')[-2]
                #tnum = int(name.rsplit('.', 1)[0].split('-')[-1])
                #names.append(tnum)
                names.append(name.rsplit('.', 1)[0])
            First = False
            continue
        ID, chr, start, end = cols[:4]
        start, end = int(start), int(end)
        if chr_choices!=None and chr not in chr_choices:
            continue
        pos = int(float(start + end)/2)
        ID = (chr, start, end)
        assert ID not in ID_pos
        ID_pos[ID] = pos
        values = [float(value) for value in cols[4:]]
        for name, value in zip(names, values):
            if name not in name_ID_value:
                name_ID_value[name] = {}
            assert ID not in name_ID_value[name]
            name_ID_value[name][ID] = value
    return ID_pos, name_ID_value

def read_bin_Chalf (fname, chr_choices=None):
    ID_pos = {}
    ID_value = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        if First:
            First = False
            continue
        binID, chr, st, ed, Chalf = cols
        if chr_choices != None and chr not in chr_choices:
            continue
        st, ed = int(st), int(ed)
        ID = (st, ed)
        st, ed = int(st), int(ed)
        Chalf = float(Chalf)
        pos = int(float(st + ed)/2)
        ID_pos[ID] = pos
        ID_value[ID] = Chalf
    return ID_pos, ID_value

# parameters
# set path
#path = ""
#path = "/home/spark159/../../media/spark159/sw/"
path = "/home/spark159/../../storage/"

agent_fullname = {'sp':'Spermine(4+)',
                  'spd':'Spermidine(3+)',
                  'CoH':'Cobalt Hexammine(3+)',
                  'PEG':'PEG 8000',
                  'HP1a':'HP1$\\alpha$',
                  'HP1bSUV':'HP1$\\beta$+tSUV',
                  'LKH':'Linker histone1',
                  'Ki67':'Ki67',
                  'FUS':'FUS',
                  'Mg':'Magnesium',
                  'Ca':'Calcium'}

agent_charge = {'sp':'4+',
                'spd':'3+',
                'CoH':'3+',
                'PEG':'',
                'Ca':'2+',
                'Mg':'2+',
                'HP1a':'',
                'HP1bSUV':'',
                'LKH':'',
                'Ki67':'',
                'FUS':''}

# experiment list (cell, sample, agent, tnum)
#exp_list = [('H1', 'NCP', 'sp', 8),
#            ('H1', 'NCP', 'spd', 6),
#            ('H1', 'NCP', 'CoH', 5),
#            ('H1', 'NCP', 'PEG', 6),
#            ('H1', 'NCP', 'Ca', 5),
#            ('H1', 'NCP', 'Mg', 5),
#            ('H1', 'NCP', 'HP1a', 3),
#            ('H1', 'NCP', 'HP1bSUV', 4),
#            ('H1', 'NCP', 'LKH', 3),
#            ('H1', 'NCP', 'Ki67', 4),
#            ('H1', 'NCP', 'FUS', 4)]

#exp_list = [('H1', 'NCP', 'sp', 8),
#            ('H1', 'NCP', 'spd', 6),
#            ('H1', 'NCP', 'CoH', 5),
#            ('H1', 'NCP', 'PEG', 6),
#            ('H1', 'NCP', 'Ca', 5),
#            ('H1', 'NCP', 'HP1a', 3),
#            ('H1', 'NCP', 'LKH', 3),
#            ('H1', 'NCP', 'Ki67', 4)]

exp_list = [('H1', 'NCP', 'sp', 8),
            ('H1', 'NCP', 'HP1a', 3),
            ('H1', 'NCP', 'LKH', 3),
            ('H1', 'NCP', 'Ki67', 4)]

# bin size
#bin_size = 10000
#bin_size = 5000
bin_size = 1000

# chromosome choice
chr_list = ['chr1']

# HMM type
HMM_type = 'Gaussian'
covariance_type = 'spherical'

#HMM_type = 'Poisson'

#HMM_type = 'Categorical'
#cate_num = 10

# HMM optimization
HMM_opt = False

# state parameters
#state_num = 12
#state_num = 11
state_num = 10
#state_num = 6
seed = 0

# data type
if HMM_type == 'Gaussian':
    dtype = 'zscore'
elif HMM_type == 'Poisson':
    dtype = 'score'
elif HMM_type == 'Categorical':
    dtype = 'score'

# HMM file name
HMM_fname = '_'.join([str(state_num), 'ConHMMstate', HMM_type, str(int(bin_size/1000.0))+'kb'])


# read data
exp_ID_pos = {}
exp_ID_score = {}
for exp in exp_list:
    cell, sample, agent, tnum = exp

    fname = '_'.join([cell, sample, agent,
                      str(int(bin_size/1000.0)) + 'kb',
                      dtype]) + '.cn'

    if dtype in ['score', 'zscore']:
        field_name = '-'.join([cell, sample, agent, str(tnum)])
        
    elif dtype in ['Chalf', 'zChalf']:
        field_name = 'Chalf'
        
    ID_pos, field_ID_value = read_bin_file(path + fname, chr_choices=chr_list) 
    ID_score = field_ID_value[field_name]
        
    exp_ID_pos[exp] = ID_pos
    exp_ID_score[exp] = ID_score


# get common IDs
ID_list = set([])
for i in range(len(exp_list)):
    exp = exp_list[i]
    if i == 0:
        ID_list |= set(exp_ID_pos[exp].keys())
        continue
    ID_list &= set(exp_ID_pos[exp].keys())
ID_list = sorted(list(ID_list))


# extract input
if HMM_type == 'Gaussian':
    X = []
    for ID in ID_list:
        row = []
        for exp in exp_list:
            score = exp_ID_score[exp][ID]
            row.append(score)
        X.append(row)

elif HMM_type == 'Poisson':
    exp_trial = {}
    for exp in exp_list:
        max_score = max([exp_ID_score[exp][ID] for ID in ID_list])
        min_prob = np.exp(-max_score)
        trial = int(1.0/min_prob)
        exp_trial[exp] = trial
    
    X = []
    for ID in ID_list:
        row = []
        for exp in exp_list:
            score = exp_ID_score[exp][ID]
            count = int(np.exp(-score)*exp_trial[exp])
            row.append(count)
        X.append(row)

elif HMM_type == 'Categorical':
    X = [[] for ID in ID_list]
    for exp in exp_list:
        min_score = min(exp_ID_score[exp].values())
        max_score = max(exp_ID_score[exp].values())
        bins = np.linspace(min_score-1, max_score+1, num=cate_num+1)
        for i in range(len(ID_list)):
            ID = ID_list[i]
            score = exp_ID_score[exp][ID]
            Find = False
            for k in range(cate_num):
                left = bins[k]
                right = bins[k+1]
                if score >= left and score < right:
                    X[i].append(k)
                    Find = True
                    break
            assert Find

# HMM score optimization
if HMM_opt:
    pass # Todo
    

# HMM training
if HMM_type == 'Gaussian':
    model = hmm.GaussianHMM(n_components=state_num,
                            covariance_type=covariance_type,
                            n_iter=1000,
                            verbose=True,
                            random_state=seed)
    model.fit(X)

    # write model parameters
    pickle.dump(model.transmat_, open(HMM_fname + '_transm.pickle', "wb"), protocol=2)
    pickle.dump(model.means_, open(HMM_fname + '_means.pickle', "wb"), protocol=2)
    pickle.dump(model.covars_, open(HMM_fname + '_covars.pickle', "wb"), protocol=2)


elif HMM_type == 'Poisson':
    model = hmm.PoissonHMM(n_components=state_num,
                           n_iter=1000,
                           verbose=True,
                           random_state=seed)
    model.fit(X)

elif HMM_type == 'Categorical':
    model = hmm.CategoricalHMM(n_components=state_num,
                               n_iter=1000,
                               verbose=True,
                               random_state=seed)
    model.fit(X)



## write HMM annotation
def write_line (sinterval, f):
    state, chr, start, end = sinterval
    s = '%s\t%d\t%d\t%s' % (chr, start, end, 'E'+str(state+1))
    print (s, end='\n', file=f)

f = open(HMM_fname + '.bed', 'w')

state_list = model.predict(X)
prev_state = None
sinterval = []
for ID, state in zip(ID_list, state_list):
    chr, start, end = ID
    if prev_state == None:
        sinterval = [state, chr, start]
        prev_state = state
        prev_chr = chr
        prev_end = end
    if prev_state != state or prev_chr != chr:
        sinterval.append(prev_end)
        write_line (sinterval, f)
        sinterval = [state, chr, start]
    prev_state = state
    prev_chr = chr
    prev_end = end
sinterval.append(prev_end)
write_line (sinterval, f)



# plot transition matrix
fig = plt.figure()
plt.imshow(model.transmat_)
plt.colorbar()
plt.xlabel("State To")
plt.ylabel("State From")
plt.xticks(range(state_num), [str(i+1) for i in range(state_num)])
plt.yticks(range(state_num), [str(i+1) for i in range(state_num)])
plt.title("Transition prob matrix")
plt.savefig("transm.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()


if HMM_type == 'Gaussian':
    # plot emision distribution
    if covariance_type in ['spherical', 'diag']:
        fig, axes = plt.subplots(figsize=(5, 10),
                                 nrows=state_num,
                                 sharex=True)

        for k in range(state_num):
            for i in range(len(exp_list)):
                cell, sample, agent, tnum = exp_list[i]
                mean = model.means_[k][i]
                std = np.sqrt(model.covars_[k][i][i])
                x = np.linspace(mean - 3*std, mean + 3*std, 100)
                axes[k].fill_between(x, stats.norm.pdf(x, mean, std),
                                     #color=agent_color[agent],
                                     alpha=0.3,
                                     label=agent)

        #leg = plt.legend(loc='upper right')
        #for lh in leg.legendHandles:
        #    lh.set_alpha(1)
        plt.savefig("emission_dist.svg", format='svg', bbox_inches='tight')
        plt.close()


    # plot emission mean matrix
    img = []
    for k in range(state_num):
        img.append(model.means_[k])

    cmap = 'bwr_r'
    #cmap = 'jet'
    fig = plt.figure()
    plt.imshow(img, cmap=cmap, vmin=-2, vmax=2)
    plt.xticks(range(len(exp_list)), [agent for cell, sample, agent, tnum in exp_list],
               ha='right', va='center', rotation_mode='anchor', rotation=45)
    plt.yticks(range(state_num), [str(i+1) for i in range(state_num)])
    plt.ylabel("State #")
    plt.title("Emission matrix")
    plt.colorbar(shrink=0.5)
    plt.savefig("emission_mean.svg", format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()


if HMM_type == 'Poisson':
    # plot emission lambda matrix
    exp_state_score = {}
    for k in range(state_num):
        for i in range(len(exp_list)):
            exp = exp_list[i]
            lamb = model.lambdas_[k][i]
            score = -np.log(float(lamb)/exp_trial[exp])
            if exp not in exp_state_score:
                exp_state_score[exp] = {}
            exp_state_score[exp][k] = score

    img = []
    for k in range(state_num):
        row = []
        for i in range(len(exp_list)):
            exp = exp_list[i]
            score = exp_state_score[exp][k]
            min_score = min(exp_state_score[exp].values())
            max_score = max(exp_state_score[exp].values())
            re_score = float(score - min_score)/max_score
            row.append(re_score)
        img.append(row)

    cmap = 'Blues'

    fig = plt.figure()
    plt.imshow(img, cmap=cmap)
    plt.xticks(range(len(exp_list)), [agent for cell, sample, agent, tnum in exp_list],
               ha='right', va='center', rotation_mode='anchor', rotation=45)
    plt.yticks(range(state_num), [str(i+1) for i in range(state_num)])
    plt.ylabel("State #")
    plt.title("Emission matrix")
    plt.colorbar(shrink=0.5)
    plt.savefig("emission_mean.svg", format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()
