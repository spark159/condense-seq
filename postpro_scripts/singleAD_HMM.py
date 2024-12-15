import matplotlib.pyplot as plt
import numpy as np
import math
import random
import matplotlib.cm as cm
import scipy
import matplotlib as mpl
from hmmlearn import hmm
import scipy.stats as stats
import pickle
from functools import cmp_to_key

def chr_cmp (chr1, chr2):
    try:
        chrnum1 = int(chr1[3:])
    except:
        chrnum1 = chr1[3:]
    try:
        chrnum2 = int(chr2[3:])
    except:
        chrnum2 = chr2[3:]
    if type(chrnum1) == type(chrnum2):
        if chrnum1 < chrnum2:
            return -1
        elif chrnum1 > chrnum2:
            return 1
        else:
            return 0
    else:
        if type(chrnum1) != int:
            return 1
        return -1
    

def read_bedgraph(fname):
    chr_win_value = {}
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        chr, start, end, value = cols
        start = int(start)
        end = int(end)
        value = float(value)
        if chr not in chr_win_value:
            chr_win_value[chr] = {}
        win = (start, end)
        chr_win_value[chr][win] = value
    return chr_win_value
    

# parameters
#path = "/home/spark159/../../storage/"
path = "./data/"

# HMM type
HMM_type = 'Gaussian'
#covariance_type = 'tied'
covariance_type = 'spherical'

# state parameters
state_num = 2
seed = 0

# bedgraph file information
name_fname = {"LAD":"Lamin_score.bedgraph",
              "SPAD":"SON_score.bedgraph",
              "NAD":"K562_NucleolarDamID.bedgraph"}


# set data names to be analyzed
#names = ['SPAD', 'NAD', 'LAD', 'eigen']
names = ['SPAD', 'NAD', 'LAD']

# set chromosome list
chr_list = None

# read data and HMM analysis
for name in names:
    print ("Processing %s" % (name))
    fname = name_fname[name]
    chr_win_value = read_bedgraph(path+fname)

    if chr_list == None:
        chr_list = sorted(chr_win_value.keys(), key=cmp_to_key(chr_cmp))

    # extract data
    data = []
    lengths = []
    for chr in chr_list:
        win_value = chr_win_value[chr]
        for win in sorted(win_value):
            start, end = win
            value = win_value[win]
            data.append([value])
        lengths.append(len(win_value))

    # data standardization
    mean = np.mean([value[0] for value in data])
    std = np.std([value[0] for value in data])
    for value in data:
        value[0] = float(value[0] - mean) / std

    # HMM training
    HMM_fname = "%s_%s_%d" % (name, HMM_type, state_num)

    model = hmm.GaussianHMM(n_components=state_num,
                            covariance_type=covariance_type,
                            #init_params='stc',
                            n_iter=1000,
                            verbose=True,
                            random_state=seed)

    #model.means_ = [[-2], [0], [2]]

    model.fit(data, lengths=lengths)
    state_list = model.predict(data)
    
    chr_win_state = {}
    pt = 0
    for chr, length in zip(chr_list, lengths):
        win_value = chr_win_value[chr]
        assert len(win_value) == length
        for win in sorted(win_value):
            state = state_list[pt]
            if chr not in chr_win_state:
                chr_win_state[chr] = {}
            chr_win_state[chr][win] = state
            pt +=1
    assert len(state_list) == pt

    # write model parameters
    pickle.dump(model.transmat_, open(HMM_fname + '_transm.pickle', "wb"), protocol=2)
    pickle.dump(model.means_, open(HMM_fname + '_means.pickle', "wb"), protocol=2)
    pickle.dump(model.covars_, open(HMM_fname + '_covars.pickle', "wb"), protocol=2)


    # plot transition matrix
    fig = plt.figure()
    plt.imshow(model.transmat_, vmin=0, vmax=1)
    plt.colorbar()
    plt.xlabel("State To")
    plt.ylabel("State From")
    plt.xticks(range(state_num), [str(i+1) for i in range(state_num)])
    plt.yticks(range(state_num), [str(i+1) for i in range(state_num)])
    plt.title("Transition prob matrix")
    plt.savefig(HMM_fname + '_transm.svg', format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()


    # plot emision distribution
    fig = plt.figure()
    plt.hist([value[0] for value in data], bins=100, alpha=0.25, density=True)

    for k in range(state_num):
        mean = model.means_[k][0]
        std = np.sqrt(model.covars_[k][0][0])
        x = np.linspace(mean - 3*std, mean + 3*std, 100)
        plt.plot(x, stats.norm.pdf(x, mean, std), alpha=1, label='E'+str(k+1))
        #plt.gca().fill_between(x, stats.norm.pdf(x, mean, std), alpha=0.3, label='E'+str(k+1))

    leg = plt.legend(loc='upper right')
    for lh in leg.legend_handles:
        lh.set_alpha(1)

    plt.savefig(HMM_fname + '_em_dist.svg', format='svg', bbox_inches='tight')
    plt.close()


    ## add state information to the input bedgraph
    f = open(fname.rsplit('.')[0] + '_state.txt', 'w')
    for chr in chr_list:
        for win in sorted(chr_win_value[chr]):
            start, end = win
            value = chr_win_value[chr][win]
            state = chr_win_state[chr][win]
            s = '%s\t%d\t%d\t%f\t%s' % (chr, start, end, value, 'E'+str(state+1))
            print (s, end='\n', file=f)
    f.close()
    

    ## write HMM annotation
    def write_line (sinterval, f):
        state, chr, start, end = sinterval
        s = '%s\t%d\t%d\t%s' % (chr, start, end, 'E'+str(state+1))
        print (s, end='\n', file=f)

    f = open(HMM_fname + '.bed', 'w')
    prev_state = None
    sinterval = []
    for chr in chr_list:
        win_state = chr_win_state[chr]
        for win in sorted(win_state):
            start, end = win
            state = win_state[win]
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
    f.close()

