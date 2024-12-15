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
import Interval_dict_python3 as Interval_dict
import copy

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
path = "/home/spark159/../../storage/"


# chromosome choice
chr_list = ['chr1']


# HMM type
HMM_type = 'Gaussian'
#covariance_type = 'tied'
covariance_type = 'spherical'


# state parameters
state_num = 5
seed = 0


# HMM fname
HMM_fname = "multiAD_%s_%d" % (HMM_type, state_num)


# bedgraph file information
name_fname = {"LAD":"./data/Lamin_score.bedgraph",
               "SPAD":"./data/SON_score.bedgraph",
               "NAD":"./data/K562_NucleolarDamID.bedgraph"}


# set data names to be analyzed
names = ['SPAD', 'NAD', 'LAD']


# set common bin_size
bin_size = 10000


# load data
name_chr_win_value = {}
for name in names:
    fname = name_fname[name]
    chr_win_value = read_bedgraph(path+fname)
    name_chr_win_value[name] = chr_win_value


# find common data range
chr_range = {}
for chr in chr_list:
    minpos_list = []
    maxpos_list = []
    for name in names:
        win_list = sorted(name_chr_win_value[name][chr].keys())
        minpos_list.append(win_list[0][0])
        maxpos_list.append(win_list[-1][1])
    minbinID = int(max(minpos_list) / bin_size)
    maxbinID = int(min(maxpos_list) / bin_size + 1)
    chr_range[chr] = (minbinID, maxbinID)
    

# binning the data
name_chr_binID_value = {}
name_values = {}
for chr in chr_list:
    minbinID, maxbinID = chr_range[chr]
    binID_interval = {i:(i*bin_size, (i+1)*bin_size) for i in range(minbinID, maxbinID)}
    bin_dict = Interval_dict.bin_hash(binID_interval, 
                                      bin_size=bin_size,
                                      bin_step=bin_size,
                                      max_pos=maxbinID*bin_size)

    for name in names:
        win_value = name_chr_win_value[name][chr]
        for win in sorted(win_value):
            start, end = win
            value = win_value[win]
            bin_dict.insert_range(start, end, value)
        binID_value = copy.deepcopy(bin_dict.get())
        bin_dict.ID_value = {} # clear out

        if name not in name_chr_binID_value:
            name_chr_binID_value[name] = {}
        name_chr_binID_value[name][chr] = binID_value

        if name not in name_values:
            name_values[name] = []
        name_values[name] += binID_value.values()


# extract the input
name_mean = {}
name_std = {}
for name in names:
    mean = np.mean(name_values[name])
    std = np.std(name_values[name])
    name_mean[name] = mean
    name_std[name] = std
del name_values
    
data = []
lengths = []
for chr in chr_list:
    minbinID, maxbinID = chr_range[chr]
    binID_list = range(minbinID, maxbinID)
    for i in range(len(binID_list)):
        binID = binID_list[i]
        row = []
        for name in names:
            try:
                value = name_chr_binID_value[name][chr][binID]
            except:
                value = 0.0
            mean = name_mean[name]
            std = name_std[name]
            zscore = float(value - mean)/std
            row.append(zscore)
        data.append(row)
    lengths.append(len(binID_list))

    
# HMM training
model = hmm.GaussianHMM(n_components=state_num,
                        covariance_type=covariance_type,
                        n_iter=1000,
                        verbose=True,
                        random_state=seed)
model.fit(data, lengths=lengths)

state_list = model.predict(data)

chr_binID_state = {}
pt = 0
for chr, length in zip(chr_list, lengths):
    minbinID, maxbinID = chr_range[chr]
    binID_list = range(minbinID, maxbinID)
    assert len(binID_list) == length
    for binID in binID_list:
        state = state_list[pt]
        if chr not in chr_binID_state:
            chr_binID_state[chr] = {}
        chr_binID_state[chr][binID] = state
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
fig, axes = plt.subplots(figsize=(5, 10),
                         nrows=state_num,
                         sharex=True)

for k in range(state_num):
    for i in range(len(names)):
        name = names[i]
        mean = model.means_[k][i]
        std = np.sqrt(model.covars_[k][i][i])
        x = np.linspace(mean - 3*std, mean + 3*std, 100)
        axes[k].plot(x, stats.norm.pdf(x, mean, std), alpha=0.3, label=name)

#leg = plt.legend(loc='upper right')
#for lh in leg.legendHandles:
#    lh.set_alpha(1)
plt.savefig(HMM_fname + '_em_dist.svg', format='svg', bbox_inches='tight')
plt.close()


# plot emission mean matrix
img = []
for k in range(state_num):
    img.append(model.means_[k])

#cmap = "Spectral_r"
cmap = 'bwr'
#cmap = 'jet'
height = 0.5*state_num
width = 0.5*len(names)
fig = plt.figure(figsize=(width, height))
plt.imshow(img, cmap=cmap, vmin=-2, vmax=2)
plt.xticks(range(len(names)), names,
           ha='right', va='center', rotation_mode='anchor', rotation=45)
plt.yticks(range(state_num), [str(i+1) for i in range(state_num)])
plt.ylabel("State #")
plt.title("Emission matrix")
#plt.colorbar(shrink=0.5)
plt.savefig(HMM_fname + '_em_mean.svg', format='svg', bbox_inches='tight')
#plt.show()
plt.close()

    
## write HMM annotation
def write_line (sinterval, f):
    state, chr, start, end = sinterval
    s = '%s\t%d\t%d\t%s' % (chr, start, end, 'E'+str(state+1))
    print (s, end='\n', file=f)

f = open(HMM_fname + '.bed', 'w')
prev_state = None
sinterval = []
for chr in chr_list:
    binID_state = chr_binID_state[chr]
    minbinID, maxbinID = chr_range[chr]
    for binID in range(minbinID, maxbinID):
        start, end = binID*bin_size, (binID+1)*bin_size
        state = binID_state[binID]
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
