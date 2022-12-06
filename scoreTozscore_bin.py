import math
import numpy as np

# parameters
path = "/home/spark159/../../media/spark159/sw/"

# sample information
#cell = 'H1'
#cell = 'GM'
cell = 'mCD8T'
#sample = 'NCP'
#sample = 'WT-NCP'
#sample = 'inht-NCP'
sample = 'KO-NCP'
agent = 'sp'

# bin information
bin_size = 1000

# set species and gender
#species = 'human'
#gender = 'male'
if cell in ['H1', 'GM']:
    species = 'human'
elif cell in ['mCD8T']:
    species = 'mouse'

if cell in ['H1']:
    gender = 'male'
elif cell in ['GM', 'mCD8T']:
    gender = 'female'

# set chromosome list
if species == 'human':
    chr_list = ['chr' + str(i) for i in range(1, 23)]
elif species == 'mouse':
    chr_list = ['chr' + str(i) for i in range(1, 20)]
chr_list += ['chrX']

if gender == 'male':
    chr_list += ['chrY']


# read all score files
names = []
name_data = {}
print "reading bin_score file"
score_fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb']) + '_score.cn'
First = True
for line in open(score_fname):
    line = line.strip()
    if not line:
        continue
    cols = line.split('\t')
    if First:
        names = cols[4:]
        First = False
        continue
    scores = [float(score) for score in cols[4:]]
    for name, score in zip(names, scores):
        if name not in name_data:
            name_data[name] = []
        name_data[name].append(score)

# get zscore
print "computing zscore"
for name, data in name_data.items():
    mean = np.mean(data)
    std = np.std(data)
    name_data[name] = [float(value-mean)/std for value in data]

# write zscore files
idx = 0
print "writing bin_zscore file"
zscore_fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb']) + '_zscore.cn'
f = open(zscore_fname, 'w')
s = ['SNP', 'Chromosome', 'Start', 'End']
First = True
for line in open(score_fname):
    line = line.strip()
    if not line:
        continue
    cols = line.split('\t')
    if First:
        s += cols[4:]
        print >> f, '\t'.join(s)
        First = False
        continue
    s = cols[:4]
    for name in names:
        s.append(str(name_data[name][idx]))
    print >> f, '\t'.join(s)
    idx +=1
f.close()

# sanity check
for name in names:
    assert len(name_data[name]) == idx
