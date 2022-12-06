import math    

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

num_fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb']) + '_num.cn'
score_fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb']) + '_score.cn'

print "reading bin_num file"
chr_range_scores = {}
field_names = []
First = True
for line in open(num_fname):
    line = line.strip()
    if not line:
        continue
    cols = line.split('\t')
    if First:
        field_names += cols[:-2]
        First = False
        continue
    if cols[-1] == '*':
        continue
    ID, chr, st, ed = cols[:4]
    counts = [int(count) for count in cols[4:]]
    if 0 in counts:
        continue
    control = counts[-1]
    scores = []
    for count in counts[:-1]:
        score = -math.log(float(count)/control)
        scores.append(score)

    if chr not in chr_range_scores:
        chr_range_scores[chr] = {}
    chr_range_scores[chr][(st, ed)] = scores

print "writing bin_score file"
f = open(score_fname, 'w')
print >> f, '\t'.join(field_names)

ID = 0
for chr in chr_list:
    for st, ed in sorted(chr_range_scores[chr]):
        s = [str(ID), chr, st, ed]
        s += [str(score) for score in chr_range_scores[chr][(st,ed)]]
        print >> f, '\t'.join(s)
        ID +=1
f.close()
