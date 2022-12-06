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


new_fname = path + '_'.join([cell, sample, agent]) + '_score.cn'
f = open(new_fname, 'w')

ID = 0
for i in range(len(chr_list)):
    print "reading %s" % (chr_list[i])
    score_fname = path + '_'.join([cell, sample, agent, chr_list[i]]) + '_score.cn'
    First = True
    for line in open(score_fname):
        line = line.strip()
        if not line:
            continue
        if First:
            if i == 0:
                print >> f, line
            First = False
            continue
        cols = line.split('\t')
        cols[0] = str(ID)
        print >> f, '\t'.join(cols)
        ID +=1

f.close()
