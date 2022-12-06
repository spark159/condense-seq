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

for chr in chr_list:
    num_fname = path + '_'.join([cell, sample, agent, chr]) + '_num.cn'
    score_fname = path + '_'.join([cell, sample, agent, chr]) + '_score.cn'
    print "reading %s" % (chr)
    f = open(score_fname, 'w')
    s = ['SNP', 'Chromosome', 'PhysicalPosition']
    ID = 0
    First = True
    for line in open(num_fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        if First:
            s += cols[3:-2]
            print >> f, '\t'.join(s)
            First = False
            continue
        if cols[-1] == '*':
            continue
        _, chr, pos = cols[:3]
        counts = [int(count) for count in cols[3:]]
        if sum(counts) <= 0:
            continue
        control = counts[-1]
        if control <= 0:
            continue
        control +=1
        scores = []
        for count in counts[:-1]:
            count +=1
            score = -math.log(float(count)/control)
            scores.append(score)
        s = [str(ID), chr, str(pos)]
        s += [str(score) for score in scores]
        print >> f, '\t'.join(s)
        ID +=1
    f.close()
