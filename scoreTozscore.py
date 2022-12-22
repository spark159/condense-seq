import math
import numpy as np

# parameters
#path = "/home/spark159/../../media/spark159/sw/"
path = "/home/spark159/../../storage/"

# sample information
#exp_list = [('H1', 'NCP', 'HP1a'),
#            ('H1', 'NCP', 'Ki67'),
#            ('H1', 'NCP', 'LKH')]

exp_list = [('H1', 'NCP', 'sp'),
            ('H1', 'NCP', 'HP1a')]


for cell, sample, agent in exp_list:

    print "Processing %s %s %s" % (cell, sample, agent)

    # set species and gender
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
    for chr in chr_list:
        print "reading %s" % (chr)
        score_fname = path + '_'.join([cell, sample, agent, chr]) + '_score.cn'
        First = True
        for line in open(score_fname):
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            if First:
                names = cols[3:]
                First = False
                continue
            scores = [float(score) for score in cols[3:]]
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
    for chr in chr_list:
        print "writing %s" % (chr)
        score_fname = path + '_'.join([cell, sample, agent, chr]) + '_score.cn'
        zscore_fname = path + '_'.join([cell, sample, agent, chr]) + '_zscore.cn'
        f = open(zscore_fname, 'w')
        s = ['SNP', 'Chromosome', 'PhysicalPosition']
        First = True
        for line in open(score_fname):
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            if First:
                s += cols[3:]
                print >> f, '\t'.join(s)
                First = False
                continue
            s = cols[:3]
            for name in names:
                s.append(str(name_data[name][idx]))
            print >> f, '\t'.join(s)
            idx +=1
        f.close()

    # sanity check
    for name in names:
        assert len(name_data[name]) == idx
    print
