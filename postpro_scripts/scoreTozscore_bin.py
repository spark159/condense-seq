import math
import numpy as np

# parameters
#path = "/home/spark159/../../media/spark159/sw/"
path = "/home/spark159/../../storage/replicates/"
#path = "/home/spark159/../../storage/"


# sample information

# experiment list (cell, sample, agent)
#exp_list = [('H1', 'NCP', 'sp'),
#            ('H1', 'NCP', 'spd'),
#            ('H1', 'NCP', 'CoH'),
#            ('H1', 'NCP', 'PEG'),
#            ('H1', 'NCP', 'Ca'),
#            ('H1', 'NCP', 'Mg'),
#            ('H1', 'NCP', 'HP1a'),
#            ('H1', 'NCP', 'HP1bSUV'),
#            ('H1', 'NCP', 'LKH'),
#            ('H1', 'NCP', 'Ki67'),
#            ('H1', 'NCP', 'FUS')]

#exp_list = [('H1', 'NCP', 'sp'),
#            ('H1', 'NCP', 'HP1a'),
#            ('H1', 'NCP', 'LKH'),
#            ('H1', 'NCP', 'Ki67')]

exp_list = [('H1', 'NCP', 'sp'),
            ('H1', 'NCP', 'spd'),
            ('H1', 'NCP', 'CoH'),
            ('H1', 'NCP', 'PEG'),
            ('H1', 'NCP', 'Ca'),
            ('H1', 'NCP', 'HP1a'),
            ('H1', 'NCP', 'LKH'),
            ('H1', 'NCP', 'Ki67'),
            ('H1', 'DNA', 'HP1a'),
            ('H1', 'DNA', 'LKH'),
            ('H1', 'DNA', 'Ki67'),
            ('mCD8T', 'WT-NCP', 'sp'),
            ('mCD8T', 'inht-NCP', 'sp'),
            ('mCD8T', 'KO-NCP', 'sp')]

exp_list = [('GM', 'NCP', 'sp')]

exp_list = [('H1', 'DNA', 'sp'),
            ('H1', 'DNA', 'HP1a'),
            ('H1', 'DNA', 'LKH'),
            ('H1', 'DNA', 'Ki67')]


exp_list = [('mCD8T', 'WT-NCP', 'sp'),
            ('mCD8T', 'inht-NCP', 'sp'),
            ('mCD8T', 'KO-NCP', 'sp')]


exp_list = [('H1', 'NCP', 'HP1a'),
            ('H1', 'NCP', 'LKH'),
            ('H1', 'NCP', 'Ki67'),
            ('H1', 'DNA', 'HP1a'),
            ('H1', 'DNA', 'LKH'),
            ('H1', 'DNA', 'Ki67')]

exp_list = [('mCD8T', 'KO-NCP', 'sp')]



# bin information
#bin_size = 1000
bin_size = 10000
#bin_size = 5000
#bin_size = 25000
score_type = 'score'
#score_type = 'Chalf'

for cell, sample, agent in exp_list:

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
    score_fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', score_type]) + '.cn'
    print "working on %s" % (score_fname)
    print "reading bin_score file"
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
    zscore_fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', 'z'+score_type]) + '.cn'
    f = open(zscore_fname, 'w')
    s = ['BinID', 'Chromosome', 'Start', 'End']
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
    print 
