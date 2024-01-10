import math

# parameters
#path = "/home/spark159/../../media/spark159/sw/"
#path = "/home/spark159/../../storage/"
path = "/home/spark159/../../storage/replicates/"

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

#exp_list = [('H1', 'NCP', 'LKH')]\

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
#bin_size = 5000
#bin_size = 1000
bin_size = 10000
#bin_size = 25000

# skip choice
skip_zero = True
skip_star = False

for cell, sample, agent in exp_list:

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

    num_fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb']) + '_num.cn'
    score_fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb']) + '_score.cn'

    print "working on %s" % (num_fname)
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
            if skip_star:
                continue
            cols = cols[:-1]
        ID, chr, st, ed = cols[:4]
        st, ed = int(st), int(ed)
        counts = [int(count) for count in cols[4:]]
        if skip_zero and 0 in counts:
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
            s = [str(ID), chr, str(st), str(ed)]
            s += [str(score) for score in chr_range_scores[chr][(st,ed)]]
            print >> f, '\t'.join(s)
            ID +=1
    f.close()
    print
