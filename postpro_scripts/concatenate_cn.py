import math

# parameters
#path = "/home/spark159/../../media/spark159/sw/"
path = "/home/spark159/../../storage/"


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

#exp_list = [('H1', 'NCP', 'HP1a')]

exp_list = [('H1', 'NCP', 'sp'),
            ('H1', 'NCP', 'HP1a'),
            ('H1', 'NCP', 'Ki67'),
            ('H1', 'NCP', 'LKH')]

#exp_list = [('H1', 'NCP', 'HP1a')]

dtype = 'zscore'

for cell, sample, agent in exp_list:

    print "Processing %s %s %s" % (cell, sample, agent)
    
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


    new_fname = path + '_'.join([cell, sample, agent, dtype]) + '.cn'
    f = open(new_fname, 'w')

    ID = 0
    for i in range(len(chr_list)):
        print "reading %s" % (chr_list[i])
        score_fname = path + '_'.join([cell, sample, agent, chr_list[i], dtype]) + '.cn'
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
    print
