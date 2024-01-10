import math

# parameters

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

#exp_list = [('H1', 'NCP', 'HP1a'),
#            ('H1', 'NCP', 'Ki67'),
#            ('H1', 'NCP', 'LKH')]

exp_list = [('H1', 'NCP', 'HP1a')]



# skip choice
skip_star = True

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
                if skip_star:
                    continue
                cols = cols[:-1]
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
    print
