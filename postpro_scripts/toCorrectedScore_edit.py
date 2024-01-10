import math
import numpy as np
import sys
import load_file

# read titration file
def read_titration (fname, bg=False):
    all_fracs = {}
    tnum_tfrac = {}
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        cols = line.strip().split()
        try:
            tnum = int(cols[-1])
        except:
            continue
        total_frac = float(cols[-3])
        assert tnum not in tnum_tfrac
        tnum_tfrac[tnum] = total_frac
    return tnum_tfrac

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

#exp_list = [('H1', 'NCP', 'HP1a'),
#            ('H1', 'NCP', 'HP1bSUV'),
#            ('H1', 'NCP', 'LKH'),
#            ('H1', 'NCP', 'Ki67'),
#            ('H1', 'NCP', 'FUS')]

#exp_list = [('H1', 'NCP', 'HP1a'),
#            ('H1', 'NCP', 'Ki67'),
#            ('H1', 'NCP', 'LKH')]

#exp_list = [('H1', 'NCP', 'Ki67'),
#            ('H1', 'NCP', 'LKH')]

#exp_list = [('H1', 'NCP', 'sp')]

exp_list = [('H1', 'NCP', 'HP1a')]



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

    # total nucleosome input number
    total_num = 1.6*(10**12)  

    # read titration data
    tfname = '_'.join([cell, sample, agent, 'titration']) + '.csv'
    tnum_tfrac = read_titration(tfname)


    # get total of all nucleosome coverage
    names = []
    total_covs = []
    total_fracs = []
    print "get total coverage"
    for chr in chr_list:
        print "reading %s" % (chr)
        fname = path + '_'.join([cell, sample, agent, chr]) + '_Ncov.cn' 
        First = True
        for line in open(fname):
            if not line.strip():
                continue
            cols = line.strip().split()
            if First:
                if not names:
                    names = [name.rsplit('/')[-1].replace('-new','') for name in cols[3:]]
                else:
                    assert names == [name.rsplit('/')[-1].replace('-new','') for name in cols[3:]]

                if not total_covs:
                    total_covs = [0.0 for i in range(len(names))]

                if not total_fracs:
                    for name in names:
                        tnum = int(name.rsplit('.', 1)[0].split('-')[-1])
                        total_fracs.append(tnum_tfrac[tnum])
                    assert total_fracs[-1] == 1 # last data is control
                    #total_covs = [0.0] * (len(names) + 1)
                First = False
                continue
            covs = cols[3:]
            for i in range(len(covs)):
                #count = 1 + float(counts[i]) # add 1 to all data
                cov = float(covs[i])
                #count = float(counts[i])
                #if count <= 0:
                #    count = 1.0
                total_covs[i] += cov
    print total_covs
    print 

    # convert to real number
    print "calculate real number"
    for chr in chr_list:
        print "reading %s" % (chr)
        input_fname = path + '_'.join([cell, sample, agent, chr]) + '_Ncov.cn'
        output_fname = path + '_'.join([cell, sample, agent, chr]) + '_num.cn'
        f = open(output_fname, 'w')
        s = 'SNP\tChromosome\tPhysicalPosition'
        First = True
        for line in open(input_fname):
            if not line.strip():
                continue
            cols = line.strip().split()
            if First:
                names == [name.rsplit('/')[-1] for name in cols[3:]]
                for name in names:
                    s += '\t' + name
                s += '\t' + 'note'
                print >> f, s
                First = False
                continue
            ID, chr, pos = cols[:3]
            s = '\t'.join([ID, chr, pos])
            num_list = []
            covs = cols[3:]
            for i in range(len(covs)):
                #count = 1 + float(counts[i])
                cov = float(covs[i])
                frac = float(cov)/total_covs[i]
                num = total_num * total_fracs[i] * frac
                #lognum = np.log(total_num) + np.log(total_fracs[i]) + np.log(count) - np.log(total_covs[i])
                #num = np.exp(lognum)
                num = int(round(num))
                num_list.append(num)
                s += '\t' + str(num)
            check = True
            for i in range(len(num_list)-1):
                if num_list[i] > num_list[-1]:
                    check = False
                    break
            if not check:
                s += '\t' + '*'
            print >> f, s
        f.close()
    print

    # Sanity check
    total_nums = [0.0 for i in range(len(names))]
    print "Sanity check"
    total, weird = 0, 0
    star = 0
    for chr in chr_list:
        First = True
        fname = path + '_'.join([cell, sample, agent, chr]) + '_num.cn'
        for line in open(fname):
            if not line.strip():
                continue
            cols = line.strip().split()
            if First:
                names = cols[3:-1]
                #total_nums = [0.0] * (len(names))
                First = False
                continue
            if cols[-1] == '*':
                counts = cols[3:-1]
                star +=1
            else:
                counts = cols[3:]
            nums = [int(value) for value in counts]
            for num in nums[:-1]:
                if num > nums[-1]:
                    assert cols[-1] == '*'
                    weird +=1
                    break
            total+=1
            #assert nums[-1] >= nums[0]
            #assert nums[-1] >= nums[1]
            for i in range(len(nums)):
                total_nums[i] += nums[i]


    print total_fracs[-1], total_nums[-1] / total_num
    for i in range(len(total_fracs)-1):
        print total_fracs[i], total_nums[i]/total_nums[-1]
    print str(100*float(star)/total) + '%', str(100*float(weird)/total) + '%' 
    print
