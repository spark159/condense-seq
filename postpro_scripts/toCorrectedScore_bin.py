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

#exp_list = [('H1', 'NCP', 'HP1a'),
#            ('H1', 'NCP', 'HP1bSUV'),
#            ('H1', 'NCP', 'LKH'),
#            ('H1', 'NCP', 'Ki67'),
#            ('H1', 'NCP', 'FUS')]

#exp_list = [('H1', 'NCP', 'LKH')]

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


exp_list = [('H1', 'NCP', 'HP1a'),
            ('H1', 'NCP', 'LKH'),
            ('H1', 'NCP', 'Ki67'),
            ('H1', 'DNA', 'HP1a'),
            ('H1', 'DNA', 'LKH'),
            ('H1', 'DNA', 'Ki67')]

#exp_list = [('mCD8T', 'WT-NCP', 'sp'),
#            ('mCD8T', 'inht-NCP', 'sp'),
#            ('mCD8T', 'KO-NCP', 'sp')]


exp_list = [('mCD8T', 'KO-NCP', 'sp')]


bin_size = 10000
#bin_size = 1000
#bin_size = 5000
#bin_size = 25000
skip = 0 # skip first titration points
#offset = 1 # only GC content
#offset = 2 # GC content and template length

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

    # set total nucleosome input number
    total_num = 1.6*(10**12) 

    if agent in ['Ca', 'Mg']:
        total_num = 5*total_num

    # read titration data
    tfname = '_'.join([cell, sample, agent, 'titration']) + '.csv'
    tnum_tfrac = read_titration(tfname)

    # get sum of all nucleosome coverage
    fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb']) + '_bin.cn'
    print "working on %s" % (fname)
    print "get total counts"

    chrnum_range_counts = {}
    First = True
    for line in open(fname):
        if not line.strip():
            continue
        cols = line.strip().split()
        if First:
            if cols[-1] == 'Meantlen':
                offset = 2
            else:
                offset = 1
            names = cols[4+skip:-offset]
            total_counts = [0.0 for i in range(len(names))]
            total_fracs = []
            for name in names:
                tnum = int(name.rsplit('.', 1)[0].split('-')[-1])
                total_frac = tnum_tfrac[tnum]
                total_fracs.append(total_frac)
            assert total_fracs[-1] == 1
            #total_covs = [0.0] * (len(names) + 1)
            First = False
            continue

        ID, chr, start, end = cols[:4]

        if chr not in chr_list:
            continue

        chrnum = chr[3:]
        try:
            chrnum = int(chrnum)
        except:
            pass

        start, end = int(start), int(end)

        counts = [float(count) for count in cols[4+skip:-offset]]

        for i in range(len(counts)):
            count = counts[i]
            total_counts[i] += count

        if chrnum not in chrnum_range_counts:
            chrnum_range_counts[chrnum] = {}
        chrnum_range_counts[chrnum][(start, end)] = counts

    print total_counts
    print


    # convert to real count
    print "calculate real number"

    fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb']) + '_num.cn'
    f = open(fname, 'w')
    s = '\t'.join(['BinID', 'Chromosome', 'Start', 'End'] + names + ['note'])
    print >> f, s

    ID = 0
    for chrnum in sorted(chrnum_range_counts):
        for start, end in sorted(chrnum_range_counts[chrnum]):
            chr = 'chr' + str(chrnum)

            counts = chrnum_range_counts[chrnum][(start, end)]

            s = '\t'.join([str(ID), chr, str(start), str(end)])
            num_list = []
            for i in range(len(counts)):
                count = counts[i]
                frac = float(count)/total_counts[i]
                num = total_num * total_fracs[i] * frac
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
            ID +=1

    f.close()

    # Sanity check
    total_nums = [0.0 for i in range(len(total_fracs))]
    print "Sanity check"
    total, weird = 0, 0
    star = 0
    First = True
    fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb']) + '_num.cn'
    for line in open(fname):
        if not line.strip():
            continue
        cols = line.strip().split()
        if First:
            names = cols[4:-1]
            #total_nums = [0.0] * (len(names))
            First = False
            continue
        if cols[-1] == '*':
            counts = cols[4:-1]
            star +=1
        else:
            counts = cols[4:]
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



"""
# convert to real count
print "calculate real number"
input_fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb']) + '_bin.cn'
output_fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb']) + '_num.cn'
f = open(output_fname, 'w')
s = 'BinID\tChromosome\tStart\tEnd'
First = True
for line in open(input_fname):
    if not line.strip():
        continue
    cols = line.strip().split()
    if First:
        names = cols[4+skip:-offset]
        for name in names:
            s += '\t' + name.rsplit('/')[-1]
        s += '\t' + 'note'
        print >> f, s
        First = False
        continue
    ID, chr, start, end = cols[:4]
    chrnum = chr[3:]
    try:
        int(chrnum)
    except:
        if chrnum not in XY_chrom:
            continue
    s = '\t'.join([ID, chr, start, end])
    num_list = []
    counts = cols[4+skip:-offset]
    for i in range(len(counts)):
        #count = 1 + float(counts[i])
        count = float(counts[i])
        frac = float(count)/total_counts[i]
        num = total_num * total_fracs[i] * frac
        #lognum = np.log(total_num) + np.log(total_fracs[i]) + np.log(count) - np.log(total_counts[i])
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
"""
