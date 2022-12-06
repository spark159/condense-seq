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
path = "/home/spark159/../../media/spark159/sw/"

#cell = 'H1'
#cell = 'GM'
cell = 'mCD8T'
#sample = 'NCP'
#sample = 'WT-NCP'
#sample = 'inht-NCP'
sample = 'KO-NCP'
agent = 'sp'

#bin_size = 10000
bin_size = 1000
skip = 0 # skip first titration points
#offset = 1 # only GC content
offset = 2 # GC content and template length
XY_chrom = ['X']
#XY_chrom = ['X', 'Y']

total_num = 1.6*(10**12) # total nucleosome input number

# read titration data
tfname = '_'.join([cell, sample, agent, 'titration']) + '.csv'
tnum_tfrac = read_titration(tfname)

## Raw nanodrop data
# H1 NCP sp4+
#total_fracs = [0.243182021, 0.104906831, 1.0] # bg corrected survival fraction for titration #4, #8, #0
# GM NCP sp4+
#total_fracs = [0.233804495, 0.140807819, 1.0] # bg corrected survival fraction for titration #4, #8, #0
# H1 NCP HP1a
#total_fracs = [0.108659824, 1.0] # bg corrected survival fraction for titration #3, #0
# mouse CD8 Tcell (WT)
#total_fracs = [0.169016035, 0.120853265, 1.0] # bg corrected survival fraction for titration #4, #8, #0
# mouse CD8 Tcell (+inht)
#total_fracs = [0.189020626, 0.12685576, 1.0] # bg corrected survival fraction for titration #4, #8, #0
# mouse CD8 Tcell (KO)
#total_fracs = [0.185366333, 0.117691975, 1.0] # bg corrected survival fraction for titration #4, #8, #0


## backgroud corrected nanodrop data
# H1 NCP sp4+
#total_fracs = [0.178899286, 0.028922932, 1.0] # bg corrected survival fraction for titration #4, #8, #0
# GM NCP sp4+
#total_fracs = [0.113430444, 0.005803754, 1.0] # bg corrected survival fraction for titration #4, #8, #0
# H1 NCP HP1a
#total_fracs = [0.108659824, 1.0] # bg corrected survival fraction for titration #3, #0
# mouse CD8 Tcell (WT)
#total_fracs = [0.0720702, 0.018194008, 1.0] # bg corrected survival fraction for titration #4, #8, #0
# mouse CD8 Tcell (+inht)
#total_fracs = [0.084142467, 0.013857276, 1.0] # bg corrected survival fraction for titration #4, #8, #0
# mouse CD8 Tcell (KO)
#total_fracs = [0.081973904, 0.005736345, 1.0] # bg corrected survival fraction for titration #4, #8, #0


# get total of all nucleosome coverage
#total_counts = [0.0 for i in range(len(total_fracs))]
#total_counts = []
print "get total counts"
fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb']) + '_bin.cn' 
First = True
for line in open(fname):
    if not line.strip():
        continue
    cols = line.strip().split()
    if First:
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
    chrnum = chr[3:]
    try:
        int(chrnum)
    except:
        if chrnum not in XY_chrom:
            continue
    counts = cols[4+skip:-offset]
    for i in range(len(counts)):
        #count = 1 + float(counts[i]) # add 1 to all data
        count = float(counts[i])
        #count = float(counts[i])
        #if count <= 0:
        #    count = 1.0
        total_counts[i] += count
print total_counts
print 

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
