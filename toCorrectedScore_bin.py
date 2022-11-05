import math
import numpy as np
import sys


# parameters
path = "/home/spark159/../../media/spark159/sw/"

cell = 'H1'
#cell = 'GM'
#cell = 'mCD8T'
sample = 'NCP'
#sample = 'WT-NCP'
#sample = 'inht-NCP'
agent = 'sp'

bin_size = 1000
#offset = 1 # only GC content
offset = 2 # GC content and template length


total_num = 1.6*(10**12) # total nucleosome input number

# H1 NCP sp4+
total_fracs = [0.178899286, 0.028922932, 1.0] # bg corrected survival fraction for titration #4, #8, #0

# GM NCP sp4+
#total_fracs = [0.113430444, 0.005803754, 1.0] # bg corrected survival fraction for titration #4, #8, #0

# H1 NCP HP1a
#total_fracs = [0.108659824, 1.0] # bg corrected survival fraction for titration #3, #0

# mouse CD8 Tcell (WT)
#total_fracs = [0.0720702, 0.018194008, 1.0] # bg corrected survival fraction for titration #4, #8, #0

# mouse CD8 Tcell (+inht)
#total_fracs = [0.084142467, 0.013857276, 1.0] # bg corrected survival fraction for titration #4, #8, #0


# get total of all nucleosome coverage
total_counts = [0.0 for i in range(len(total_fracs))]
print "get total counts"
fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb']) + '_bin.cn' 
First = True
for line in open(fname):
    if not line.strip():
        continue
    cols = line.strip().split()
    if First:
        names = cols[4:-offset]
        #total_covs = [0.0] * (len(names) + 1)
        First = False
        continue
    ID, chr, start, end = cols[:4]
    chrnum = chr[3:]
    try:
        int(chrnum)
    except:
        if chrnum not in ['X', 'Y']:
            continue
    counts = cols[4:-offset]
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
        names = cols[4:-offset]
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
        if chrnum not in ['X', 'Y']:
            continue
    s = '\t'.join([ID, chr, start, end])
    num_list = []
    counts = cols[4:-offset]
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
    else:
        counts = cols[4:]
    nums = [int(value) for value in counts]
    if nums[-1] < nums[0] or nums[-1] < nums[1]:
        weird +=1
        #print line
    total+=1
    #assert nums[-1] >= nums[0]
    #assert nums[-1] >= nums[1]
    for i in range(len(nums)):
        total_nums[i] += nums[i]

print total_fracs[-1], total_nums[-1] / total_num
print total_fracs[0], total_nums[0]/total_nums[-1]
print total_fracs[1], total_nums[1]/total_nums[-1]
print str(100*float(weird)/total) + ' %' 
print
