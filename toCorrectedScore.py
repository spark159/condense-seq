import math
import numpy as np
import sys

fname = "H1_NCP_sp_chr1_Ncov.cn"

# get total of all nucleosome coverage
First = True
for line in open(fname):
    if not line.strip():
        continue
    cols = line.strip().split()
    if First:
        names = cols[3:-1]
        total_covs = [0.0] * (len(names) + 1)
        First = False
        continue
    counts = cols[3:]
    for i in range(len(counts)):
        count = 1 + float(counts[i]) # add 1 to all data
        #count = float(counts[i])
        #if count <= 0:
        #    count = 1.0
        total_covs[i] += count
print total_covs

total_num = 1.6*(10**12)
#total_fracs = [0.274827749, 0.152922919, 1.0] # survival fraction for titration #4, #8, #0
total_fracs = [0.150847925,0.008715373, 1.0] # bg corrected survival fraction for titration #4, #8, #0

f = open("H1_NCP_sp_chr1_num.cn", 'w')
s = 'SNP\tChromosome\tPhysicalPosition'
First = True
for line in open(fname):
    if not line.strip():
        continue
    cols = line.strip().split()
    if First:
        names = cols[3:]
        for name in names:
            s += '\t' + name.rsplit('/')[-1]
        s += '\t' + 'note'
        print >> f, s
        First = False
        continue
    ID, chr, pos = cols[:3]
    s = ID + '\t' + chr + '\t' + pos
    num_list = []
    counts = cols[3:]
    for i in range(len(counts)):
        count = 1 + float(counts[i])
        #count = float(counts[i])
        frac = float(count)/total_covs[i]
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
f.close()

# Sanity check
fname = "H1_NCP_sp_chr1_num.cn"
total, weird = 0, 0
First = True
for line in open(fname):
    if not line.strip():
        continue
    cols = line.strip().split()
    if First:
        names = cols[3:-1]
        total_nums = [0.0] * (len(names))
        First = False
        continue
    if cols[-1] == '*':
        counts = cols[3:-1]
    else:
        counts = cols[3:]
    nums = [int(value) for value in counts]
    if nums[-1] < nums[0] or nums[-1] < nums[1]:
        weird +=1
        #print line
    total+=1
    #assert nums[-1] >= nums[0]
    #assert nums[-1] >= nums[1]
    for i in range(len(nums)):
        total_nums[i] += nums[i]
print total_nums[-1] / total_num
print total_nums[0]/total_nums[-1]
print total_nums[1]/total_nums[-1]
print str(100*float(weird)/total) + ' %' 

sys.exit(1)

total_fracs = [0.274827749, 0.152922919] # survival fraction for titration #4, #8

# covert coverage to score
def metric (test, control):
    rcount = float(test)/float(control)
    return -math.log(rcount)

f = open("H1_NCP_sp_chr1_crscore.cn", 'w')
s = 'SNP\tChromosome\tPhysicalPosition'
First = True
for line in open(fname):
    if not line.strip():
        continue
    cols = line.strip().split()
    if First:
        names = cols[3:-1]
        for name in names:
            s += '\t' + name
        print >> f, s
        First = False
        continue
    ID, chr, pos = cols[:3]
    s = ID + '\t' + chr + '\t' + pos
    control = 1 + float(cols[-1])
    control = control/total_covs[-1]
    counts = cols[3:-1]
    for i in range(len(counts)):
        count = 1 + float(counts[i])
        count = count/total_covs[i]
        #score = round(metric(count, control), 5)
        score = -math.log(total_fracs[i]) + round(metric(count, control), 5)
        s += '\t' + str(score)
    print >> f, s
f.close()
