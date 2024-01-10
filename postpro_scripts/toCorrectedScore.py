import math
import numpy as np
import sys

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

cell = 'H1'
#cell = 'GM'
#cell = 'mCD8T'
sample = 'NCP'
#sample = 'WT-NCP'
#sample = 'inht-NCP'
#sample = 'KO-NCP'
agent = 'sp'
#note = 'Ncov'
#note = '1001win501step_cov_Bsig'

chr_list = ['chr' + str(i) for i in range(1, 23)] #human
#chr_list = ['chr' + str(i) for i in range(1, 20)] #mouse
#chr_list += ['chrX']
chr_list += ['chrX', 'chrY']
#chr_list = ['chr' + str(i) for i in range(1, 2)]


total_num = 1.6*(10**12) # total nucleosome input number

## Raw nanodrop data
# H1 NCP sp4+
total_fracs = [0.243182021, 0.104906831, 1.0] # survival fraction for titration #4, #8, #0
# GM NCP sp4+
#total_fracs = [0.233804495, 0.140807819, 1.0] # survival fraction for titration #4, #8, #0
# H1 NCP HP1a
#total_fracs = [0.108659824, 1.0] # survival fraction for titration #3, #0
# mouse CD8 Tcell (WT)
#total_fracs = [0.169016035, 0.120853265, 1.0] # survival fraction for titration #4, #8, #0
# mouse CD8 Tcell (+inht)
#total_fracs = [0.189020626, 0.12685576, 1.0] # survival fraction for titration #4, #8, #0
# mouse CD8 Tcell (KO)
#total_fracs = [0.185366333, 0.117691975, 1.0] # survival fraction for titration #4, #8, #0


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
total_covs = [0.0 for i in range(len(total_fracs))]
print "get total coverage"
for chr in chr_list:
    print "reading %s" % (chr)
    #fname = path + '_'.join([cell, sample, agent, chr, note]) + '.cn'
    fname = path + '_'.join([cell, sample, agent, chr, 'Ncov']) + '.cn'
    #fname = "H1_NCP_sp_%s_Ncov.cn" % (chr)
    First = True
    for line in open(fname):
        if not line.strip():
            continue
        cols = line.strip().split()
        if First:
            names = cols[3:-1]
            #total_covs = [0.0] * (len(names) + 1)
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
print 

# convert to real count
print "calculate real number"
for chr in chr_list:
    print "reading %s" % (chr)
    #input_fname = "H1_NCP_sp_%s_Ncov.cn" % (chr)
    #output_fname = "H1_NCP_sp_%s_num.cn" % (chr)
    input_fname = path + '_'.join([cell, sample, agent, chr]) + '_Ncov.cn'
    output_fname = path + '_'.join([cell, sample, agent, chr]) + '_num.cn'
    #input_fname = path + '_'.join([cell, sample, agent, chr]) + '.cn'
    #output_fname = path + '_'.join([cell, sample, agent, chr, note]) + '_num.cn'
    #output_fname = path + '_'.join([cell, sample, agent, chr]) + '_num.cn'
    f = open(output_fname, 'w')
    s = 'SNP\tChromosome\tPhysicalPosition'
    First = True
    for line in open(input_fname):
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
total_nums = [0.0 for i in range(len(total_fracs))]
print "Sanity check"
for chr in chr_list:
    print "reading %s" % (chr)
    total, weird = 0, 0
    star = 0
    First = True
    #fname = "H1_NCP_sp_%s_num.cn" % (chr)
    fname = path + '_'.join([cell, sample, agent, chr]) + '_num.cn'
    #fname = path + '_'.join([cell, sample, agent, chr]) + '_num.cn'
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
