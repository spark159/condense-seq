import math

fname = "hg19_chr1_cov.cn"

"""
# get total counts
First = True
for line in open(fname):
    if not line.strip():
        continue
    cols = line.strip().split()
    if First:
        names = cols[3:-1]
        total_list = [0.0] * (len(names) + 1)
        First = False
        continue
    counts = cols[3:]
    for i in range(len(counts)):
        count = float(counts[i])
        if count <= 0:
            count = 1.0
        total_list[i] += count
print total_list
"""

total_list = [4027035549.0, 4106574224.0, 3867424201.0]

# covert coverage to score
def metric (test, control):
    if test <=0:
        test = 1.0
    if control <= 0:
        control = 1.0
    rcount = float(test)/float(control)
    return -math.log(rcount)
f = open("score.cn", 'w')
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
    control = float(cols[-1])
    if control <= 0:
        control = 1.0
    control = control/total_list[-1]
    counts = cols[3:-1]
    for i in range(len(counts)):
        count = float(counts[i])
        if count <= 0:
            count = 1.0
        count = count/total_list[i]
        score = round(metric(count, control), 5)
        s += '\t' + str(score)
    print >> f, s
f.close()


