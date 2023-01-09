import math
import subprocess
import random

def write_temp_para (fname, energy_list, valency_list):
    f = open(fname + '_para.cn', 'w')
    print >> f, "%s\t%s\t%s\t%s\t%s" % ("SNP", "Chromosome", "PhysicalPosition", "Score", "Valency")
    for i in range(len(energy_list)):
        print >> f, "%s\t%s\t%s\t%s\t%s" % (str(i), "chr1", str(10*i), energy_list[i], valency_list[i])
    f.close()

def gel_point (N, f, volume):
    return math.log(N) + math.log(f) + 2*math.log(f-2) - math.log(volume) - math.log(f-1)

### parameters
volume = 10**7 

### test simulation
N = 100
f = 10
E = gel_point (N, f, volume)

energy_list = []
valency_list = []
for i in range(N):
    #energy = E + random.gauss(0,5)
    energy = E
    energy_list.append(energy)
    #valency = f + int(random.gauss(0,5))
    valency = f
    valency_list.append(valency)
fname = 'hetero_%s_%s_%s' % (N, f, E)
write_temp_para (fname, energy_list, valency_list)
subprocess.call(['python',
                 'rgs_submit.py',
                 fname + '_para.cn',
                 '--bethe',
                 '-o', fname,
                 '--replica', str(10),
                 '--max-root', str(100000),
                 '--cycle', str(100000)])

