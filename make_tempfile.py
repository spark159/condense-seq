import subprocess
import random

def write_temp_para (fname, energy_list, valency_list):
    f = open(fname + '_para.cn', 'w')
    print >> f, "%s\t%s\t%s\t%s\t%s" % ("SNP", "Chromosome", "PhysicalPosition", "Score", "Valency")
    for i in range(len(energy_list)):
        print >> f, "%s\t%s\t%s\t%s\t%s" % (str(i), "chr1", str(10*i), energy_list[i], valency_list[i])
    f.close()


"""
N = 100
energy_list = []
valency_list = []
for i in range(N):
    energy_list.append(-100)
    valency_list.append(23)
write_temp_para ("test_N100_V23_E-100", energy_list, valency_list)
"""



"""
N = 100
energy_list = []
valency_list = []
for i in range(N):
    if i == 0:
        energy = -13
        valency = 56
    else:
        energy = -10
        valency = 23
    #energy = -10 + random.gauss(0,5)
    #energy = -10
    #valency = 23
    energy_list.append(energy)
    #valency = 23 + int(random.randint(-5,5))
    #valency = 23 + int(random.gauss(0,5))
    valency_list.append(valency)
write_temp_anot2 ("hetero_100_23_-10", energy_list, valency_list)
#write_temp_anot2 ("Just", energy_list, valency_list)
"""

"""
N = 100
f = 23
E = -10
for v in range(10):     
    f_single = 3 + 10*v
    for i in range(10):
        E_single = - 2*i
        energy_list = [E_single] + [E]*(N-1)
        valency_list = [f_single] + [f]*(N-1)
        fname = 'single_' + str(N) + '_' + str(f) + '_' + str(E) + ':' + str(f_single) + '_' + str(E_single)
        write_temp_anot2 (fname, energy_list, valency_list)
        subprocess.call(['python', 'rca_submit_edit.py', fname + '_anot.cn', '--model', 'Free', '--metric', 'raw', '--num', str(N), '--bethe', '-o', fname, '--cycle', '100000'])



N = 100
for v in range(10):
    valency = 3 + 10*v
    for i in range(10):
        energy = -2*i
        fname = "homo_" + str(N) + '_' + str(valency) + '_' + str(energy)
        write_temp_para (fname, [energy]*N, [valency]*N)
        subprocess.call(['python', 'rgs_submit.py', fname + '_para.cn', '--bethe', '-o', fname])
"""

"""
N = 100
valency = 23
for i in [4]:
    energy = -2*i
    print energy
    fname = "homo_" + str(N) + '_' + str(valency) + '_' + str(energy)
    write_temp_para (fname, [energy]*N, [valency]*N)
    subprocess.call(['python', 'rgs_submit.py', fname + '_para.cn', '--bethe', '-o', fname, '--cycle', str(20000), '--replica', str(10), '--max-root', str(100), '--sigma', str(1), '--burn', str(12000)])

"""

N = 100
f = 23
E = -10
for v in [3]:     
    f_single = 3 + 10*v
    for i in range(10)[::-1]:
        E_single = - 2*i
        energy_list = [E_single] + [E]*(N-1)
        valency_list = [f_single] + [f]*(N-1)
        fname = 'single_' + str(N) + '_' + str(f) + '_' + str(E) + ':' + str(f_single) + '_' + str(E_single)
        write_temp_para (fname, energy_list, valency_list)
        subprocess.call(['python', 'rgs_submit.py', fname + '_para.cn', '--bethe', '-o', fname, '--cycle', str(100000), '--replica', str(10), '--max-root', str(100), '--sigma', str(1)])
