import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import math
import pickle

def AT_content(seq):
    seq = seq.upper()
    output = 0.0
    for nt in seq:
        if nt in "AT":
            output += 1.0
    return output/len(seq)


def read_ref (ref_fname):
    id_seq = {}
    for line in open(ref_fname):
        line = line.strip()
        if line.startswith('>'):
            id = int(line[4:])
            continue
        if line:
            assert id not in id_seq
            id_seq[id] = line
    return id_seq

def read_ncount(fname):
    id_ncount = {}
    for line in open(fname):
        if line.strip():
            id, seq, ncount = line.strip().split()
            id = int(id)
            assert id not in id_ncount
            id_ncount[id] = float(ncount)
    return id_ncount

def read_DNAshape(fname):
    names = ['MGW', 'HelT', 'ProT', 'Roll']
    dic_list = [{} for i in range(len(names))]
    for i in range(len(names)):
        data = []
        for line in open(fname+"."+names[i]):
            line = line.strip()
            if line.startswith('>'):
                if data:
                    assert id not in dic_list[i]
                    dic_list[i][id] = data
                id = int(line[1:].strip().split('_')[1])
                data =[]
                continue
            if line:
                temp = line.split(',')
                for k in range(len(temp)):
                    try:
                        temp[k] = float(temp[k])
                    except Exception:
                        pass
                data += temp
        assert id not in dic_list[i]
        dic_list[i][id] = data
    return dic_list

def read_titration (fname):
    conc_list = []
    mean_list = []
    std_list = []
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        if not line:
            continue
        cols = line.strip().split()
        conc, mean, std = float(cols[0]), float(cols[1]), float(cols[2])
        conc_list.append(conc)
        mean_list.append(mean)
        std_list.append(std)
    return conc_list, mean_list, std_list

def GC_content(seq):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    return (num/float(len(seq)))*100

def average(l):
    sum = 0.0
    count =0
    for e in l:
        if type(e) != str:
            count +=1
            sum += e
    return sum/count

id_seq = read_ref('Blib.ref')
id_MGW, id_HelT, id_ProT, id_Roll = read_DNAshape('php6POcc7')
_, mean_list, _ = read_titration("AakashDNAlib_spermine.csv")
fidx_to_tnum = [None, 4, 7, 8, [9,10], 11, 12]
id_ncount3 = read_ncount('ncount_sp4.txt')
id_ncount5 = read_ncount('ncount_sp6.txt')

out_fname = 'YWlib'
f = open(out_fname + '_score.txt', 'w')
s = 'ID\tSequence\tATcontent\tMGW\tHelT\tProT\tRoll\t'
s += 'Sp3\tSp5'
print >> f, s

for id in sorted(id_seq.keys()):
    seq = id_seq[id]
    ATcontent = AT_content(seq)
    MGW, HelT, ProT, Roll = id_MGW[id], id_HelT[id], id_ProT[id], id_Roll[id]

    data = [seq, str(AT_content), str(MGW), str(HelT), str(ProT), str(Roll)]

    mean_fract3 = mean_list[fidx_to_tnum[3]]
    mean_fract5 = mean_list[fidx_to_tnum[5]]
    try:
        score3 = -np.log2(id_ncount3[id]) -np.log2(mean_fract3)
    except:
        print id

    score5 = -np.log2(id_ncount5[id]) -np.log2(mean_fract5)

    data.append(str(score3))
    data.append(str(score5))
        
    print >> f, '\t'.join(data)

f.close()
