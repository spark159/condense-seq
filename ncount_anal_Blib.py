import matplotlib.pyplot as plt
import numpy as np
import math

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
    seq_ncount = {}
    for line in open(fname):
        if line.strip():
            id, seq, ncount = line.strip().split()
            assert id not in seq_ncount
            seq_ncount[seq] = float(ncount)
    return seq_ncount

def read_DNAshape(fname, id_seq):
    names = ['MGW', 'HelT', 'ProT', 'Roll']
    dic_list = [{} for i in range(len(names))]
    for i in range(len(names)):
        data = []
        for line in open(fname+"."+names[i]):
            line = line.strip()
            if line.startswith('>'):
                if data:
                    assert seq not in dic_list[i]
                    dic_list[i][seq] = data
                id = int(line[1:].strip().split('_')[1])
                seq = id_seq[id]
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
        assert seq not in dic_list[i]
        dic_list[i][seq] = data
    return dic_list

def average(l):
    sum = 0.0
    count =0
    for e in l:
        if type(e) != str:
            count +=1
            sum += e
    return sum/count

def get_corr(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = np.average(x)
    avg_y = np.average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff
    return diffprod / math.sqrt(xdiff2 * ydiff2)

id_seq = read_ref('Blib.ref')
seq_MGW, seq_HelT, seq_ProT, seq_Roll = read_DNAshape('php6POcc7', id_seq)

fname_list = ['ncount_sp' + str(i+1) + '.txt' for i in range(7)]

for i in range(len(fname_list)):
    fname = fname_list[i]
    seq_ncount = read_ncount(fname)
    X1, X2, X3, X4, Y = [], [], [], [], []
    
    for seq in seq_ncount:
        X1.append(average(seq_MGW[seq]))
        X2.append(average(seq_HelT[seq]))
        X3.append(average(seq_ProT[seq]))
        X4.append(average(seq_Roll[seq]))
        Y.append(seq_ncount[seq])

    names = ['MGW', 'HelT', 'ProT', 'Roll']
    Xlist = [X1,X2,X3,X4]
    f = open("corr_sp" + str(i+1) + ".txt", 'w')
    for j in range(len(Xlist)):
        X = Xlist[j]
        fig = plt.figure()
        plt.scatter(X,Y)
        plt.xlabel(names[j])
        plt.ylabel('Norm. counts')
        plt.savefig(names[j] + 'VSncount_sp' +str(i+1) + ".png")
        plt.close()
        corr = get_corr(X,Y)
        print >> f, names[j], corr
    f.close()

"""
x, y, z = [], [], []
for GC in GC_ncount:
    x.append(GC)
    y.append(np.mean(GC_ncount[GC]))
    z.append(np.std(GC_ncount[GC]))
fig = plt.figure()
plt.plot(x,y,'.')
plt.errorbar(x,y,yerr=z,fmt='o')
plt.xlabel('GC contents')
plt.ylabel('Norm. counts')
#$plt.show()
plt.savefig('GCVScount.png'+ str(i+1) + ".png")
plt.close()
"""
