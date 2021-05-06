import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
from sklearn import linear_model
import random
import matplotlib as mpl


def read_rlen(fname, chr_list=None):
    rlen_count = {}
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        line = line.strip()
        if not line:
            continue
        cols = line.strip().split()
        chr, rlen, count = cols
        if chr_list and chr not in chr_list:
            continue
        rlen, count = int(rlen), int(count)
        if rlen not in rlen_count:
            rlen_count[rlen] = 0
        rlen_count[rlen] += count
    return rlen_count

#path = "/home/spark159/../../media/spark159/sw/sp_spd_tests_detail/"
path = ""
#fnames = ["H1_DNA_spd_H1-DNA-spd-0_rlen.txt", "H1_DNA_sp_H1-DNA-sp-0_rlen.txt", "H1_NCP-new_spd_H1-NCP-new-spd-0_rlen.txt", "H1_NCP-new_sp_H1-NCP-new-sp-0_rlen.txt"]
#names = ["H1-DNA-spd-0", "H1-DNA-sp-0", "H1-NCP-spd-0", "H1-NCP-sp-0"]

#fnames = ["H1_NCP-new_sp_H1-NCP-new-sp-%s_rlen.txt" % (i) for i in range(6)]
#names = ["H1-NCP-sp-%s" % (i) for i in range(6)]

#fnames = ["H1_NCP-new_spd_H1-NCP-new-spd-%s_rlen.txt" % (i) for i in range(6)]
#names = ["H1-NCP-spd-%s" % (i) for i in range(6)] 

#fnames = ["H1_DNA_sp_H1-DNA-sp-%s_rlen.txt" % (i) for i in range(6)]
#names = ["H1-DNA-sp-%s" % (i) for i in range(6)] 

#fnames = ["H1_DNA_spd_H1-DNA-spd-%s_rlen.txt" % (i) for i in range(6)]
#names = ["H1-DNA-sp-%s" % (i) for i in range(6)] 


rlen_count_list = []
X_list, Y_list = [], []
for fname in fnames:
    rlen_count = read_rlen(path+fname)
    X = sorted(rlen_count.keys())
    Y = [ rlen_count[x] for x in X ]
    rlen_count_list.append(rlen_count)
    X_list.append(X)
    Y_list.append(Y)

color_list = np.linspace(0.01, 0.7, num=len(names))
alpha_list = np.linspace(0.3, 1, num=len(names))[::-1]
cmap = mpl.cm.get_cmap("copper")

fig = plt.figure()
for i in range(len(names)):
    name = names[i]
    X, Y = X_list[i], Y_list[i]
    plt.plot(X, Y, color=cmap(color_list[i]), label=name, alpha=1)
    #plt.plot(X, Y, color='b', label=name, alpha=alpha_list[i])
plt.xticks(range(0,500,50))
plt.grid(True)
plt.xlabel("Read length (bp)")
plt.ylabel("Read counts")
plt.title("Read length distribution")
plt.legend()
#plt.savefig(path+"rlen.png")
plt.savefig(fname.split('.')[0]+".png")
#plt.show()
plt.close()

"""
rlen_count1 = read_rlen(path+"hg19_sp1_rlen.txt")
rlen_count2 = read_rlen(path+"hg19_sp7_rlen.txt")
rlen_count3 = read_rlen(path+"hg19_sp8_rlen.txt")
    
X1 = sorted(rlen_count1.keys())
Y1 = [ rlen_count1[x] for x in X1 ]
X2 = sorted(rlen_count2.keys())
Y2 = [ rlen_count2[x] for x in X2 ]
X3 = sorted(rlen_count3.keys())
Y3 = [ rlen_count3[x] for x in X3 ]

fig = plt.figure()
plt.plot(X1, Y1, label='Input')
plt.plot(X2, Y2, label='NCP Sp6')
plt.plot(X3, Y3, label='NCP Sp7')
plt.xticks(range(0,500,50))
plt.grid(True)
plt.xlabel("Read length (bp)")
plt.ylabel("Read counts")
plt.title("Read length distribution")
plt.legend()
#plt.show()
plt.savefig(path+"rlen.png")
plt.close()
"""
