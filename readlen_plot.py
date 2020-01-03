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

path = "/home/spark159/../../media/spark159/sw/sp_spd_tests_detail/"

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

