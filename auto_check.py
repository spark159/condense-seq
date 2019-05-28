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


def check_acf (x):
    mean = np.mean(x)
    std = np.std(x)
    for lag in range(len(x)-10):
        X, Y = [], []
        temp = 0.0
        for i in range(len(x)-lag):
            X.append(x[i])
            Y.append(x[i+lag])
            temp += (x[i]-mean)*(x[i+lag]-mean)
        temp = float(temp)/(len(x)-lag)
        temp = temp/(std**2)
        fig = plt.figure()
        plt.plot(X, Y, '.')
        plt.plot(mean,mean, 'x')
        plt.title(str(lag) + " " + str(temp))
        plt.show()
        plt.close()
    return 

x1 = [ random.randint(0,1) for i in range(100)]
#x2 = [0]*100 + [0]*1 + [0]*100
#x2 = [0,0,1]*10
#x2 = [1]*10 + [0.0]*10
x2 = [i for i in range(100)]

check_acf(x2)

acf1 = statis.acf(x1)
acf2 = statis.acf(x2)

fig = plt.figure()
plt.plot(acf1, alpha=0.5, label='random')
plt.plot(acf2, alpha=0.5, label='test')
plt.legend()
plt.show()
plt.close()


"""
def tuple_cmp (a, b):
    if a[0] < b[0]:
        return -1
    else:
        return 1

ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("data/hg19_chr1_167win25step_new_anot.cn", num_max=100000)
ID_value = name_ID_value["work/condense_seq/sp10_hg19_chr1"]

posID = [[pos, ID] for ID, pos in ID_pos.items()]
posID = sorted(posID, cmp=tuple_cmp)
ID_list = [ value[1] for value in posID ]

value_list = [ID_value[ID] for ID in ID_list]
acf_list = statis.acf(value_list)
X = np.asarray([i*25 for i in range(len(acf_list))])
random_list = [random.random() for i in range(len(value_list))]
acf_list2 = statis.acf(random_list)

fig = plt.figure()
plt.plot(np.log10(X+1), acf_list, label="condensability", alpha=0.8)
plt.plot(np.log10(X+1), acf_list2, label="random", alpha=0.8)
plt.legend()
plt.show()
plt.close()
"""
