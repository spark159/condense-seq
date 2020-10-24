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

"""
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

#check_acf(x2)

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

    
path = "./data/"
jump = 10
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path + "hg19_chr1_167win25step_anot.cn", jump=jump)
ID_score1 = name_ID_value['data/sp_spd_tests_detail/sp7']

#names = ['data/sp_spd_tests_detail/sp7', 'ATcontent', 'CpGNumber', 'random']
names = ['data/sp_spd_tests_detail/sp7', 'ATcontent', 'CpGNumber']
step_size = 25*jump

# select the longest data list with constant step size and wihtout np.nan
posID = [[pos, ID] for ID, pos in ID_pos.items()]
posID = sorted(posID, cmp=tuple_cmp)
IDs = [ value[1] for value in posID ]

IDs_list = []
i = 0
while i < len(IDs):
    ID = IDs[i]
    score = ID_score1[ID]
    if not np.isnan(score):
        temp = []
        j = 0
        while i + j < len(IDs):
            ID = IDs[i+j]
            score = ID_score1[ID]
            if np.isnan(score):
                break
            if j > 0:
                prev_ID = temp[-1]
                pos_change = ID_pos[ID] - ID_pos[prev_ID] 
                if pos_change != step_size:
                    break
            temp.append(ID)
            j += 1
        IDs_list.append(temp)
        i += j
    else:
        i +=1

max_num = 0
for IDs in IDs_list:
    if len(IDs) > max_num:
        max_num = len(IDs)
        selected_IDs = IDs

# get acfs
name_acf = {}
for name in names:
    if name == 'random':
        value_list = [ID_score1[ID] for ID in selected_IDs]
        random.shuffle(value_list)
        acf = statis.acf(value_list)
    else:
        value_list = [name_ID_value[name][ID] for ID in selected_IDs]
        acf = statis.acf(value_list)
    name_acf[name] = acf

# plot acfs
def symlog (value):
    if value >= 0:
        sign = 1.0
    else:
        sign = -1.0
    return sign*np.log(1+value)

fig = plt.figure()
for name in names:
    acf = name_acf[name]
    X = np.asarray([i*step_size for i in range(len(acf))])
    if name == 'data/sp_spd_tests_detail/sp7':
        name = "Condensability"
        alpha = 1
    elif name == 'ATcontent':
        name = 'AT content'
        alpha = 0.5
    elif name == 'CpGNumber':
        name = 'CpG number'
        alpha = 0.5
    #plt.plot(X+1, acf, label=name, alpha=alpha)
    #acf = np.asarray([symlog(value) for value in acf])
    #acf = np.exp(acf)
    plt.plot(X[1:], acf[1:], label=name, alpha=alpha)
plt.xlabel("Distance (bp)")
plt.ylabel("Correlation")
plt.xscale("log")
#plt.yscale("symlog", basey=10)
leg = plt.legend()
for lh in leg.legendHandles:
    lh._legmarker.set_alpha(1)
plt.show()
plt.close()
