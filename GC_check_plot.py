import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import seaborn as sns
import copy
import math
from scipy import signal

def read_GC (fname):
    data_list = []
    Find = False
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        if line.startswith("@GCcontent"):
            Find = True
            continue
        if not Find:
            continue
        if Find and line.startswith('@'):
            break
        cols = line.strip().split(',')
        data_list += [float(value) for value in cols]
    return data_list


for k in range(4):
    title = ""
    if k % 2 == 0:
        title += "NCP_"
    else:
        title += "DNA_"
    if k < 2:
        title += "Spermine(4+)"
    else:
        title += "Spermidine(3+)"
    
    Xrange = range(8*k+1, 8*(k+1)+1)
    mean_list, std_list = [], []
    for i in Xrange:
        filename = "/home/spark159/Downloads/sp_spd_test" + str(i) + "_qc.txt"
        data_list = read_GC (filename)
        mean_list.append(np.mean(data_list))
        std_list.append(np.std(data_list)/np.sqrt(len(data_list)))

    fig = plt.figure()
    plt.plot(range(1,9), mean_list, '.')
    plt.errorbar(range(1,9), mean_list, yerr=std_list, fmt='.')
    plt.ylim([0.25, 0.6])
    plt.title(title)
    plt.xlabel("Titration point")
    plt.ylabel("mean GC content")
    plt.savefig(title+'.png')
    #plt.show()
    plt.close()


