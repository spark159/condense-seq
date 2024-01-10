import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import seaborn as sns
import copy
import math
from scipy import signal

def read_input(fname):
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        if line.startswith('@'):
            fields = line[1:].split()
            for i in range(len(fields)):
                try:
                    fields[i] = float(fields[i])
                except:
                    pass
            data = {}
            for field in fields:
                data[field] = []
            continue
        cols = line.split()
        for i in range(len(cols)):
            field = fields[i]
            value = cols[i]
            try:
                value = float(value)
            except:
                pass
            data[field].append(value)
    return data

data = read_input("phase_input.csv")

X = data['Spermine'] * 4
#X = list(np.log(data['Spermine'])) * 4
Y = [12.5]*16 + [50]*16 + [250]*16 + [500]*16
Z = data[12.5] + data[50] + data[250] + data[500]

fig = plt.figure()
plt.scatter(X,Y,c=Z)
#plt.xscale("log")
#plt.xlim([0.005, 1.5])
plt.xlim([-0.05, 1.05])
plt.ylim([-50, 550])
plt.xlabel('Spermine (mM)')
plt.ylabel('DNA (ng/ul)')
plt.colorbar()
plt.savefig("phase_map.png", bbox_inches='tight')
#plt.show()
plt.close()
