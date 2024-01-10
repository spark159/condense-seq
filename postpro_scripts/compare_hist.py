import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
import Interval_dict
from scipy import stats
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.cm as cm


path = '/home/spark159/../../media/spark159/sw/'

field_ID_value1 = load_file.read_tabular_file(path+"mCD8T_WT-NCP_sp_chr1_score.cn", mode='col')
field_ID_value2 = load_file.read_tabular_file(path+"mCD8T_inht-NCP_sp_chr1_score.cn", mode='col')
field_ID_value3 = load_file.read_tabular_file(path+"mCD8T_KO-NCP_sp_chr1_score.cn", mode='col')

ID_value1 = field_ID_value1['mCD8T-WT-NCP-sp-8']
ID_value2 = field_ID_value2['mCD8T-inht-NCP-sp-8']
ID_value3 = field_ID_value3['mCD8T-KO-NCP-sp-8']

fig = plt.figure()
plt.hist(ID_value1.values(), bins=1000, alpha=0.5, label='WT')
plt.hist(ID_value2.values(), bins=1000, alpha=0.5, label='+inht')
plt.hist(ID_value3.values(), bins=1000, alpha=0.5, label='KO')
plt.xlim([-1, 6])
plt.legend()
plt.savefig("comp8.png", bbox_inches='tight')
#plt.show()
plt.close()

IDs = list(set(ID_value1) & set(ID_value2) & set(ID_value3))

X, Y, Z = [], [], []
for ID in IDs:
    X.append(ID_value1[ID])
    Y.append(ID_value2[ID])
    Z.append(ID_value3[ID])

fig = plt.figure()
plt.plot(X, Y, ',')
plt.show()
plt.close()
