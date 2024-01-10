import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import seaborn as sns
import copy
import math
from scipy import signal

A = 10
B = 5

X, Y = [], []
dZ, uZ = [], []
for x in range(1,10):
    for y in range(1,10):
      X.append(x)
      Y.append(y)
      dZ.append(x*y*min(1.0/x, 1.0/y))
      uZ.append(x*y*min(float(A)/x, float(B)/y))

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(X, Y, dZ)
ax.scatter(X, Y, uZ)
plt.xlabel("x")
plt.ylabel("y")
plt.show()
