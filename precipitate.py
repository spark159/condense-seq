import numpy as np
import math
import copy
import matplotlib.pyplot as plt

N = 100
C = 20
p_list = np.linspace(0.0, 1.0, 100)

fract_list = []
for p in p_list:
    fract = 1.0
    temp = 0.0
    for a in range(C, N+1):
        temp += a*(p**a)*((1-p)**(N-a))*float(math.factorial(N))/(math.factorial(N-a)*math.factorial(a))
    fract -= temp/N
    fract_list.append(fract)

fig = plt.figure()
plt.plot(p_list, fract_list)
plt.show()
plt.close()

