import numpy as np
import matplotlib.pyplot as plt

Y = [np.exp(-(i)**2) for i in range(-10,10)]
Y = np.asarray(Y)
for i in range(5):
    plt.plot(Y*i)
plt.show()
