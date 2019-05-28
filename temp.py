import numpy as np
import matplotlib.pyplot as plt


def read(fname):
    count = 0
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        if First:
            total = [0.0] * len(cols[6:])
            First = False
            continue
        profile = cols[6:]
        for i in range(len(profile)):
            total[i] += float(profile[i])
        count +=1
    return [total[i]/count for i in range(len(total))]

mean_profile = read("data/sp1_chr1_profile.txt")

fig = plt.figure()
plt.plot(mean_profile)
plt.show()
