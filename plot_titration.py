import numpy as np
import matplotlib.pyplot as plt

def read_data (fname):
    conc_list = []
    mean_list = []
    std_list = []
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        if not line:
            continue
        cols = line.strip().split()
        conc, mean, std = float(cols[0]), float(cols[1]), float(cols[2])
        conc_list.append(conc)
        mean_list.append(mean)
        std_list.append(std)
    return conc_list, mean_list, std_list
conc_list, mean_list, std_list = read_data("/home/spark159/script/condense-seq/BC1_6_NCP_spermine.csv")


fig = plt.figure()
plt.plot(conc_list, mean_list, 'bo-')
plt.errorbar(conc_list, mean_list, yerr=std_list, fmt='.')
plt.xlabel("Spermine concentration (mM)")
plt.ylabel("Soluble fraction")
plt.xscale("log")
plt.title("H1 NCP spermine(4+)")
plt.show()
plt.close()

