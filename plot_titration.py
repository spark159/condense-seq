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
#conc_list, mean_list, std_list = read_data("/home/spark159/script/condense-seq/BC1_6_NCP_spermine.csv")
#conc_list, mean_list, std_list = read_data("/home/spark159/script/condense-seq/BC1_6_NCP_HP1a_with_PEG.csv")
conc_list, mean_list, std_list = read_data("/home/spark159/script/condense-seq/H1_NCP_spermidine.csv")

#del conc_list[0]
#del mean_list[0]
#del std_list[0]

fig = plt.figure()
plt.plot(conc_list, mean_list, 'bo-')
plt.errorbar(conc_list, mean_list, yerr=std_list, fmt='.')
plt.xlabel("Spermidine concentration (mM)")
#plt.xlabel("HP1$\\alpha$ concentration ($\\mu$M)")
plt.ylabel("Soluble fraction")
plt.xscale("log")
#plt.xscale('log', basex=2)
#plt.title("BC1~6 NCP mix spermine ($4+$)")
#plt.title("BC1~6 NCP mix HP1$\\alpha$")
plt.show()
plt.close()

