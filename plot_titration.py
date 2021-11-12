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

# H1 ionic condensing agent
#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/H1_NCP_spermidine_new.csv")
#conc_list2, mean_list2, std_list2 = read_data("/home/spark159/Projects/condense-seq/H1_DNA_spermidine.csv")
#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/H1_NCP_PEG.csv")
#conc_list2, mean_list2, std_list2 = read_data("/home/spark159/Projects/condense-seq/H1_DNA_PEG.csv")
#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/H1_NCP_CoHex.csv")
#conc_list2, mean_list2, std_list2 = read_data("/home/spark159/Projects/condense-seq/H1_DNA_CoHex.csv")
#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/H1_NCP_Mg.csv")
#conc_list2, mean_list2, std_list2 = read_data("/home/spark159/Projects/condense-seq/H1_NCP_Ca2.csv")
#conc_list1, mean_list1, std_list1 = read_data("NCP_PstI.csv")
#conc_list2, mean_list2, std_list2 = read_data("NCPDNA_PstI.csv")

# H1 protein condensing agent
#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/H1_NCP_HP1a_with_PEG.csv")
#conc_list2, mean_list2, std_list2 = read_data("/home/spark159/Projects/condense-seq/H1_DNA_HP1a_with_PEG.csv")

#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/H1_NCP_HP1b+tSUV39H1_with_PEG.csv")
#conc_list2, mean_list2, std_list2 = read_data("/home/spark159/Projects/condense-seq/H1_DNA_HP1b+tSUV39H1_with_PEG.csv")

#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/H1_NCP_HP1b+tTRIM28_with_PEG.csv")
#conc_list2, mean_list2, std_list2 = read_data("/home/spark159/Projects/condense-seq/H1_DNA_HP1b+tTRIM28_with_PEG.csv")

#conc_list2, mean_list2, std_list2 = read_data("/home/spark159/Projects/condense-seq/H1_NCP_HP1b+tSUV39H1_with_PEG.csv")
#conc_list3, mean_list3, std_list3 = read_data("/home/spark159/Projects/condense-seq/H1_NCP_HP1b+tTRIM28_with_PEG.csv")

#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/H1_NCP_Ki67_with_PEG.csv")
#conc_list2, mean_list2, std_list2 = read_data("/home/spark159/Projects/condense-seq/H1_DNA_Ki67_with_PEG.csv")

#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/H1_NCP_LKH1_redoredo.csv")
#conc_list2, mean_list2, std_list2 = read_data("/home/spark159/Projects/condense-seq/H1_DNA_LKH1.csv")

conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/H1_NCP_WT_FUS_with_PEG.csv")
#conc_list2, mean_list2, std_list2 = read_data("/home/spark159/Projects/condense-seq/H1_DNA_Ki67_with_PEG.csv")

# Synthetic mono-nucleosome library
#conc_list, mean_list, std_list = read_data("/home/spark159/Projects/condense-seq/BC1_6_NCP_spermine.csv")
#conc_list, mean_list, std_list = read_data("/home/spark159/Projects/condense-seq/BC1_6_NCP_HP1a_with_PEG.csv")
#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/BC1_6_NCP_HP1a_with_PEG_1%_spike-in.csv")

#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/Oncolib_spermine.csv")
#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/Oncolib_spermidine.csv")
#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/Oncolib_CoHex.csv")
#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/Oncolib_PEG.csv")

#PTM library
#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/PTMlib_spermine_corrected.csv")
#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/PTMlib_spermidine_corrected.csv")
#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/PTMlib_CoHex_corrected.csv")
#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/PTMlib_PEG.csv")
#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/PTMlib_HP1alpha.csv")
#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/PTMlib_spermine_filter.csv")
#conc_list1, mean_list1, std_list1 = read_data("/home/spark159/Projects/condense-seq/PTMlib_spermidine_filter.csv")


del conc_list1[0]
del mean_list1[0]
del std_list1[0]

#del conc_list2[0]
#del mean_list2[0]
#del std_list2[0]

#del conc_list3[0]
#del mean_list3[0]
#del std_list3[0]

#del conc_list2[0:3]
#del mean_list2[0:3]
#del std_list2[0:3]


fig = plt.figure()
#plt.plot(conc_list1, mean_list1, 'bo-', label='NCP')
#plt.plot(conc_list1, mean_list1, 'bo-', label='Mg2+')
#plt.errorbar(conc_list1, mean_list1, yerr=std_list1, fmt='.', alpha=0.8)

plt.plot(conc_list1, mean_list1, 'bo-', label='NCP')
plt.errorbar(conc_list1, mean_list1, yerr=std_list1, fmt='.', alpha=0.8)
#plt.plot(conc_list2, mean_list2, 'ro-', label='DNA')
#plt.errorbar(conc_list2, mean_list2, yerr=std_list2, fmt='.', alpha=0.8)
#plt.plot(conc_list3, mean_list3, 'bo-', label='HP1$\\beta$ + tTRIM28')
#plt.errorbar(conc_list3, mean_list3, yerr=std_list3, color='b', fmt='.', alpha=0.8)

#plt.xlabel("Spermine concentration (mM)")
#plt.xlabel("Spermidine concentration (mM)")
#plt.xlabel("HP1$\\alpha$ concentration ($\\mu$M)")
#plt.xlabel("CoHex concentration (mM)")
#plt.xlabel("PEG 8000 (%)")
#plt.xlabel("Concentration (mM)")
#plt.xlabel("Ki-67 concentration ($\\mu$M)")
#plt.xlabel("HP1 concentration ($\\mu$M)")
#plt.xlabel("Linker histone H1.0 concentraion ($\\mu$M)")
plt.xlabel("WT FUS concentraion ($\\mu$M)")
plt.ylabel("Soluble fraction")
#plt.xscale("log")
plt.xscale("log", basex=2)
plt.ylim([0,1])
#plt.ylim([-0.1,1.1])
#plt.ylim([-0.07,1.2])
#plt.ylim([0,1])
#plt.ylim([-0.07,1.1])
#plt.ylim([0,1.2])
#plt.ylim([-0.05,1.5])
#plt.ylim([-0.05, 0.9])
#plt.ylim([0.2,1.1])
#plt.xscale('log', basex=2)
#plt.title("BC1~6 NCP mix 1% spike-in")
#plt.title("BC1~6 NCP mix HP1$\\alpha$")
#plt.title("H1 spermidine($3+$)")
#plt.title("PstI treatment test")
#plt.title("H1 CoHex$3+$")
#plt.title("H1 PEG 8000")
#plt.title("H1 NCP Mg2+/Ca2+")
#plt.title("H1 Ki-67")
#plt.title("H1 HP1$\\alpha$")
#plt.title("PTM MN library 1% spike-in")
#plt.title("H1 Linker histone H1.0")
plt.title("H1 WT FUS")
#plt.title("H1 HP1$\\beta+$tSUV39H1")
#plt.title("H1 HP1$\\beta+$tTRIM28")
plt.legend()
#plt.tight_layout()
plt.show()
plt.close()
