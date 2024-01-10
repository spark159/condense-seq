import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math

path = "./data/"
ratio = 0.181

ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+"hg19_chr1_NCP_ics_anot.cn")
ID_score = name_ID_value["data/sp_spd_tests_detail/sp7"]

f = open("chr1_condensability.txt", 'w')

print >> f, "%s\t%s\t%s\t%s" % ("ID", "Chromosome", "PhysicalPosition", "Condensability")

for ID in sorted(ID_pos.keys()):
    chr = ID_chr[ID]
    pos = ID_pos[ID]
    new_score = ID_score[ID] - np.log(ratio)
    print >> f, "%d\t%s\t%d\t%.5f" % (ID, chr, pos, new_score)

f.close()
