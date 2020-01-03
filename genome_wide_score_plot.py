import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
from sklearn import linear_model
import random

path = "/home/spark159/../../media/spark159/sw/sp_spd_tests_detail/"
field_ID_value = load_file.read_tabular_file (path + "hg19_chr1_167win25step_anot.cn", mode='col', jump=10)

ID_pos = field_ID_value['PhysicalPosition']
ID_AT = field_ID_value['ATcontent']
ID_k27ac = field_ID_value['k27ac']
ID_score1 = field_ID_value["data/sp_spd_tests_detail/sp7"]
ID_score2 = field_ID_value["data/sp_spd_tests_detail/sp8"]

for field in field_ID_value:
    if field.startswith("data/"):
        continue
    if field in ['PhysicalPosition', 'SNP', 'Chromosome']:
        continue
    print field
    ID_value = field_ID_value[field]
    graphics.draw_along_genome_pair (ID_pos, ID_score1, ID_value, win=10000, ylabel1='NCP sp6 score', ylabel2=field, xlabel="Genomic coordinate (bp)", title='Chromosome1', note='Sp6 VS ' + field)
    graphics.draw_along_genome_pair (ID_pos, ID_score2, ID_value, win=10000, ylabel1='NCP sp7 score', ylabel2=field, xlabel="Genomic coordinate (bp)", title='Chromosome1', note='Sp7 VS ' + field)


#ID_cov_list = [ID_cov1, ID_cov2, ID_cov3]
#graphics.draw_along_genome (ID_pos, ID_cov_list, win=1000, labels=['Input', 'NCP Sp6', 'NCP Sp7'], ylabel='Coverage (X)', ylim=[0,100], title='Chr1 Coverage', scatt=False, note="cov")

