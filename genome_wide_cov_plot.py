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

#path = "/home/spark159/../../media/spark159/sw/sp_spd_tests_detail/"
path = ""
field_ID_value = load_file.read_tabular_file (path + "H1_NCP_sp_chr1_cov.cn", mode='col', jump=1000)

ID_pos = field_ID_value['PhysicalPosition']
ID_cov1 = field_ID_value["work/2021_06_07_H1_sp_detail/H1-NCP-sp-0"]
ID_cov2 = field_ID_value["work/2021_06_07_H1_sp_detail/H1-NCP-sp-4"]
ID_cov3 = field_ID_value["work/2021_06_07_H1_sp_detail/H1-NCP-sp-8"]

ID_cov_list = [ID_cov1, ID_cov2, ID_cov3]
graphics.draw_along_genome (ID_pos, ID_cov_list, win=1000, labels=['H1-NCP-sp-0', 'H1-NCP-sp-4', 'H1-NCP-sp8'], ylabel='Coverage (X)', ylim=[-0.5,40], title='Chr1 Coverage', scatt=False, note="cov")
