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
field_ID_value = load_file.read_tabular_file (path + "hg19_chr1_cov.cn", mode='col', jump=1000)

ID_pos = field_ID_value['PhysicalPosition']
ID_cov1 = field_ID_value["data/sp_spd_tests_detail/sp1"]
ID_cov2 = field_ID_value["data/sp_spd_tests_detail/sp7"]
ID_cov3 = field_ID_value["data/sp_spd_tests_detail/sp8"]

ID_cov_list = [ID_cov1, ID_cov2, ID_cov3]
graphics.draw_along_genome (ID_pos, ID_cov_list, win=1000, labels=['Input', 'NCP Sp6', 'NCP Sp7'], ylabel='Coverage (X)', ylim=[0,100], title='Chr1 Coverage', scatt=False, note="cov")

