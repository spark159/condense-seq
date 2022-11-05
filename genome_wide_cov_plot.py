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
path = "/home/spark159/../../media/spark159/sw/"
field_ID_value = load_file.read_tabular_file (path + "mCD8T_inht-NCP_sp_chr1_cov.cn", mode='col', jump=1000)

ID_pos = field_ID_value['PhysicalPosition']
ID_cov1 = field_ID_value["/home/spark159/scratch4-tha4/sangwoo/MouseCD8Tcell_detail/mCD8T-inht-NCP-sp-0"]
ID_cov2 = field_ID_value["/home/spark159/scratch4-tha4/sangwoo/MouseCD8Tcell_detail/mCD8T-inht-NCP-sp-4"]
ID_cov3 = field_ID_value["/home/spark159/scratch4-tha4/sangwoo/MouseCD8Tcell_detail/mCD8T-inht-NCP-sp-8"]

ID_cov_list = [ID_cov1, ID_cov2, ID_cov3]
#ID_cov_list = [ID_cov1, ID_cov2]
#graphics.draw_along_genome (ID_pos, ID_cov_list, win=1000, labels=['GM-NCP-sp-0', 'GM-NCP-sp-4', 'GM-NCP-sp-8'], ylabel='Coverage (X)', ylim=[-0.5,50], title='Chr1 Coverage', scatt=False, note="cov")
#graphics.draw_along_genome (ID_pos, ID_cov_list, win=1000, labels=['H1-DNA-HP1a-0', 'H1-DNA-HP1a-3'], ylabel='Coverage (X)', ylim=[-0.5,50], title='Chr1 Coverage', scatt=False, note="cov")

graphics.draw_along_genome (ID_pos, ID_cov_list, win=1000, labels=['mCD8T:inht-NCP-sp-0', 'mCD8T:inht-NCP-sp-4', 'mCD8T:inht-NCP-sp-8'], ylabel='Coverage (X)', ylim=[-0.5,50], title='Chr1 Coverage', scatt=False, note="cov")



