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

ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("hg19_chr1_171_everything_anot.cn")
ID_score1 = name_ID_value['work/condense_seq/sp9_hg19_chr1']
ID_score2 = name_ID_value['work/condense_seq/sp10_hg19_chr1']
ID_AT_raw = name_ID_value['ATcontent']
#ID_CpG = name_ID_value['CpGNumber']
#ID_me = name_ID_value['meGCNumber']

ID_chip_list, chip_names = [], []
for name in name_ID_value:
    if name.startswith('k'):
        ID_value = name_ID_value[name]
        chip_names.append(name)
        ID_chip_list.append(ID_value)

ID_AT = {}
for ID in ID_AT_raw:
    ID_AT[ID] = ID_AT_raw[ID]*100

new_ID_score2 = statis.neutralize_score_by_target(ID_score2, ID_AT)


graphics.draw_along_genome (ID_pos, [ID_AT, ID_score2], 1, labels = ["AT", "raw"], ylabel="")




"""
scores2 = ID_score2.values()
new_scores2 = new_ID_score2.values()
feature_list, target_list = [], []
for ID in ID_AT:
    AT = ID_AT[ID]
    score2 = ID_score2[ID]
    feature_list.append([AT])
    target_list.append([score2])

reg = linear_model.Ridge(alpha=0.5)
reg.fit (feature_list, target_list)
reg = reg
print reg.coef_
rsquare = reg.score(feature_list, target_list)
print "r-square: " + str(rsquare)

X = [[i] for i in range(100)]
Y = reg.predict(X)
Y = [y[0] for y in Y]


#graphics.Scatter_plot(ID_AT, ID_score2, note='sp10')

fig = plt.figure()
plt.hist(scores2, bins=200)
plt.show()

fig = plt.figure()
plt.hist(new_scores2, bins=200)
plt.show()

"""
