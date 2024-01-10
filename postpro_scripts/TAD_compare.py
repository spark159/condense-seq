import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
import random
import matplotlib.cm as cm
from scipy.optimize import curve_fit
from sklearn import linear_model
from scipy.special import expit

# read annotation file
#path = "./data/"
path = ""
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+"H1_NCP_sp_chr1_anot.cn")
ID_score1 = name_ID_value["work/2021_06_07_H1_sp_detail/H1-NCP-sp-4"]
ID_score2 = name_ID_value["work/2021_06_07_H1_sp_detail/H1-NCP-sp-8"]
#ID_seq = name_ID_value['Sequence']
ID_AT = name_ID_value['ATcontent']
#ID_CG = name_ID_value['CpGNumber']
#ID_me = name_ID_value['meGCNumber']


# read TAD domain annotation file
def read_TADfile (fname):
    id = 0
    TAD_interval = {}
    for line in open(fname):
        if line.startswith('#'):
            continue
        cols = line.strip().split()
        chr, start, end = cols[:3]
        start, end = int(start), int(end)
        TAD_interval[id] = (start, end)
        id +=1
    return TAD_interval

#TAD_interval = read_TADfile(path + "WT_contact_domain_list/5000_blocks.bedpe"
TAD_interval = read_TADfile(path + "TAD_chr1_5kb.bedpe")

TADinterval_dict = Interval_dict.double_hash(TAD_interval, 100000, 250000000)

OutTAD_IDs, InTAD_IDs = [], []
for ID in ID_pos:
    pos = ID_pos[ID]
    finds = TADinterval_dict.find(pos)
    if len(finds) <= 0:
        OutTAD_IDs.append(ID)
    else:
        InTAD_IDs.append(ID)

#box_data = [[ID_score2[ID] for ID in OutTAD_IDs], [ID_score2[ID] for ID in InTAD_IDs]]
box_data = [[ID_AT[ID] for ID in OutTAD_IDs], [ID_AT[ID] for ID in InTAD_IDs]]
fig = plt.figure(figsize=(2.8, 2.4))
plt.boxplot(box_data, 0, "", positions=[1, 1.5])
#plt.xlabel('Partitions')
plt.ylabel('Condensability(A.U.)', fontsize=8)
plt.xticks([1,1.5], ['Outside of TAD', 'Inside of TAD'], fontsize=8)
plt.xlim([1-0.3,1.5+0.3])
plt.title('Single-nucleosome condensability', fontsize=8)
plt.gca().tick_params(axis='both', which='major', labelsize=5)
plt.gca().tick_params(axis='both', which='minor', labelsize=5)
plt.savefig('TADcompare_condensability.svg', format='svg', bbox_inches='tight')
#plt.savefig("ATcontent"+"_pbox.png")
#plt.show()
plt.close()



# read A/B compartment annotation file
def read_eigenfile (fname, bin_size=1000000):
    eigen_list = []
    interval_list = []
    i = 0
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        try:
            value = float(line)
        except:
            value = np.nan
        st = i*bin_size
        ed = (i+1)*bin_size
        eigen_list.append(value)
        interval_list.append((st,ed))
        i +=1
    return eigen_list, interval_list

#fnames = ['eigen_WT_50kb.txt', 'eigen_CohesinKO_50kb.txt']
fnames = ['eigen_H1_100kb.txt']
box_data_list = []
for fname in fnames:
    eigen_list, interval_list = read_eigenfile(path+fname, bin_size=100000)        
    eID_interval = {i:interval_list[i] for i in range(len(interval_list))}
    eID_dict = Interval_dict.double_hash(eID_interval, 100000, 250000000)
    A_IDs, B_IDs = [], []
    for ID in ID_pos:
        pos = ID_pos[ID]
        eIDs = eID_dict.find(pos)
        assert len(eIDs) <=1
        for eID in eIDs:
            evalue = eigen_list[eID]
            if evalue >= 0:
                A_IDs.append(ID)
            else:
                B_IDs.append(ID)

    box_data = [[ID_score2[ID] for ID in A_IDs], [ID_score2[ID] for ID in B_IDs]]
    box_data_list.append(box_data)

fig = plt.figure(figsize=(2.8, 2.4))
plt.boxplot(box_data_list[0], 0, "", positions=[1, 1.5])
#plt.xlabel('Partitions')
plt.ylabel('Condensability(A.U.)', fontsize=8)
plt.xticks([1,1.5], ['A-compartment', 'B-compartment'], fontsize=8)
plt.xlim([1-0.3,1.5+0.3])
plt.title('Single-nucleosome condensability', fontsize=8)
plt.gca().tick_params(axis='both', which='major', labelsize=5)
plt.gca().tick_params(axis='both', which='minor', labelsize=5)
plt.savefig('ABcompare_condensability.svg', format='svg', bbox_inches='tight')
#plt.savefig("ATcontent"+"_pbox.png")
#plt.show()
plt.close()


fig = plt.figure(figsize=(2.8, 2.4))
#box1 = plt.boxplot([box_data_list[0][0],box_data_list[1][0]], 0, "", positions=[1-0.1, 1.75-0.1], patch_artist=True)
box1 = plt.boxplot([box_data_list[0][0]], 0, "", positions=[1-0.15], patch_artist=True)
for patch in box1['boxes']:
    patch.set_facecolor('tab:red')
#box2 = plt.boxplot([box_data_list[0][1],box_data_list[1][1]], 0, "", positions=[1+0.1, 1.75+0.1], patch_artist=True)
box2 = plt.boxplot([box_data_list[0][1]], 0, "", positions=[1+0.15], patch_artist=True)
for patch in box2['boxes']:
    patch.set_facecolor('tab:blue')
#plt.xlabel('Partitions')
#plt.xticks([1,1.75],['WT', 'Cohesin KO'])
plt.xlim([0,2])
plt.xticks([],[])
plt.ylabel('Condensability(A.U.)', fontsize=8)
plt.title('Single-nucleosome condensability', fontsize=8)
#plt.title('')
#plt.savefig("ATcontent"+"_pbox.png")
plt.legend([box1["boxes"][0], box2["boxes"][0]], ['A compartment', 'B compartment'], loc='upper right', fontsize=8)
#plt.savefig('ABcompare_condensability.svg', format='svg', bbox_inches='tight')
#plt.show()
plt.close()

# read A/B subcompartment annotation file
def read_subcompartment (fname, chr_choice='chr1'):
    subID_interval = {}
    i = 0
    for line in open(fname):
        cols = line.strip().split()
        chr, st, ed, sub = cols[:4]
        if chr != chr_choice:
            continue
        st, ed = int(st), int(ed)
        subID = str(i) + ':' + sub
        subID_interval[subID] = (st, ed)
        i +=1
    return subID_interval

subID_interval = read_subcompartment (path+"GSE63525_GM12878_subcompartments.bed")
subID_dict = Interval_dict.double_hash(subID_interval, 100000, 250000000)
sub_IDs = {}
for ID in ID_pos:
    pos = ID_pos[ID]
    subIDs = subID_dict.find(pos)
    assert len(subIDs) <=1
    for subID in subIDs:
        _, sub = subID.split(':')
        if sub not in sub_IDs:
            sub_IDs[sub] = []
        sub_IDs[sub].append(ID)

sub_list = ['A1', 'A2', 'B1', 'B2', 'B3', 'NA']
box_data = []
for sub in sub_list:
    data = [ID_score1[ID] for ID in sub_IDs[sub]]
    box_data.append(data)

fig = plt.figure()
plt.boxplot(box_data, 0, "")
plt.show()
plt.close()

