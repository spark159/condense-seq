import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
from sklearn.decomposition import NMF
import scipy.sparse as sparse
import pickle
import random

def all_path (N, states='ATCG'):
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output

def get_density (ID_CNum, ID_meNum):
    ID_mefrac = {}
    for ID in ID_CNum:
        CNum = ID_CNum[ID]
        if CNum <= 0:
            ID_mefrac[ID] = 0.0
            continue
        meNum = ID_meNum[ID]
        mefrac = float(meNum) / (CNum)
        ID_mefrac[ID] = mefrac
    return ID_mefrac


# load annotation data
print "Data reading start"

path = "/home/spark159/../../storage/"

fname = 'H1_NCP_sp_10kb_chip-only_anot.txt'
field = "H1-NCP-sp-8"

#fname = 'H1_NCP_HP1a_10kb_chip-only_anot.txt'
#field = "H1-NCP-HP1a-3"

#fname = 'H1_NCP_Ki67_10kb_chip-only_anot.txt'
#field = "H1-NCP-Ki67-4"

name_ID_value = load_file.read_tabular_file(path+fname, mode='col')
ID_score = name_ID_value[field]


# set domain names
domain_names = {'NSpeckle':['SON'],
                 'Trnx':['POLR2A', 'POLR2AphosphoS5', 'H3K9ac', 'H3K4me3', 'H3K27ac'],
                 'Polycomb':['CBX8', 'EZH2', 'RNF2', 'SUZ12'],
                 'Hetero':['H3K9me3', 'CBX5'],
                 'Nucleolus':['Nucleolar'],
                 'Lamin':['LaminB1'],
                 'Compartment':['eigen'],
                 'other':['ATcontent']}

domains = ['NSpeckle', 'Trnx', 'Polycomb', 'Hetero', 'Nucleolus', 'Lamin', 'Compartment', 'other']

names = []
for domain in domains:
    names += domain_names[domain]

IDs = sorted(ID_score.keys())

#sys.exit(1)

X = [[] for i in range(len(IDs))]
for name in names:
    print name
    values = [name_ID_value[name][ID] for ID in IDs]
    min_value = min(values)
    max_value = max(values)
    for i in range(len(IDs)):
        value = values[i]
        re_value = float(value-min_value)/max_value
        X[i].append(re_value)
    del values
    del min_value
    del max_value
    
#with open("NMF_input.pickle", "wb") as f:
#    pickle.dump(X, f)

X = sparse.csr_matrix(X)
#values = np.asarray([ID_score[ID] for ID in IDs])

del ID_pos
del ID_chr
del name_ID_value

print "Data reading is done"


# None-negative matrix factorization
class_num = 10

try:
    with open("W_domain.pickle", "rb") as f:
        W = pickle.load(f)
    with open("H_domain.pickle", "rb") as f:
        H = pickle.load(f)

except:
    print "NMF start"
    model = NMF(n_components=class_num, init='random', random_state=0, verbose=True)
    W = model.fit_transform(X)
    H = model.components_
    print "NMF is done"

    with open("W_domain.pickle", "wb") as f:
        pickle.dump(W, f)
    with open("H_domain.pickle", "wb") as f:
        pickle.dump(H, f)

# post-analysis of NMF
cID_prog = []
for i in range(class_num):
    cID_prog.append(H[i])

ID_cID = {}
cID_IDs = [[] for i in range(class_num)]
for i in range(len(IDs)):
    ID = IDs[i]
    cID = np.argmax(W[i])
    ID_cID[ID] = cID
    cID_IDs[cID].append(ID)

cID_scores = [[] for i in range(class_num)]
for i in range(len(cID_IDs)):
    for ID in cID_IDs[i]:
        score = ID_score[ID]
        cID_scores[i].append(score)

# sort according to condensability
score_cID = sorted([(np.median(cID_scores[cID]), cID) for cID in range(len(cID_scores))])
cID_list = [cID for score, cID in score_cID]
cID_newcID = {cID_list[i]:i for i in range(len(cID_list))}
ID_newcID = {}
for ID in IDs:
    cID = ID_cID[ID]
    newcID = cID_newcID[cID]
    ID_newcID[ID] = newcID
    
with open("NMF_sorted_cID_domain.pickle", "wb") as f:
    pickle.dump(ID_newcID, f)


# plot property matrix
img_list = []
for i in range(len(names)):
    img = np.zeros((len(names), class_num))
    img[:] = np.nan
    for k in range(len(cID_list)):
        cID = cID_list[k]
        img[i][k] = cID_prog[cID][i]
    img_list.append(img)

cmap_list = ['Reds', 'Blues', 'Greens', 'Purples', 'Oranges',  'Greys']*10


fig = plt.figure(figsize=(5, 6))

for img, cmap in zip(img_list, cmap_list):
    plt.imshow(img, cmap=cmap, aspect='auto')

for i in range(len(cID_list)):
    cID = cID_list[i]
    for j in range(len(names)):
        mean, std = np.mean(H[:,j]), np.std(H[:,j])
        value = cID_prog[cID][j]
        if value > mean + std:
            color = 'white'
        else:
            color = 'black'
        plt.text(i, j, str(round(value, 2)), ha="center", va="center", fontsize=8, color=color)

plt.xticks(range(len(cID_prog)), range(1, len(cID_prog)+1), fontsize=8)
plt.yticks(range(len(names)), names, fontsize=8)
plt.xlabel("Property class", fontsize=10)
#plt.title("NMF property matrix", fontsize=12)
plt.savefig("NMF_property_matrix.svg", format='svg', bbox_inches='tight')
#plt.tight_layout()
#plt.show()
plt.close()

# plot condensability of each property class
fig = plt.figure(figsize=(5,3))
plt.boxplot([cID_scores[cID] for cID in cID_list], 0, "")
plt.title("Condensability by property class", fontsize=12)
#plt.xlabel('Property class', fontsize=10)
plt.ylabel('Condensability (A.U.)', fontsize=10)
plt.gca().tick_params(axis='both', which='major', labelsize=8)
plt.gca().tick_params(axis='both', which='minor', labelsize=8)
plt.savefig("condensability_by_class.svg", format='svg', bbox_inches='tight')
#plt.xticks(range(1, len(cID_scores)+1), range(1, len(cID_scores)+1))
#plt.tight_layout()
#plt.savefig("anatomy_pbox.png")
#plt.show()
plt.close()

# print NMF result
f = open("NMF_property_class_domain.txt", 'w')
print >> f, "Class#" + "\t" + "\t".join(names)
for i in range(len(cID_prog)):
    print >> f, str(i+1) + "\t" + "\t".join([str(value) for value in cID_prog[i]])
f.close()

f = open("NMF_NCPClass_domain.txt", 'w')
print >> f, "ID" + "\t" + "Class#"
for ID in ID_cID:
    print >> f, str(ID) + "\t" + str(ID_cID[ID]+1)
f.close()
