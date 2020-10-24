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

def all_path (N, states='ATCG'):
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output


def poly_score (seq, nts='AT', pos=False):
    num = []
    num_pos = {}
    i = 0
    while i < len(seq):
        if seq[i] in nts:
            nt = seq[i]
            count = 1
            j = i + 1
            while j < len(seq):
                if seq[j] != nt:
                    break
                count +=1
                j +=1
            num.append(count)
            if count not in num_pos:
                num_pos[count] = []
            num_pos[count].append(i)
            i = j
        else:
            i +=1
    if pos:
        return num_pos
    if len(num) == 0:
        return 0
    return max(num)

def get_dincount(seq, din=None):
    if din:
        count = 0
        for i in range(len(seq)-1):
            if seq[i:i+2].upper() == din.upper():
                count +=1
        return count
    din_count = {}
    seq = seq.upper()
    for i in range(len(seq)-1):
        din = seq[i:i+2]
        if 'N' in din:
            continue
        if din not in din_count:
            din_count[din] = 0
        din_count[din] += 1
    return din_count


# load annotation data
path = './data/'

ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+"hg19_chr1_NCP_ics_anot.cn")
ID_score1 = name_ID_value['data/sp_spd_tests_detail/sp7']
ID_seq = name_ID_value['Sequence']
ID_AT = name_ID_value['ATcontent']
ID_CpG = name_ID_value['CpGNumber']
ID_me = name_ID_value['meGCNumber']

IDs = sorted(ID_pos.keys())

#for ID in ID_AT:
#    ID_AT[ID] = ID_AT[ID]*100

ID_mefrac = {}
for ID in ID_CpG:
    CpG = ID_CpG[ID]
    if CpG <= 0:
        #ID_mefrac[ID] = np.NaN
        ID_mefrac[ID] = 0.5
        continue
    me = ID_me[ID]
    mefrac = float(me) / (2*CpG)
    ID_mefrac[ID] = mefrac


# get other sequence features
ID_polyAT, ID_polyGC = {}, {}
for ID in ID_seq:
    seq = ID_seq[ID]
    num_pos = poly_score(seq, nts='AT', pos=True)
    score = 0.0
    for num, pos in num_pos.items():
        if num < 3:
            continue
        score += len(pos)*(num)
        #score += len(pos)*(num**2)
    ID_polyAT[ID] = score
    num_pos = poly_score(seq, nts='GC', pos=True)
    score = 0.0
    for num, pos in num_pos.items():
        if num < 3:
            continue
        score += len(pos)*(num)
        #score += len(pos)*(num**2)
    ID_polyGC[ID] = score

ID_TA = {}
for ID in ID_seq:
    seq = ID_seq[ID]
    count = get_dincount(seq, din="TA")
    ID_TA[ID] = count

# collect all features
#names = ["AT content", "Poly-G", "CpG count", "TpA count", "meCpG density", "k4me3", "k27ac", "k9ac", "k36me3", "k9me2", "k9me3", "k27me3"]
#ID_value_list = [ID_AT, ID_polyGC, ID_CpG, ID_TA, ID_mefrac, name_ID_value['k4me3'], name_ID_value['k27ac'], name_ID_value['k9ac'], name_ID_value['k36me3_2'], name_ID_value['k9me2_2'], name_ID_value['k9me3_2'], name_ID_value['k27me3a_2']]

#names = ["AT content", "meCpG number", "k4me3", "k27ac", "k9ac", "k36me3", "k9me2", "k9me3", "k27me3"]
#ID_value_list = [ID_AT, ID_CpG, ID_me, name_ID_value['k4me3'], name_ID_value['k27ac'], name_ID_value['k9ac'], name_ID_value['k36me3_2'], name_ID_value['k9me2_2'], name_ID_value['k9me3_2'], name_ID_value['k27me3a_2']]

names = ["AT content", "meCpG density", "k4me3", "k27ac", "k9ac", "k36me3", "k9me2", "k9me3", "k27me3"]
ID_value_list = [ID_AT, ID_mefrac, name_ID_value['k4me3'], name_ID_value['k27ac'], name_ID_value['k9ac'], name_ID_value['k36me3_2'], name_ID_value['k9me2_2'], name_ID_value['k9me3_2'], name_ID_value['k27me3a_2']]

max_value_list = [float(max(ID_value.values())) for ID_value in ID_value_list] # for rescaling feature

ID_state = {}
for ID in IDs:
    state = []
    for i in range(len(ID_value_list)):
        ID_value, max_value = ID_value_list[i], max_value_list[i]
        state.append(ID_value[ID]) 
        #state.append(ID_value[ID]/max_value) # rescaling the feature
    ID_state[ID] = state

print "Data reading is done"


# None-negative matrix factorization
class_num = 7

try:
    with open("W.pickle", "rb") as f:
        W = pickle.load(f)
    with open("H.pickle", "rb") as f:
        H = pickle.load(f)

except:
    X = []
    for ID in IDs:
        X.append(ID_state[ID])
    X = sparse.csr_matrix(X)

    print "NMF start"
    model = NMF(n_components=class_num, init='random', random_state=0, verbose=True)
    W = model.fit_transform(X)
    H = model.components_
    print "NMF is done"

    with open("W.pickle", "wb") as f:
        pickle.dump(W, f)
    with open("H.pickle", "wb") as f:
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

# plot property matrix
img_list = []
for i in range(len(names)):
    img = np.zeros((class_num, len(names)))
    img[:] = np.nan
    for k in range(class_num):
        img[k][i] = H[k][i]
    img_list.append(img)

fig = plt.figure()

#cmap_list = ['hot', 'cool', 'YlGn', 'OrRd', 'Purples', 'magma', 'Greys', 'Greys', 'jet']
#cmap_list = ['Reds']*len(img_list)
cmap_list = ['Reds', 'Greens', 'Blues', 'Purples', 'Oranges', 'Greens', 'Blues', 'Purples', 'Oranges', 'Greens', 'Blues', 'Greys']
for img, cmap in zip(img_list, cmap_list):
    plt.imshow(img, cmap=cmap)
#plt.imshow(H)
for i in range(len(H)):
    for j in range(len(H[i])):
        mean, std = np.mean(H[:,j]), np.std(H[:,j])
        if H[i][j] > mean + std:
            color = 'white'
        else:
            color = 'black'
        plt.text(j, i, str(round(H[i][j], 2)), ha="center", va="center", fontsize=10, color=color)
plt.yticks(range(len(cID_prog)), range(1, len(cID_prog)+1))
plt.xticks(range(len(names)), names, rotation=70)
#cbar = plt.colorbar()
#cbar.ax.set_ylabel('$log_2$', rotation=-90, va="bottom")
plt.ylabel("Property class")
#plt.xlabel("Features")
plt.title("NMF property matrix")
plt.tight_layout()
plt.show()
plt.close()

# plot condensability of each property class
cID_scores = [[] for i in range(len(H))]

for i in range(len(cID_IDs)):
    for ID in cID_IDs[i]:
        score = ID_score1[ID]
        cID_scores[i].append(score)

fig = plt.figure()
plt.boxplot(cID_scores, 0, "")
plt.title("Condensability by property class")
plt.xlabel('Property class')
plt.ylabel('Condensability (A.U.)')
#plt.xticks(range(1, len(cID_scores)+1), range(1, len(cID_scores)+1))
plt.tight_layout()
#plt.savefig("anatomy_pbox.png")
plt.show()
plt.close()



