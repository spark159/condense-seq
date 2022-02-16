import sys
import copy
import random
import pickle

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sparse
from sklearn.decomposition import NMF
import umap
import umap.plot

def read_anot_file(fname, target_names=None, jump=None, num_max=sys.maxsize):
    ID_chr, ID_pos = {}, {}
    name_ID_value = {}
    First = True
    count = 0
    counter = -1
    for line in open(fname):
        if count > num_max:
            break
        cols = line.strip().split()
        if First:
            names = cols[3:]
            First = False
            continue
        counter += 1
        if jump and counter % jump != 0:
            continue
        ID, chr, pos = int(cols[0]), cols[1], int(cols[2])
        ID_chr[ID] = chr
        ID_pos[ID] = pos
        cols = cols[3:]
        for i in range(len(cols)):
            name = names[i]
            if target_names and name not in target_names:
                continue
            if name not in name_ID_value:
                name_ID_value[name] = {}
            assert ID not in name_ID_value[name]
            try:
                value = float(cols[i])
            except:
                value = cols[i]
            if value == 'NA':
                value = np.NaN
            name_ID_value[name][ID] = value
        count += 1
    return ID_chr, ID_pos, name_ID_value

"""
# load annotation data
print("Data reading start", file=sys.stderr)

path = ''
ID_chr, ID_pos, name_ID_value = read_anot_file(path+"H1_NCP_sp_chr1_anot.cn")
ID_score = name_ID_value['work/2021_06_07_H1_sp_detail/H1-NCP-sp-8']

ID_CpG = name_ID_value['CNumber(CpG)']
ID_meCpG = name_ID_value['meCNumber(CpG)']


ID_meCpGfrac = {}
for ID in ID_CpG:
    CpG = ID_CpG[ID]
    if CpG <= 0:
        #ID_meCpGfrac[ID] = np.NaN
        ID_meCpGfrac[ID] = 0.0
        continue
    meCpG = ID_meCpG[ID]
    meCpGfrac = float(meCpG) / (CpG)
    ID_meCpGfrac[ID] = meCpGfrac

name_ID_value['meCpGfrac'] = ID_meCpGfrac

del ID_CpG
del ID_meCpG
del ID_meCpGfrac

names = ['ATcontent', 'meCpGfrac', 'H2AFZ', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3k27ac', 'H3K9ac', 'H3K36me3', 'H3K9me3', 'H3K27me3', 'H4k20me1']
#IDs = sorted(ID_pos.keys())
IDs = random.sample(sorted(ID_pos.keys()), 10000)

#sys.exit(1)

X = [[] for i in range(len(IDs))]
for name in names:
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
    
X = sparse.csr_matrix(X)
values = np.asarray([ID_score[ID] for ID in IDs])

del ID_score
del ID_pos
del ID_chr
del name_ID_value

print("Data reading is done", file=sys.stderr)
"""

# NMF decomposition
class_num=12
print("NMF start", file=sys.stderr)
try:
    with open("W.pickle", "rb") as f:
        W = pickle.load(f)
    with open("H.pickle", "rb") as f:
        H = pickle.load(f)

except:
    model = NMF(n_components=class_num, init='random', random_state=0, verbose=True)
    W = model.fit_transform(X)
    H = model.components_

    with open("W.pickle", "wb") as f:
        pickle.dump(W, f)
    with open("H.pickle", "wb") as f:
        pickle.dump(H, f)
        
print("NMF end", file=sys.stderr)

cIDs = []
newX = []
for i in range(len(IDs)):
    cID = np.argmax(W[i])
    cIDs.append(cID)
    #total = sum(W[i])
    #newX.append([float(comp)/total for comp in W[i]])
    newX.append(W[i])
cIDs = np.asarray(cIDs)
newX = sparse.csr_matrix(newX)



# UMAP embedding
print("UMAP start", file=sys.stderr)
mapper = umap.UMAP(metric='cosine',
                   random_state=42,
                   n_neighbors=15,
                   low_memory=True,
                   verbose=True).fit(newX)
print("UMAP done", file=sys.stderr)

fig = plt.figure()
umap.plot.points(mapper, labels=cIDs, cmap='jet')
#umap.plot.points(mapper, values=values, cmap='jet')
plt.savefig("UMAP.png", dpi=1000, bbox_inches='tight')
#plt.show()
plt.close()

