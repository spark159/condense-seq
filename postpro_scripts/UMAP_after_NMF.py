import sys
import copy
import random
import pickle

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sparse
import umap
import umap.plot


# laod NMF decomposition result
#with open("W.pickle", "rb") as f:
#    W = pickle.load(f, encoding='latin1')
with open("NMF_input.pickle", "rb") as f:
    W = pickle.load(f, encoding='latin1')
with open("NMF_sorted_cID.pickle", "rb") as f:
    ID_cID = pickle.load(f, encoding='latin1')

# random sampling for each cluster
scale = 100

cID_IDs = {}
for ID, cID in ID_cID.items():
    if cID not in cID_IDs:
        cID_IDs[cID] = []
    cID_IDs[cID].append(ID)

choosen_IDs = []
for cID in cID_IDs:
    IDs = cID_IDs[cID]
    sample_num = int(len(IDs)/scale)
    if sample_num < 10:
        #sample_num = min(10, len(IDs))
        choosen_IDs += IDs
        continue
    random.seed = 1234
    choosen_IDs += random.sample(IDs, sample_num)
    
choosen_IDs = sorted(choosen_IDs)
X = []
for ID in choosen_IDs:
    X.append(W[ID])
X = sparse.csr_matrix(X)
labels = np.asarray([ID_cID[ID]+1 for ID in choosen_IDs])

del W
del ID_cID
del cID_IDs
del choosen_IDs

# UMAP embedding
print("UMAP start", file=sys.stderr)
mapper = umap.UMAP(metric='cosine',
                   random_state=42,
                   n_neighbors=15,
                   n_epochs=500, 
                   low_memory=True,
                   verbose=True).fit(X)
print("UMAP done", file=sys.stderr)

fig = plt.figure()
umap.plot.points(mapper, labels=labels, background='black')
plt.savefig("UMAP.png", dpi=1000, bbox_inches='tight')
#plt.show()
plt.close()

