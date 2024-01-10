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
with open("W.pickle", "rb") as f:
    W = sparse.csr_matrix(pickle.load(f, encoding='latin1'))
with open("NMF_sorted_cID.pickle", "rb") as f:
    ID_cID = pickle.load(f, encoding='latin1')
IDs = sorted(ID_cID.keys())
labels = np.asarray([ID_cID[ID]+1 for ID in IDs])

# UMAP embedding
print("UMAP start", file=sys.stderr)
mapper = umap.UMAP(metric='cosine',
                   random_state=42,
                   n_neighbors=15,
                   low_memory=True,
                   verbose=True).fit(W)
print("UMAP done", file=sys.stderr)

fig = plt.figure()
umap.plot.points(mapper, labels=labels, cmap='jet')
plt.savefig("UMAP.png", dpi=1000, bbox_inches='tight')
#plt.show()
plt.close()

