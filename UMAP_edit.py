import sys
import copy
import random
import pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sparse
#from sklearn.decomposition import NMF
import umap
import umap.plot

def read_tabular_file (fname,
                       mode='row',
                       delim='\t',
                       header=True,
                       rowID=True,
                       jump=None):
    if rowID:
        col_st = 1
    else:
        col_st = 0
        
    ID_field_value = {}
    First = True
    counter = -1
    for line in open(fname):
        cols = line.strip().split(delim)
        if First and header:
            field_names = cols[col_st:]
            First = False
            continue
        elif First and not header:
            field_names = range(len(cols[col_st:]))
            First = False
            pass

        counter += 1
        if jump and counter % jump != 0:
            continue

        if rowID:
            ID = cols[0]
        else:
            ID = counter

        if ID not in ID_field_value:
            ID_field_value[ID] = {}

        cols = cols[col_st:]
        #print cols
        for i in range(len(cols)):
            field = field_names[i]
            try:
                value = float(cols[i])
            except:
                value = cols[i]
            if value == 'NA':
                value = np.NaN
            if field not in ID_field_value[ID]:
                ID_field_value[ID][field] = value
            else:
                if type(ID_field_value[ID][field]) != list:
                    ID_field_value[ID][field] = [ID_field_value[ID][field]]
                ID_field_value[ID][field].append(value)

    if mode == 'row':
        return ID_field_value

    if mode == 'col' or mode == 'both':
        field_ID_value = {}
        for ID in ID_field_value:
            field_value = ID_field_value[ID]
            for field in field_value:
                value = field_value[field]
                if field not in field_ID_value:
                    field_ID_value[field] = {}
                field_ID_value[field][ID] = value

    if mode == 'col':
        return field_ID_value

    if mode == 'both':
        return ID_field_value, field_ID_value



# load annotation data
print("Data reading start", file=sys.stderr)

path = "/home/spark159/../../storage/"

#fname = 'H1_NCP_sp_10kb_chip-only_anot.txt'
#field = "H1-NCP-sp-8"

#fname = 'H1_NCP_HP1a_10kb_chip-only_anot.txt'
#field = "H1-NCP-HP1a-3"

fname = 'H1_NCP_Ki67_10kb_chip-only_anot.txt'
field = "H1-NCP-Ki67-4"



field_ID_value = read_tabular_file(path+fname, mode='col')
ID_score = field_ID_value[field]


# set field names
domain_fields = {'NSpeckle':['SON'],
                 'Trnx':['POLR2A', 'POLR2AphosphoS5', 'H3K9ac', 'H3K4me3', 'H3K27ac'],
                 'Polycomb':['CBX8', 'EZH2', 'RNF2', 'SUZ12'],
                 'Hetero':['H3K9me3', 'CBX5'],
                 'Nucleolus':['Nucleolar'],
                 'Lamin':['LaminB1'],
                 'Compartment':['eigen'],
                 'other':['ATcontent']}

domains = ['NSpeckle', 'Trnx', 'Polycomb', 'Hetero', 'Nucleolus', 'Lamin', 'Compartment', 'other']
fields = []
for domain in domains:
    fields += domain_fields[domain]


    
# make input data
#IDs = ID_score.keys()
IDs = random.sample(sorted(ID_score.keys()), 10000)

X = [[] for i in range(len(IDs))]
for field in fields:
    values = [field_ID_value[field][ID] for ID in IDs]

    # nonegative rescaling
    #min_value = min(values)
    #max_value = max(values)
    #for i in range(len(IDs)):
        #value = values[i]
        #re_value = float(value-min_value)/max_value
        #X[i].append(re_value)

    # standardization
    mean = np.mean(values)
    std = np.std(values)
    for i in range(len(IDs)):
        value = values[i]
        if std > 0:
            re_value = float(value-mean)/std
        else:
            re_value = value
        X[i].append(re_value)


X = sparse.csr_matrix(X)
values = np.asarray([ID_score[ID] for ID in IDs])

del field_ID_value
print("Data reading is done", file=sys.stderr)


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
    with open("H_doamin.pickle", "wb") as f:
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


newcIDs = []
newX = []
for i in range(len(IDs)):
    ID = IDs[i]
    newcID = ID_newcID[ID]
    newX


# UMAP embedding
print("UMAP start", file=sys.stderr)
mapper = umap.UMAP(metric='cosine',
                   random_state=42,
                   n_neighbors=15,
                   low_memory=True,
                   verbose=True).fit(X)
print("UMAP done", file=sys.stderr)

fig = plt.figure()
#umap.plot.points(mapper, labels=cIDs, cmap='jet')
umap.plot.points(mapper, values=values, cmap='jet')
plt.savefig("UMAP.png", dpi=1000, bbox_inches='tight')
#plt.show()
plt.close()

