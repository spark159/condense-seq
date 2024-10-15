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

path = ''
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path+"H1_NCP_sp_chr1_extended_anot.cn")
ID_score = name_ID_value['work/2021_06_07_H1_sp_detail/H1-NCP-sp-8']
name_ID_value['AT content'] = name_ID_value['ATcontent']
name_ID_value['meCpG density'] = get_density(name_ID_value['CNumber(CpG)'],
                                             name_ID_value['meCNumber(CpG)'])
name_ID_value['meCHG density'] = get_density(name_ID_value['CNumber(CHG)'],
                                             name_ID_value['meCNumber(CHG)'])
name_ID_value['meCHH density'] = get_density(name_ID_value['CNumber(CHH)'],
                                             name_ID_value['meCNumber(CHH)'])

ID_seq = name_ID_value['Sequence']
ID_polyGC = {}
for ID in ID_seq:
    seq = ID_seq[ID]
    num_pos = poly_score(seq.upper(), nts='GC', pos=True)
    mean_len, count = 0.0, 0.0
    for num, pos in num_pos.items():
        mean_len += len(pos)*num
        count += len(pos)
    ID_polyGC[ID] = mean_len/count

name_ID_value['poly-G/C length'] = ID_polyGC

del ID_seq
del ID_polyGC

names = ['AT content', 'poly-G/C length', 'meCpG density', 'meCHG density', 'meCHH density', 'H2AFZ', 'H2AK5ac', 'H2BK120ac', 'H2BK12ac', 'H2BK15ac', 'H2BK20ac', 'H2BK5ac', 'H3K14ac', 'H3K18ac', 'H3K23ac', 'H3K23me2', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K56ac', 'H3K79me1', 'H3K79me2', 'H3K9ac', 'H3K9me3', 'H4K20me1', 'H4K5ac', 'H4K8ac', 'H4K91ac']


#names = ['ATcontent', 'meCpGfrac', 'H2AFZ', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3k27ac', 'H3K9ac', 'H3K36me3', 'H3K9me3', 'H3K27me3', 'H4k20me1']
IDs = sorted(ID_pos.keys())
#IDs = random.sample(sorted(ID_pos.keys()), 1000)

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
    with open("W.pickle", "rb") as f:
        W = pickle.load(f)
    with open("H.pickle", "rb") as f:
        H = pickle.load(f)

except:
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
    
with open("NMF_sorted_cID.pickle", "wb") as f:
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
f = open("NMF_property_class.txt", 'w')
print >> f, "Class#" + "\t" + "\t".join(names)
for i in range(len(cID_prog)):
    print >> f, str(i+1) + "\t" + "\t".join([str(value) for value in cID_prog[i]])
f.close()

f = open("NMF_NCPClass.txt", 'w')
print >> f, "ID" + "\t" + "Class#"
for ID in ID_cID:
    print >> f, str(ID) + "\t" + str(ID_cID[ID]+1)
f.close()
