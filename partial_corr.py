import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
import Interval_dict
from sklearn import linear_model
from scipy.stats import norm
from matplotlib_venn import venn3
from npeet import entropy_estimators as ee

def all_path (N, states='ATCG'):
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output

def Amer_len(seq, pos=False):
    num = []
    num_pos = {}
    i = 0
    while i < len(seq):
        if seq[i] in 'AT':
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

ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("/home/spark159/../../media/spark159/sw/sp_spd_tests_detail/hg19_chr1_NCP_ics_anot.cn")
ID_score2 = name_ID_value['data/sp_spd_tests_detail/sp7']
ID_seq = name_ID_value['Sequence']
ID_AT = name_ID_value['ATcontent']
ID_CpG = name_ID_value['CpGNumber']
ID_me = name_ID_value['meGCNumber']

for ID in ID_AT:
    ID_AT[ID] = ID_AT[ID]*100

ID_Alen = {}
for ID in ID_seq:
    seq = ID_seq[ID]
    num_pos = Amer_len(seq, pos=True)
    score = 0.0
    for num, pos in num_pos.items():
        if num < 3:
            continue
        score += len(pos)*(num**2)
    ID_Alen[ID] = score

ID_TA = {}
for ID in ID_seq:
    seq = ID_seq[ID]
    count = get_dincount(seq, din="TA")
    ID_TA[ID] = count

ID_mefrac = {}
for ID in ID_CpG:
    CpG = ID_CpG[ID]
    if CpG <= 0:
        ID_mefrac[ID] = np.NaN
        continue
    me = ID_me[ID]
    mefrac = float(me) / (2*CpG)
    ID_mefrac[ID] = mefrac

# correlation and p-correlation matrix
names = ["Condensability", "AT content", "Poly-A", "CpG count", "TpA count", "meCpG density", "k4me3", "k27ac", "k9ac", "k36me3", "k9me2", "k9me3", "k27me3"]
ID_value_list = [ID_score2, ID_AT, ID_Alen, ID_CpG, ID_TA, ID_me, name_ID_value['k4me3'], name_ID_value['k27ac'], name_ID_value['k9ac'], name_ID_value['k36me3_2'], name_ID_value['k9me2_2'], name_ID_value['k9me3_2'], name_ID_value['k27me3a_2']]

IDs = ID_seq.keys()
values_list = []
for ID_value in ID_value_list:
    temp = [ ID_value[ID] for ID in IDs ]
    values_list.append(temp)
print "reading done"
print

print "correlation"
corr_matrix = [ [ np.nan for i in range(len(names)) ] for i in range(len(names)) ]
for i in range(len(names)-1):
    for j in range(i+1, len(names)):
        corr = statis.get_corr (values_list[i], values_list[j])
        corr_matrix[i][j] = corr
        corr_matrix[j][i] = corr
        print names[i], names[j], corr
print

fig = plt.figure(figsize=(15, 15))
im = plt.matshow(corr_matrix, cmap='bwr', fignum=1, vmin=-0.6, vmax=0.6)

for i in range(len(corr_matrix)):
    for j in range(len(corr_matrix[i])):
        plt.text(i, j, str(round(corr_matrix[i][j], 2)), ha="center", va="center", fontsize=10, color="k")

plt.xticks(range(len(names)), names, rotation=90, fontsize=15)
plt.yticks(range(len(names)), names, fontsize=15)
plt.colorbar(im,fraction=0.046, pad=0.04)
#plt.colorbar()
plt.savefig("corr_matrix.png", bbox_inches='tight')
#plt.show()
plt.close()

print "p-correlation"
pcorr_matrix = [ [ np.nan for i in range(len(names)) ] for i in range(len(names)) ]
for i in range(len(names)-1):
    for j in range(i+1, len(names)):
        Zs = values_list[:i] + values_list[i+1:j] + values_list[j+1:]
        assert len(Zs) == len(names)-2
        ID_newvalue1 = statis.neutralize (ID_value_list[i], Zs)
        ID_newvalue2 = statis.neutralize (ID_value_list[j], Zs)
        A, B = [], []
        for ID in ID_newvalue1:
            newvalue1 = ID_newvalue1[ID]
            newvalue2 = ID_newvalue2[ID]
            A.append(newvalue1)
            B.append(newvalue2)
        pcorr = statis.get_corr(A, B)
        #pcorr = statis.partial_corr (values_list[i], values_list[j], Zs)
        pcorr_matrix[i][j] = pcorr
        pcorr_matrix[j][i] = pcorr
        print names[i], names[j], pcorr
print

fig = plt.figure(figsize=(15, 15))
im = plt.matshow(pcorr_matrix, cmap='bwr', fignum=1, vmin=-0.6, vmax=0.6)

for i in range(len(pcorr_matrix)):
    for j in range(len(pcorr_matrix[i])):
        plt.text(i, j, str(round(pcorr_matrix[i][j], 2)), ha="center", va="center", fontsize=10, color="k")

plt.xticks(range(len(names)), names, rotation=90, fontsize=15)
plt.yticks(range(len(names)), names, fontsize=15)
#plt.colorbar()
plt.colorbar(im,fraction=0.046, pad=0.04)
plt.savefig("pcorr_matrix.png", bbox_inches='tight')
#plt.show()
plt.close()


"""
ptm_pcorr = {}
# CpG methylation dependence
print "CpG me"
# get partial correlation by recursive algorithm
X, Y, Z = [], [], []
IDs = ID_AT.keys()
for ID in IDs:
    X.append(ID_me[ID])
    Y.append(ID_score2[ID])
    Z.append(ID_CpG[ID])

partial_corr = statis.partial_corr (X, Y, [Z])
ptm_pcorr["CpGme"] = partial_corr
print partial_corr

# get partial correlation by linear regression
ID_newme = statis.neutralize (ID_me, [ID_CpG])
ID_newscore2 = statis.neutralize (ID_score2, [ID_CpG])

A, B = [], []
for ID in ID_newme:
    newme = ID_newme[ID]
    newscore2 = ID_newscore2[ID]
    A.append(newme)
    B.append(newscore2)

print statis.get_corr(A, B)

feature_list = [ [a] for a in A ]
target_list = [ [b] for b in B ]

reg = linear_model.Ridge(alpha=0.5)
reg.fit (feature_list, target_list)
Ypred = reg.predict(feature_list)

fig = plt.figure()
plt.plot(A, B, '.', alpha=0.01)
plt.plot([x[0] for x in feature_list], [y[0] for y in Ypred], '--')
plt.title("Partial-correlation " + "(" + "CpGme" + ")")
plt.xlabel("CpGme" + ' - CG count')
plt.ylabel('Condensability - CG count')
#plt.savefig("pcorr_" + "CpGme" + "_scatter.apng")
plt.show()
plt.close()

# Venn diagram representation of variance    
X = statis.standardize(X)
Y = statis.standardize(Y)
Z = statis.standardize(Z)

area3 = (statis.semi_partial_corr(X, Y, [Z])**2)
area5 = (statis.semi_partial_corr(Z, Y, [X])**2)
area6 = (statis.semi_partial_corr(Z, X, [Y])**2)
area7 = (statis.get_corr(X, Y)**2) - area3
area1 = 1 - (area3 + area7 + area5)
area2 = 1 - (area3 + area7 + area6)
area4 = 1 - (area5 + area6 + area7)

print np.sqrt(area3 / (area1 + area3))

area1 = round(area1, 2)
area2 = round(area2, 2)
area3 = round(area3, 2)
area4 = round(area4, 2)
area5 = round(area5, 2)
area6 = round(area6, 2)
area7 = round(area7, 2)

fig = plt.figure()
venn3(subsets = (area1, area2, area3, area4, area5, area6, area7), set_labels = ('Condensability', "CpG me", 'CG count'))
#plt.savefig("pcorr_" + "CpGme" + "_Venn.png")
plt.show()
plt.close()

print


# histone PTM dependence
for name in name_ID_value:
    if name in ['work/condense_seq/sp9_hg19_chr1', 'work/condense_seq/sp10_hg19_chr1', 'Sequence', 'ATcontent', 'CpGNumber', 'meGCNumber' ]:
        continue

    print name
    ID_value = name_ID_value[name]
    # get partial correlation by recursive algorithm
    X, Y, Z = [], [], []
    IDs = ID_AT.keys()
    for ID in IDs:
        X.append(ID_value[ID])
        Y.append(ID_score2[ID])
        Z.append(ID_AT[ID])

    partial_corr = statis.partial_corr (X, Y, [Z])
    ptm_pcorr[name] = partial_corr
    print partial_corr

    # get partial correlation by linear regression
    ID_newvalue = statis.neutralize (ID_value, [ID_AT])
    ID_newscore2 = statis.neutralize (ID_score2, [ID_AT])

    A, B = [], []
    for ID in ID_newvalue:
        newvalue = ID_newvalue[ID]
        newscore2 = ID_newscore2[ID]
        A.append(newvalue)
        B.append(newscore2)
    
    print statis.get_corr(A, B)

    feature_list = [ [a] for a in A ]
    target_list = [ [b] for b in B ]

    reg = linear_model.Ridge(alpha=0.5)
    reg.fit (feature_list, target_list)
    Ypred = reg.predict(feature_list)

    fig = plt.figure()
    plt.plot(A, B, '.', alpha=0.01)
    plt.plot([x[0] for x in feature_list], [y[0] for y in Ypred], '--')
    plt.title("Partial-correlation " + "(" + name + ")")
    plt.xlabel(name + ' - AT content')
    plt.ylabel('Condensability - AT content')
    #plt.savefig("pcorr_" + name + "_scatter.apng")
    plt.show()
    plt.close()

    # Venn diagram representation of variance    
    X = statis.standardize(X)
    Y = statis.standardize(Y)
    Z = statis.standardize(Z)

    area3 = (statis.semi_partial_corr(X, Y, [Z])**2)
    area5 = (statis.semi_partial_corr(Z, Y, [X])**2)
    area6 = (statis.semi_partial_corr(Z, X, [Y])**2)
    area7 = (statis.get_corr(X, Y)**2) - area3
    area1 = 1 - (area3 + area7 + area5)
    area2 = 1 - (area3 + area7 + area6)
    area4 = 1 - (area5 + area6 + area7)

    print np.sqrt(area3 / (area1 + area3))

    area1 = round(area1, 2)
    area2 = round(area2, 2)
    area3 = round(area3, 2)
    area4 = round(area4, 2)
    area5 = round(area5, 2)
    area6 = round(area6, 2)
    area7 = round(area7, 2)
    
    fig = plt.figure()
    venn3(subsets = (area1, area2, area3, area4, area5, area6, area7), set_labels = ('Condensability', name, 'AT content'))
    #plt.savefig("pcorr_" + name + "_Venn.png")
    plt.show()
    plt.close()

    print

# dinucleotide dependence
# dinucleotide counting
ID_din_count = {}
for ID in ID_seq:
    seq = ID_seq[ID]
    din_count = get_dincount(seq)
    ID_din_count[ID] = din_count
print "dinucleotide counting is done"
print

din_pcorr = {}
dY = []
all_din = all_path(2)
for din in all_din:
    ID_din = {}
    for ID in ID_seq:
        try:
            ID_din[ID] = ID_din_count[ID][din]
        except:
            ID_din[ID] = 0

    print din
    # get partial correlation by recursive algorithm
    X, Y, Z = [], [], []
    IDs = ID_AT.keys()
    for ID in IDs:
        X.append(ID_din[ID])
        Y.append(ID_score2[ID])
        Z.append(ID_AT[ID])

    partial_corr = statis.partial_corr (X, Y, [Z])
    din_pcorr[din] = partial_corr
    print partial_corr

    # get partial correlation by linear regression
    ID_newdin = statis.neutralize (ID_din, [ID_AT])
    ID_newscore2 = statis.neutralize (ID_score2, [ID_AT])

    A, B = [], []
    for ID in ID_newdin:
        newdin = ID_newdin[ID]
        newscore2 = ID_newscore2[ID]
        A.append(newdin)
        B.append(newscore2)
    
    print statis.get_corr(A, B)

    feature_list = [ [a] for a in A ]
    target_list = [ [b] for b in B ]

    reg = linear_model.Ridge(alpha=0.5)
    reg.fit (feature_list, target_list)
    Ypred = reg.predict(feature_list)

    fig = plt.figure()
    plt.plot(A, B, '.', alpha=0.01)
    plt.plot([x[0] for x in feature_list], [y[0] for y in Ypred], '--')
    plt.title("Partial-correlation " + "(" + din + ")")
    plt.xlabel(din + ' count - AT content')
    plt.ylabel('Condensability - AT content')
    plt.savefig("pcorr_" + din + "_scatter.png")
    #plt.show()
    plt.close()

    
    ## check statistical significance
    #def p_value (sample_size, Z_num, partial_corr):
    #    def Fisher_transform(rho):
    #        return 0.5*np.log(float(1+rho)/(1-rho))
    #    return 2*(1-norm.cdf((np.sqrt(sample_size-Z_num-3))*abs(Fisher_transform(partial_corr))))
#
#    print 'p-value: ', p_value(len(X), 1, partial_corr)
    
    
    # Venn diagram representation of variance    
    X = statis.standardize(X)
    Y = statis.standardize(Y)
    Z = statis.standardize(Z)

    area3 = (statis.semi_partial_corr(X, Y, [Z])**2)
    area5 = (statis.semi_partial_corr(Z, Y, [X])**2)
    area6 = (statis.semi_partial_corr(Z, X, [Y])**2)
    area7 = (statis.get_corr(X, Y)**2) - area3
    area1 = 1 - (area3 + area7 + area5)
    area2 = 1 - (area3 + area7 + area6)
    area4 = 1 - (area5 + area6 + area7)

    print np.sqrt(area3 / (area1 + area3))

    area1 = round(area1, 2)
    area2 = round(area2, 2)
    area3 = round(area3, 2)
    area4 = round(area4, 2)
    area5 = round(area5, 2)
    area6 = round(area6, 2)
    area7 = round(area7, 2)
    
    fig = plt.figure()
    venn3(subsets = (area1, area2, area3, area4, area5, area6, area7), set_labels = ('Condensability', din + " count", 'AT content'))
    plt.savefig("pcorr_" + din + "_Venn.png")
    #plt.show()
    plt.close()

    
   
    # Compute mutual information and Van diagram representation
    if len(dY) <=0:
        dY = statis.binning(Y, 1000) # binning continuous variable
    
    H_X = ee.entropyd(X)
    H_dY = ee.entropyd(dY)
    H_Z = ee.entropyd(Z)

    X_list = [ [x] for x in X ]
    dY_list = [ [dy] for dy in dY ]
    Z_list = [ [z] for z in Z ]
    
    area3 = ee.cmidd(X_list, dY_list, Z_list)
    area5 = ee.cmidd(Z_list, dY_list, X_list)
    area6 = ee.cmidd(X_list, Z_list, dY_list)
    area7 = ee.midd(X_list, dY_list) - area3
    area1 = H_dY - (area3 + area5 + area7)
    area2 = H_X - (area3 + area6 + area7)
    area4 = H_Z - (area5 + area6 + area7)

    area1 = round(area1, 2)
    area2 = round(area2, 2)
    area3 = round(area3, 2)
    area4 = round(area4, 2)
    area5 = round(area5, 2)
    area6 = round(area6, 2)
    area7 = round(area7, 2)

    print area3

    fig = plt.figure()
    venn3(subsets = (area1, area2, area3, area4, area5, area6, area7), set_labels = ('Condensability', din + " count", 'AT content'))
    plt.title("Mutual Information")
    #plt.show()
    plt.close()


    print

fig = plt.figure()

X, Y = [], []
for din in all_din:
    X.append(din)
    Y.append(din_pcorr[din])

plt.bar(range(1, len(X)+1), Y)
plt.xticks(range(1, len(X)+1), X)
plt.ylabel("partial-correlation")
plt.savefig("pcorr_din_bar.png")
plt.show()
plt.close()
"""
