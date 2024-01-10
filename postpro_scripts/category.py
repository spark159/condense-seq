import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import copy
import math
import Interval_dict

#parameters
win_size = 171
chr_target = 'chr1'
chr_name = "chromosome 1"

def norm(L):
    total = sum(L)
    return [L[i]/float(total) for i in range(len(L))]

def read_genome(fname):
    genome_size = {}
    for line in open(fname):
        line = line.strip()
        if line.startswith('>'):
            key = line[1:]
            assert key not in genome_size
            genome_size[key] = 0
            continue
        genome_size[key] += len(line)
    return genome_size
genome_size = read_genome("data/hg19.fa")

def read_NCPcov_file (fname, chr_target):
    ID_pos = {}
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            names = cols[3:]
            ID_score_list = [{} for i in range(len(names)-1)]
            #bin_GC = {}
            First = False
            continue
        ID, chr, pos = int(cols[0]), cols[1], int(cols[2])
        if chr != chr_target:
            continue
        ID_pos[ID] = pos
        counts = cols[3:]
        if float(counts[-1]) <= 0:
            counts[-1] = 1.0
        for i in range(len(ID_score_list)):
            ID_score = ID_score_list[i]
            #bin_score[BinID] = np.log(float(counts[i])/float(counts[-1]))
            if float(counts[i]) <= 0 :
                #score = np.nan
                rcount = (float(counts[-1])/float(1.0))
            else:
                rcount = (float(counts[-1])/float(counts[i]))
            ID_score[ID] = math.log(rcount)
        #bin_GC[BinID] = float(GC)
    return ID_score_list, ID_pos, names
ID_score_list, ID_pos, names = read_NCPcov_file("data/hg19_chr1_171_Ncov.cn", chr_target)
ID_score1, ID_score2 = ID_score_list

def quantile (ID_score, num, frac=None):
    def value_cmp(a, b):
        if a[1] <= b[1]:
            return -1
        else:
            return 1
    IDscore = [[ID, score] for ID, score in ID_score.items()]
    IDscore = sorted(IDscore, cmp=value_cmp)
    if frac == None:
        size_list = [int(math.ceil(len(IDscore) / float(num)))]*num
    else:
        if sum(frac) != 1:
            frac = norm(frac)
        num = len(frac)
        size_list = [int(round(f*len(IDscore))) for f in frac]
    size_list[-1] += len(IDscore) - sum(size_list)
    if size_list[-1] == 0:
        size_list[-2] -= 1
        size_list[-1] += 1
    assert sum(size_list) == len(IDscore)
    output = []
    ed = 0
    for i in range(num):
        #st = i*size
        #ed = min((i+1)*size,len(IDscore))
        size = size_list[i]
        st = ed             
        ed = st + size
        temp = [IDscore[j][0] for j in range(st,ed)]
        output.append(temp)
    return output
#frac=[(4**i) for i in range(1,11)]
#frac=[(2**i) for i in range(1,11)]
frac=[(5*i) for i in range(1,11)]
frac = norm(frac)[::-1]
group1 = quantile(ID_score1, len(frac), frac=frac)
group2 = quantile(ID_score2, len(frac), frac=frac)

fig = plt.figure()
offset = 0
for i in range(len(group1)):
    IDs = group1[i]
    X, Y = [], []
    for j in range(len(IDs)):
        if i - 1 >= 0:
            offset += len(group1[i-1])
        else:
            offset += 0
        X.append(offset + j)
        ID = IDs[j]
        score1 = ID_score1[ID]
        Y.append(score1)
    plt.fill_between(X,0,Y, label='Partition ' + str(i+1))
plt.xlabel("Nucleosome peaks")
plt.ylabel("Condensibility (A.U.)")
plt.legend()
plt.show()


def read_hgtable(fname, chr_target):
    TSS_count = 0
    TTS_count = 0
    Body_count = 0
    Prom_count = 0

    ID_interval = {}
    ID_info = {}
    
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            First = False
        _, geneID, chr, strand = cols[:4]
        #if geneID.startswith("NR"):
            #continue
        #if geneID.startswith("XR"):
            #continue
        if chr != chr_target:
            continue
        if strand == "+":
            TSS_st = max(int(cols[4])-500, 0)
            TSS_ed = min(int(cols[4])+501, genome_size[chr])
            TTS_st = max(int(cols[5])-500, 0)
            TTS_ed = min(int(cols[5])+501, genome_size[chr])
            Body_st = min(TSS_ed+1, genome_size[chr])
            Body_ed = max(TTS_st, 0)
            Prom_st = max(TSS_st-5000, 0)
            Prom_ed = max(TSS_st, 0)
        else:
            TSS_st = max(int(cols[5])-500, 0)
            TSS_ed = min(int(cols[5])+501, genome_size[chr])
            TTS_st = max(int(cols[4])-500, 0)
            TTS_ed = min(int(cols[4])+501, genome_size[chr])
            Body_st = min(TTS_ed+1, genome_size[chr])
            Body_ed = max(TSS_st, 0)
            Prom_st = min(TSS_ed+1, genome_size[chr])
            Prom_ed = min(TSS_ed+5001, 0)
        if TSS_ed - TSS_st > 0:
            ID_interval["TSS:" + str(TSS_count)] = [TSS_st, TSS_ed]
            TSS_count +=1
        if TTS_ed - TTS_st > 0:
            ID_interval["TTS:" + str(TTS_count)] = [TTS_st, TTS_ed]
            TTS_count +=1
        if Body_ed - Body_st > 0:
            ID_interval["Body:" + str(Body_count)] = [Body_st, Body_ed]
            Body_count +=1
        if Prom_ed - Prom_st > 0:
            ID_interval["Prom:" + str(Prom_count)] = [Prom_st, Prom_ed]
            Prom_count +=1
    return ID_interval
ID_interval = read_hgtable("data/hgTables", chr_target)
Region_category = Interval_dict.double_hash(ID_interval, 10000, genome_size[chr_target])

def categorize (IDs, ID_pos, Region_category):
    category_IDs = {"inter":[], "multi":[], "TSS":[], "TTS":[], "Prom":[], "Body":[]}
    category_count = {}
    for ID in IDs:
        pos = ID_pos[ID]
        find_Regions = Region_category.find(pos)
        if len(find_Regions) == 0:
            category_IDs["inter"].append(ID)
            continue
        temp = set([])
        for Region in find_Regions:
            typename = Region.split(':')[0]
            temp |= set([typename])
        find_Types = list(temp)
        if len(find_Types) > 1 :
            category_IDs["multi"].append(ID)
            continue
        typename = find_Types[0]
        if typename.startswith("TSS"):
            category_IDs["TSS"].append(ID)
        elif typename.startswith("TTS"):
            category_IDs["TTS"].append(ID)
        elif typename.startswith("Body"):
            category_IDs["Body"].append(ID)
        elif typename.startswith("Prom"):
            category_IDs["Prom"].append(ID)
    for category in category_IDs:
        category_count[category] = len(category_IDs[category])
    return category_IDs, category_count

category_IDs_list1 = []
category_count_list1 = []
for IDs in group1:
    category_IDs, category_count = categorize(IDs, ID_pos, Region_category)
    category_IDs_list1.append(category_IDs)
    category_count_list1.append(category_count)

category_IDs_list2 = []
category_count_list2 = []
for IDs in group2:
    category_IDs, category_count = categorize(IDs, ID_pos, Region_category)
    category_IDs_list2.append(category_IDs)
    category_count_list2.append(category_count)

    
cat_names = ["Prom", "TSS", "Body", "TTS", "inter", "multi"]
X = []
Y_list = [ [] for i in range(len(cat_names)) ]
for i in range(len(category_count_list1)):
    X.append(i+1)
    category_count = category_count_list1[i]
    for j in range(len(cat_names)):
        category = cat_names[j]
        count = category_count[category]
        rcount = float(count) / sum(category_count.values())
        Y_list[j].append(rcount)
fig = plt.figure()
for i in range(len(Y_list)):
    Y = Y_list[i]
    plt.plot(X, Y , '--o', label=cat_names[i])
plt.xticks(range(1, len(frac)+1))
plt.xlabel("Partition by condensibility")
plt.ylabel("Fraction")
plt.legend()
plt.show()
    
X = []
Y_list = [ [] for i in range(len(cat_names)) ]
for i in range(len(category_count_list2)):
    X.append(i+1)
    category_count = category_count_list2[i]
    for j in range(len(cat_names)):
        category = cat_names[j]
        count = category_count[category]
        rcount = float(count) / sum(category_count.values())
        Y_list[j].append(rcount)
fig = plt.figure()
for i in range(len(Y_list)):
    Y = Y_list[i]
    plt.plot(X, Y, '--o',  label=cat_names[i])
plt.xticks(range(1, len(frac)+1))
plt.xlabel("Partition by condensibility")
plt.ylabel("Fraction")
plt.legend()
plt.show()
