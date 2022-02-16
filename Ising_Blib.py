import sys
import math
import copy
import numpy as np
from scipy.optimize import curve_fit
from sklearn import linear_model
from sklearn.utils import shuffle
import matplotlib.pyplot as plt
#import graph_edit
import random
#from SliderClass import Slider
import LinModel

def average(l):
    sum = 0.0
    count =0
    for e in l:
        if type(e) != str:
            count +=1
            sum += e
    return sum/count

def get_corr(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = np.average(x)
    avg_y = np.average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff
    return diffprod / np.sqrt(xdiff2 * ydiff2)


def rev_comp (seq):
    rev_seq=seq[::-1]; new_seq=''
    for nt in rev_seq:
        nt = nt.upper()
        if nt == 'A':
            new_seq += 'T'
        if nt == 'T':
            new_seq += 'A'
        if nt == 'C':
            new_seq += 'G'
        if nt == 'G':
            new_seq += 'C'
    return new_seq

class Ising_model (object):

    def __init__ (self,
                 seq_list,
                 score_list):
        
        self.seq_list = seq_list
        self.score_list = score_list

        self.reg = None
        self.coeff = None
        self.rsquare = None
        self.corr = None
        self.freq = None

    def _extract_features (self,
                           seq_list,
                           nts,
                           pairs,
                           rev):

        idx_nt = sorted(nts)
        idx_pair = sorted(pairs)
        
        nt_idx = {idx_nt[i]:i for i in range(len(idx_nt))}
        pair_idx = {idx_pair[i]:i for i in range(len(idx_pair))}

        # handling the rev complementary strand
        if rev:
            for nt in nts:
                rev_nt = rev_comp(nt)
                if rev_nt not in nt_idx:
                    nt_idx[rev_nt] = nt_idx[nt]
                    
            for pair in pairs:
                rev_pair = rev_comp(pair)
                if rev_pair not in pair_idx:
                    pair_idx[rev_pair] = pair_idx[pair]


        feature_list = []
        for seq in seq_list:
            hrow = [0 for k in range(len(nts))]
            jrow = [0 for k in range(len(pairs))]

            for i in range(len(seq)):
                nt1 = seq[i]
                try:
                    idx = nt_idx[nt1]
                    hrow[idx] += 1
                except:
                    pass

                if i < len(seq)-1:
                    pair = ''.join(sorted(seq[i:i+2]))
                    try:
                        idx = pair_idx[pair]
                        jrow[idx] += 1
                    except:
                        pass

            feature_list.append(hrow+jrow)
        
        return feature_list

    def train (self,
               nts,
               pairs,
               rev=True,
               alpha=0.5,
               k_fold=10,
               graph=False):

        # extract features
        print >> sys.stderr, "reading data values"
        feature_list = self._extract_features(self.seq_list, nts, pairs, rev)
        target_list = [[value] for value in self.score_list]

        # k-fold corss validation
        print >> sys.stderr, "k-fold cross validation"
        feature_list, target_list = shuffle(feature_list, target_list)
        part_size = (len(feature_list) + 1) / k_fold
        corrs = []
        for i in range(k_fold):
            st = i*part_size
            ed = min((i+1)*part_size, len(feature_list))
            train_flist = feature_list[:st] + feature_list[ed:] 
            train_tlist = target_list[:st] + target_list[ed:]
            test_flist, test_tlist = feature_list[st:ed], target_list[st:ed]
            
            print str(i+1) + '-fold'
            reg = linear_model.Ridge(alpha=alpha)
            reg.fit (train_flist, train_tlist)
            #print reg.coef_
            rsquare = reg.score(train_flist, train_tlist)
            print 'r-square: ' + str(rsquare)

            Yexp = [ value[0] for value in test_tlist]
            Ypred = reg.predict(test_flist)
            Ypred = [ value[0] for value in Ypred]
            corr = get_corr(Yexp, Ypred)
            corrs.append(corr)
            print 'correlation: ' + str(corr)

            if graph:
                fig = plt.figure()
                low = min(Yexp + Ypred)
                up = max(Yexp + Ypred)
                mid = np.median(Yexp + Ypred)
                plt.plot(Yexp, Ypred, '.')
                plt.plot([low-mid,up+mid],[low-mid,up+mid], '--')
                plt.xlim([low - mid*0.1, up + mid*0.1])
                plt.ylim([low - mid*0.1, up + mid*0.1])
                plt.show()
                plt.close()

        self.corr = corrs
        print >> sys.stderr, "Mean correlation: " + str(np.mean(self.corr))
        
        # linear fitting all data
        print >> sys.stderr, "linear regression of full data"
        reg = linear_model.Ridge(alpha=alpha)
        reg.fit (feature_list, target_list)
        self.reg = reg
        print reg.coef_
        self.rsquare = reg.score(feature_list, target_list)
        print "r-square: " + str(self.rsquare)


        # retrieve all coefficient 
        idx_parm = sorted(nts) + sorted(pairs)

        pt = 0
        nt_h = {}
        while pt < len(nts):
            nt = idx_parm[pt]
            nt_h[nt] = reg.coef_[0][pt]
            pt +=1
        assert pt == len(nts)

        pair_j = {}
        while pt < len(nts) + len(pairs):
            pair = idx_parm[pt]
            pair_j[pair] = reg.coef_[0][pt]
            pt +=1

        assert pt == len(reg.coef_[0])

        print >> sys.stderr, "Done"

        return nt_h, pair_j


class ScalarShapeModel (object):
    def __init__ (self,
                  shape_list,
                  score_list):

        self.shape_list = shape_list
        self.score_list = score_list

        self.reg = None
        self.coeff = None
        self.rsquare = None
        self.corr = None
        self.freq = None

        return


    def train (self,
               alpha=0.5,
               k_fold=10,
               graph=False):

        # extract features
        print >> sys.stderr, "reading data values"
        feature_list = self.shape_list
        target_list = [[value] for value in self.score_list]

        # k-fold corss validation
        print >> sys.stderr, "k-fold cross validation"
        feature_list, target_list = shuffle(feature_list, target_list)
        part_size = (len(feature_list) + 1) / k_fold
        corrs = []
        for i in range(k_fold):
            st = i*part_size
            ed = min((i+1)*part_size, len(feature_list))
            train_flist = feature_list[:st] + feature_list[ed:] 
            train_tlist = target_list[:st] + target_list[ed:]
            test_flist, test_tlist = feature_list[st:ed], target_list[st:ed]
            
            print str(i+1) + '-fold'
            reg = linear_model.Ridge(alpha=alpha)
            reg.fit (train_flist, train_tlist)
            #print reg.coef_
            rsquare = reg.score(train_flist, train_tlist)
            print 'r-square: ' + str(rsquare)

            Yexp = [ value[0] for value in test_tlist]
            Ypred = reg.predict(test_flist)
            Ypred = [ value[0] for value in Ypred]
            corr = get_corr(Yexp, Ypred)
            corrs.append(corr)
            print 'correlation: ' + str(corr)

            if graph:
                fig = plt.figure()
                low = min(Yexp + Ypred)
                up = max(Yexp + Ypred)
                mid = np.median(Yexp + Ypred)
                plt.plot(Yexp, Ypred, '.')
                plt.plot([low-mid,up+mid],[low-mid,up+mid], '--')
                plt.xlim([low - mid*0.1, up + mid*0.1])
                plt.ylim([low - mid*0.1, up + mid*0.1])
                plt.show()
                plt.close()

        self.corr = corrs
        print >> sys.stderr, "Mean correlation: " + str(np.mean(self.corr))
        
        # linear fitting all data
        print >> sys.stderr, "linear regression of full data"
        reg = linear_model.Ridge(alpha=alpha)
        reg.fit (feature_list, target_list)
        self.reg = reg
        print reg.coef_
        self.rsquare = reg.score(feature_list, target_list)
        print "r-square: " + str(self.rsquare)


        # retrieve all coefficient
        print >> sys.stderr, "Done"

        return


    


# read score file
def read_score (fname):
    First = True
    tname_ID_score = {} 
    for line in open(fname):
        cols = line.strip().split()
        if First:
            _, tname_list = cols[0], cols[1:]
            First = False
            continue
        
        ID, score_list = cols[0], cols[1:]
        ID = int(ID)

        for i in range(len(tname_list)):
            tname = tname_list[i]
            try:
                tname = int(tname)
            except:
                pass
            score = float(score_list[i])

            if tname not in tname_ID_score:
                tname_ID_score[tname] = {}
            tname_ID_score[tname][ID] = score
    return tname_ID_score
tname_ID_score = read_score("YWlib_sp_score.txt")
ID_score = tname_ID_score[3]

# read reference sequence
def read_ref (ref_fname):
    id_seq = {}
    for line in open(ref_fname):
        line = line.strip()
        if line.startswith('>'):
            id = int(line[4:])
            continue
        if line:
            assert id not in id_seq
            id_seq[id] = line
    return id_seq
ID_seq = read_ref('Blib.ref')

# read shape data
def read_DNAshape(fname):
    names = ['MGW', 'HelT', 'ProT', 'Roll']
    dic_list = [{} for i in range(len(names))]
    for i in range(len(names)):
        data = []
        for line in open(fname+"."+names[i]):
            line = line.strip()
            if line.startswith('>'):
                if data:
                    assert id not in dic_list[i]
                    dic_list[i][id] = data
                id = int(line[1:].strip().split('_')[1])
                data =[]
                continue
            if line:
                temp = line.split(',')
                for k in range(len(temp)):
                    try:
                        temp[k] = float(temp[k])
                    except Exception:
                        pass
                data += temp
        assert id not in dic_list[i]
        dic_list[i][id] = data
    return dic_list
ID_MGW, ID_HelT, ID_ProT, ID_Roll = read_DNAshape('php6POcc7')

# parameters
#poly_st, poly_ed = 1, 7 # homopolymer length range
clip = 0
NCPlen = 101


# get sequence and score list
seq_list, score_list = [], []
sshape_list = []
MGW_list, HelT_list, ProT_list, Roll_list = [], [], [], []
count_list = []
for ID in ID_seq:
    seq = ID_seq[ID]
    middle = len(seq)/2
    seq = seq[middle-NCPlen/2:middle+NCPlen/2+1]
    assert len(seq) == NCPlen
    seq = seq[clip:len(seq)-clip].upper()
    seq_list.append(seq)
    score_list.append(ID_score[ID])
    count_list.append(1)
    MGW, HelT, ProT, Roll = ID_MGW[ID], ID_HelT[ID], ID_ProT[ID], ID_Roll[ID]
    MGW = MGW[middle-NCPlen/2:middle+NCPlen/2+1]
    HelT = HelT[middle-NCPlen/2:middle+NCPlen/2]
    ProT = ProT[middle-NCPlen/2:middle+NCPlen/2+1]
    Roll = Roll[middle-NCPlen/2:middle+NCPlen/2]
    sshape_list.append([np.mean(MGW), np.mean(HelT), np.mean(ProT), np.mean(Roll)])
    #scalar_shape_list.append([average(MGW), average(HelT), average(ProT), average(Roll)])
    MGW_list.append(MGW)
    HelT_list.append(HelT)
    ProT_list.append(ProT)
    Roll_list.append(Roll)
shape_list = [MGW_list, HelT_list, ProT_list, Roll_list]

mname_list = []
corrs_list = []
regmodel_list = []
varnum_list = []

IsingModel = Ising_model(seq_list, score_list)


# AT content only
nt_h, _ = IsingModel.train(nts=['A'], pairs=[])
mname_list.append('A / T only')
corrs_list.append(IsingModel.corr)
regmodel_list.append(IsingModel.reg)
varnum_list.append(1)

# Ising model
nt_h, pair_h = IsingModel.train(nts=['A'], pairs=['AA', 'CC'])
mname_list.append('Ising model')
corrs_list.append(IsingModel.corr)
regmodel_list.append(IsingModel.reg)
varnum_list.append(3)

def MM_varnum(order, seqlen):
    ntsnum = (4**(order+1))
    
    if order % 2 == 0:
        palnum = 0
    else:
        palnum = 4**((order+1)/2)

    varnum = ntsnum * ((seqlen - order) / 2)

    if (seqlen - order) % 2 != 0:
        varnum += (ntsnum + palnum) / 2

    return varnum

SeqModel = LinModel.SeqLinearModel(seq_list, score_list, count_list)
# 0-th order Markov
order = 0
SeqModel.train(MM_orders=[order], Kmer_k_b=None, PolyA_b=None, GC_b=None, Harmonic=None, sym=True, graph=False)
mname_list.append('Mono-nt steps')
corrs_list.append(SeqModel.corr)
regmodel_list.append(SeqModel.reg)
varnum_list.append(MM_varnum(order, NCPlen))

# 1-th order Markov
order = 1
SeqModel.train(MM_orders=[order], Kmer_k_b=None, PolyA_b=None, GC_b=None, Harmonic=None, sym=True, graph=False)
mname_list.append('Di-nt steps')
corrs_list.append(SeqModel.corr)
regmodel_list.append(SeqModel.reg)
varnum_list.append(MM_varnum(order, NCPlen))

def Kmer_varnum (klen, seqlen):
    ntsnum = 4**klen
    
    if order % 2 == 0:
        palnum = 0
    else:
        palnum = 4**(klen/2)

    varnum = (ntsnum + palnum) / 2

    return varnum

# 4-mer model
kmer = 4
SeqModel.train(MM_orders=[], Kmer_k_b=[kmer,1], PolyA_b=None, GC_b=None, Harmonic=None, sym=True, graph=False)
mname_list.append('4-mers')
corrs_list.append(SeqModel.corr)
regmodel_list.append(SeqModel.reg)
varnum_list.append(Kmer_varnum(kmer, NCPlen))

# scalar DNA shape model
sShapeModel = ScalarShapeModel(sshape_list, score_list)
sShapeModel.train()
mname_list.append('DNA shape')
corrs_list.append(sShapeModel.corr)
regmodel_list.append(sShapeModel.reg)
varnum_list.append(4)

# full shape model
ShapeModel = LinModel.ShapeLinearModel(shape_list, score_list, count_list)
ShapeModel.train()
mname_list.append('full DNA shape')
corrs_list.append(ShapeModel.corr)
regmodel_list.append(ShapeModel.reg)
coeff_list = []
for coeff in ShapeModel.coeff.values():
    coeff_list += coeff
varnum_list.append(len(coeff_list))

color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray']

# check correlation
fig = plt.figure(figsize=(2.8, 2))
for i in range(len(varnum_list)):
    bp = plt.boxplot(corrs_list[i], positions=[i], patch_artist=True, showfliers=False, widths=[0.4])
    for patch in bp['boxes']:
        patch.set_facecolor(color_list[i])
    for element in ['whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color='black')
#plt.boxplot(corrs_list, positions=range(len(mname_list)))
plt.xticks(range(len(corrs_list)), mname_list, rotation=45, fontsize=8,  ha="right", rotation_mode="anchor")
plt.gca().tick_params('y', labelsize=6)
plt.ylabel("Pearson correlation", fontsize=8)
plt.title("10-fold cross validation", fontsize=10)
plt.xlim([-0.5, len(mname_list)-0.5])
plt.savefig('Pearson_model.svg', format='svg', bbox_inches='tight')
#plt.show()
plt.close()

# check parmeters number
fig = plt.figure(figsize=(2.8, 2))
for i in range(len(varnum_list)):
    varnum = varnum_list[i]
    plt.bar([i], [varnum+1], width=0.5, ec ='black', color=color_list[i])
    plt.gca().text(i, varnum+5, str(varnum+1), fontsize=12, va='center', color='white', alpha=0)
    #plt.annotate(str(varnum), (i, varnum+10), va='center', ha='center', fontsize=8)
#plt.bar(range(len(varnum_list)), varnum_list
plt.ylim([1, 10**(int(np.log10(max(varnum_list)))+1.5)])
plt.gca().tick_params('y', labelsize=6)
plt.xticks(range(len(varnum_list)), mname_list, rotation=45, fontsize=8, ha="right", rotation_mode="anchor")
plt.ylabel("Number of coefficients", fontsize=8)
plt.yscale("log")
plt.title("Model complexity", fontsize=10)
plt.savefig('Numcoeff_model.svg', format='svg', bbox_inches='tight')
#plt.show()
plt.close()
