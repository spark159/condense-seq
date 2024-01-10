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

def norm(L):
    total = 0.0
    for value in L:
        total += value
    return [value/total for value in L]

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

def is_pal (seq):
    if len(seq) % 2 != 0:
        return False
    for i in range(len(seq)/2):
        if nt[i] != rev_comp(nt[len(seq)-1-i]):
            return False
    return True

def all_path(N, states, arrange=True):
    temp = []
    if N==1:
        temp = list(states)
    else:
        for path in all_path(N-1, states):
            for state in states:
                temp.append(path+state)
    if not arrange:
        return temp
    unique, non_pals = [], []
    for seq in temp:
        if rev_comp(seq) not in unique:
            unique.append(seq)
        else:
            non_pals.append(seq)
    output = unique + non_pals
    assert len(temp) == len(output)
    return output

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

def Amer_len(seq, pos=True):
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

def GC_content(seq):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    return (num/float(len(seq)))*100

class ShapeLinearModel:
    def __init__ (self, shape_list, score_list, count_list):
        for i in range(len(shape_list)):
            assert len(shape_list[i]) == len(shape_list[0])
        self.MGW_list, self.HelT_list, self.ProT_list, self.Roll_list = shape_list
        self.WINlen = len(self.MGW_list[0])
        for i in range(len(self.MGW_list)):
            assert len(self.MGW_list[i]) == len(self.ProT_list[i]) == self.WINlen
            assert len(self.HelT_list[i]) == len(self.Roll_list[i]) == self.WINlen-1
        assert len(score_list) == len(count_list) == len(shape_list[0])
        self.score_list = score_list
        self.count_list = count_list
        self.reg = None
        self.coeff = None
        self.rsquare = None
        self.corr = None
        self.freq = None
    
    def _even_sampling(self, sym=True):
        MGW_samples, HelT_samples, ProT_samples, Roll_samples = [], [], [], []
        for i in range(len(self.MGW_list)):
            MGW_samples.append(self.MGW_list[i])
            HelT_samples.append(self.HelT_list[i])
            ProT_samples.append(self.ProT_list[i])
            Roll_samples.append(self.Roll_list[i])
            if sym:
                MGW_samples.append(self.MGW_list[i][::-1])
                HelT_samples.append(self.HelT_list[i][::-1])
                ProT_samples.append(self.ProT_list[i][::-1])
                Roll_samples.append(self.Roll_list[i][::-1])
        return [MGW_samples, HelT_samples, ProT_samples, Roll_samples]


    def _bias_sampling(self, scale=100, sym=True):
        MGW_samples, HelT_samples, ProT_samples, Roll_samples = [], [], [], []
        for i in range(len(self.MGW_list)):
            raw_count = self.count_list[i]
            count = int(round(raw_count*scale))
            for i in range(count):
                MGW_samples.append(self.MGW_list[i])
                HelT_samples.append(self.HelT_list[i])
                ProT_samples.append(self.ProT_list[i])
                Roll_samples.append(self.Roll_list[i])
                if sym:
                    MGW_samples.append(self.MGW_list[i][::-1])
                    HelT_samples.append(self.HelT_list[i][::-1])
                    ProT_samples.append(self.ProT_list[i][::-1])
                    Roll_samples.append(self.Roll_list[i][::-1])
        return [MGW_samples, HelT_samples, ProT_samples, Roll_samples]

    
    def _stat_shape(self, shape_list_list):
        freqs_list = []
        for shape_list in shape_list_list:
            freqs = np.sum(shape_list, axis=0) / float(len(shape_list))
            freqs_list.append(freqs)
        return freqs_list

    def report (self,
                scale=100,
                sym=True,
                norm=True,
                graph=False):
        
        freq = {}
        print >> sys.stderr, "sampling data values"        
        bias_samples_list = self._bias_sampling(scale=scale, sym=sym)
        bias_freqs_list = self._stat_shape(bias_samples_list)
        if norm:
            even_samples_list = self._even_sampling(sym=sym)
            even_freqs_list = self._stat_shape(even_samples_list)
        names = ["MGW", "HelT", "ProT", "Roll"]
        for i in range(len(even_samples_list)):
            if norm:
                freq[names[i]] = bias_freqs_list[i] - even_freqs_list[i]
            else:
                freq[names[i]] = bias_freqs_list[i]
            self.freq = freq
        print >> sys.stderr, "Done"
        return self.freq

    def spectrum(self,
                 gnum=5,
                 reverse=False,
                 sym=True,
                 norm=True):

        score_idx = [ (self.score_list[i], i) for i in range(len(self.score_list)) ]
        score_idx = sorted(score_idx, reverse=reverse)

        MGW_list, HelT_list, ProT_list, Roll_list = [], [], [], []
        for score, idx in score_idx:
            MGW_list.append(self.MGW_list[idx])
            HelT_list.append(self.HelT_list[idx])
            ProT_list.append(self.ProT_list[idx])
            Roll_list.append(self.Roll_list[idx])

        group_shape_list_list = []
        
        size = (len(self.score_list)) / gnum
        for i in range(gnum):
            st = i*size
            ed = st + size
            if i == gnum-1:
                ed = len(self.score_list)
            shape_list_list = []
            shape_list_list.append(MGW_list[st:ed])
            shape_list_list.append(HelT_list[st:ed])
            shape_list_list.append(ProT_list[st:ed])
            shape_list_list.append(Roll_list[st:ed])
            group_shape_list_list.append(shape_list_list)

        names = ["MGW", "HelT", "ProT", "Roll"]
        group_freq = [{} for i in range(gnum)]
        for i in range(gnum):
            shape_list_list = group_shape_list_list[i]
            freqs_list = self._stat_shape(shape_list_list)
            for j in range(len(names)):
                name = names[j]
                freqs = freqs_list[j]
                group_freq[i][name] = freqs

        if norm:
            norm_freqs_list = self._stat_shape([MGW_list, HelT_list, ProT_list, Roll_list])
            for i in range(gnum):
                for j in range(len(names)):
                    name = names[j]
                    group_freq[i][name] = group_freq[i][name] - norm_freqs_list[j]
            
        return group_freq
            

    def _var_shape(self, shape_list_list, sym=True):
        var_list = []
        MGW_list, HelT_list, ProT_list, Roll_list = shape_list_list
        for i in range(len(MGW_list)):
            row = []
            row += MGW_list[i]
            row += HelT_list[i]
            row += ProT_list[i]
            row += Roll_list[i]
            var_list.append(row)
        return var_list

    def train (self,
               alpha=0.5,
               k_fold=10,
               sym=True,
               graph=False):
        
        # reading variables
        print >> sys.stderr, "reading data values"
        shape_list_list = [self.MGW_list, self.HelT_list, self.ProT_list, self.Roll_list]
        var_list = self._var_shape(shape_list_list, sym=sym)

        feature_list = var_list
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
        pt = 0
        self.coeff = {}
        names = ["MGW", "HelT", "ProT", "Roll"]
        for i in range(len(names)):
            name = names[i]
            temp = []
            for j in range(len(shape_list_list[i][0])):
                temp.append(reg.coef_[0][pt])
                pt +=1
            self.coeff[name] = temp
            
        assert len(reg.coef_[0]) == pt
        print >> sys.stderr, "Done"
        return self.coeff


    def _predict(self, target_shape_profile, bound, sym):
        MGW_profile, HelT_profile, ProT_profile, Roll_profile = target_shape_profile
        assert len(MGW_profile) == len(HelT_profile) + 1 == len(ProT_profile) == len(Roll_profile) + 1
        MGW_list, HelT_list, ProT_list, Roll_list = [], [], [], []
        pos_list = []
        for i in range(bound, len(MGW_profile) - bound - self.WINlen):
            pos_list.append(bound+i+self.WINlen/2)
            MGW_list.append(MGW_profile[i:i+self.WINlen])
            HelT_list.append(HelT_profile[i:i+self.WINlen-1])
            ProT_list.append(ProT_profile[i:i+self.WINlen])
            Roll_list.append(Roll_profile[i:i+self.WINlen-1])
        shape_list = [MGW_list, HelT_list, ProT_list, Roll_list]
        return shape_list, pos_list        
    
    def predict(self, target_shape_list, bound, sym):
        shape_list, pos_list = _predict(target_shape_list, bound=bound, sym=sym)
        var_list = self._var_shape(shape_list)
        Ypred = self.reg.predict(var_list)
        Ypred = [ value[0] for value in Ypred]
        pred_score = Ypred
        return pred_score, pos_list
    

    def display(self,
                data_dict,
                vmin=None,
                vmax=None):

        names = ["MGW", "HelT", "ProT", "Roll"]
        fig = plt.figure()
        for name in names:
            freq = data_dict[name]
            plt.plot(freq, label=name)
        plt.legend()
        plt.show()
        plt.close()
        

class SeqLinearModel:
    def __init__(self, seq_list, score_list, count_list, pos_list=None, shape_list=None):
        assert len(seq_list) == len(score_list) == len(count_list)
        self.NCPlen = len(seq_list[0])
        for seq in seq_list:
            assert len(seq) == self.NCPlen
        self.seq_list = seq_list
        self.score_list = score_list
        self.count_list = count_list
        self.pos_list = pos_list
        if shape_list:
            self.shape = True
            self.MGW_list, self.HelT_list, self.ProT_list, self.Roll_list = shape_list
        else:
            self.shape = False
            self.MGW_list, self.HelT_list, self.ProT_list, self.Roll_list = [], [], [], []
        self.reg = None
        self.coeff = None
        self.rsquare = None
        self.corr = None
        self.freq = None

    def _even_sampling(self, sym=True, shape=False):
        if not shape:
            seq_samples = []
            for i in range(len(self.seq_list)):
                seq = self.seq_list[i]
                seq_samples.append(seq)
                if sym:
                    seq_samples.append(rev_comp(seq))
            return seq_samples
        if shape:
            MGW_samples, HelT_samples, ProT_samples, Roll_samples = [], [], [], []
            for i in range(len(self.seq_list)):
                MGW_samples.append(self.MGW_list[i])
                HelT_samples.append(self.HelT_list[i])
                ProT_samples.append(self.ProT_list[i])
                Roll_samples.append(self.Roll_list[i])
                if sym:
                    MGW_samples.append(self.MGW_list[i][::-1])
                    HelT_samples.append(self.HelT_list[i][::-1])
                    ProT_samples.append(self.ProT_list[i][::-1])
                    Roll_samples.append(self.Roll_list[i][::-1])
            return [MGW_samples, HelT_samples, ProT_samples, Roll_samples]


    def _bias_sampling(self, scale=100, sym=True, shape=False):
        if not shape:
            seq_samples = []
            for i in range(len(self.seq_list)):
                seq = self.seq_list[i]
                raw_count = self.count_list[i]
                count = int(round(raw_count*scale))
                for i in range(count):
                    seq_samples.append(seq)
                    if sym:
                        seq_samples.append(rev_comp(seq))
            return seq_samples
        if shape:
            MGW_samples, HelT_samples, ProT_samples, Roll_samples = [], [], [], []
            for i in range(len(self.seq_list)):
                raw_count = self.count_list[i]
                count = int(round(raw_count*scale))
                for i in range(count):
                    MGW_samples.append(self.MGW_list[i])
                    HelT_samples.append(self.HelT_list[i])
                    ProT_samples.append(self.ProT_list[i])
                    Roll_samples.append(self.Roll_list[i])
                    if sym:
                        MGW_samples.append(self.MGW_list[i][::-1])
                        HelT_samples.append(self.HelT_list[i][::-1])
                        ProT_samples.append(self.ProT_list[i][::-1])
                        Roll_samples.append(self.Roll_list[i][::-1])

            return [MGW_samples, HelT_samples, ProT_samples, Roll_samples]

        
    def _stat_Markov(self, seq_list, order):
        ntdic = {}
        for nt in all_path(order+1, 'ATCG'):
            ntdic[nt] = 0.0

        sample_num = len(seq_list)
        
        freq = [ copy.deepcopy(ntdic) for i in range(self.NCPlen - order) ]
        for seq in seq_list:
            for i in range(len(seq) - order):
                nt = seq[i:i+1+order]
                freq[i][nt] += 1.0 / sample_num
        
        mean, std = [], []
        for ntdic in freq:
            mean.append(np.mean(ntdic.values()))
            std.append(np.std(ntdic.values()))

        stdz_freq = []
        for i in range(len(freq)):
            ntdic = freq[i]
            temp = {}
            for nt, value in ntdic.items():
                temp[nt] = (value - mean[i]) / std[i]
            stdz_freq.append(temp)

        return freq, sample_num, mean, std, stdz_freq
        
    def _stat_Kmer(self, seq_list, knum, bnum):
        seqlen = self.NCPlen / bnum
        assert seqlen >= knum
        
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2

        ntdic = {}
        for nt in all_path(knum, 'ATCG'):
            ntdic[nt] = 0.0
        freq = [ copy.deepcopy(ntdic) for i in range(bnum) ]

        sample_num = len(seq_list)*(seqlen-knum+1)
        
        for seq in seq_list:
            for k in range(bnum):
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq = seq[st:st+seqlen]                
                for i in range(seqlen - knum + 1):
                    nt = bseq[i:i+knum]
                    freq[k][nt] += 1.0 / sample_num

        mean, std = [], []
        for ntdic in freq:
            mean.append(np.mean(ntdic.values()))
            std.append(np.std(ntdic.values()))

        stdz_freq = []
        for i in range(len(freq)):
            ntdic = freq[i]
            temp = {}
            for nt, value in ntdic.items():
                temp[nt] = (value - mean[i]) / std[i]
            stdz_freq.append(temp)

        return freq, sample_num, mean, std, stdz_freq

    def _stat_PolyA(self, seq_list, bnum):
        seqlen = self.NCPlen / bnum
        assert seqlen >= 1
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2

        freq = [ [] for i in range(bnum) ]
        sample_num = len(seq_list)
        
        for seq in seq_list:
            for k in range(bnum):
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq = seq[st:st+seqlen]
                num_pos = Amer_len(bseq)
                score = 0.0
                for num, pos in num_pos.items():
                    if num >= 5:
                        score += num*len(pos)
                freq[k].append(score)
                
        mean, std = [], []
        for k in range(len(freq)):
            mean.append(np.mean(freq[k]))
            std.append(np.std(freq[k]))

        return freq, sample_num, mean, std

    def _stat_GC(self, seq_list, bnum):
        seqlen = self.NCPlen / bnum
        assert seqlen >= 1
        
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2

        freq = [ [] for i in range(bnum) ]

        sample_num = len(seq_list)
        
        for seq in seq_list:
            for k in range(bnum):
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq = seq[st:st+seqlen]
                GC = GC_content(bseq)
                freq[k].append(GC)

        mean, std = [], []
        for k in range(len(freq)):
            mean.append(np.mean(freq[k]))
            std.append(np.std(freq[k]))
        return freq, sample_num, mean, std

    def _stat_shape(self, shape_list):
        shape_list = np.asarray(shape_list)
        #print shape_list.shape
        freq = np.sum(shape_list, axis=0) / len(shape_list)
        return freq

    def report (self,
                MM_orders,
                Kmer_k_b,
                PolyA_b,
                GC_b,
                Harmonic,
                scale=100,
                shape=False,
                sym=True,
                norm=True,
                graph=False):
        
        freq = {}
        if shape:
            print >> sys.stderr, "sampling data values"        
            even_samples_list = self._even_sampling(sym=sym, shape=True)
            bias_samples_list = self._bias_sampling(scale=scale, sym=sym, shape=True)
            names = ["MGW", "HelT", "ProT", "Roll"] 
            for i in range(len(even_samples_list)):
                even_samples, bias_samples = even_samples_list[i], bias_samples_list[i]
                freq[names[i]] = self._stat_shape(bias_samples) / self._stat_shape(even_samples)
                self.freq = freq
        else:
            # sampling data
            print >> sys.stderr, "sampling data values"        
            even_samples = self._even_sampling(sym=sym)
            bias_samples = self._bias_sampling(scale=scale, sym=sym)
        
        # frequency counts
        print >> sys.stderr, "counting sequence features"
        if MM_orders:
            for order in sorted(MM_orders):
                name = 'MM' + str(order)
                freq1, sample_num, mean, std, stdz_freq = self._stat_Markov(even_samples, order)
                freq2, sample_num, mean, std, stdz_freq = self._stat_Markov(bias_samples, order)
                freq_fold = [{} for i in range(self.NCPlen-order)]
                nts = all_path(order+1, 'ATCG')
                for i in range(self.NCPlen-order):
                    for nt in nts:
                        freq_fold[i][nt] = freq2[i][nt] / freq1[i][nt]
                freq[name] = freq_fold        
        if Kmer_k_b:
            knum, bnum = Kmer_k_b
            freq1, sample_num, mean, std, stdz_freq = self._stat_Kmer(even_samples, knum, bnum)
            freq2, sample_num, mean, std, stdz_freq = self._stat_Kmer(bias_samples, knum, bnum)
            nts = all_path(knum, 'ATCG')
            for i in range(bnum):
                name = 'Kmer' + str(i)
                freq_fold = {}
                for nt in nts:
                    freq_fold[nt] = freq2[i][nt] / freq1[i][nt]
                freq[name] = freq_fold
        if PolyA_b:
            bnum = PolyA_b
            freq1, sample_num, mean1, std = self._stat_PolyA(even_samples, bnum)
            freq2, sample_num, mean2, std = self._stat_PolyA(bias_samples, bnum)
            for i in range(bnum):
                name = 'PolyA' + str(i)
                freq[name] = mean2[i]/mean1[i]
        if GC_b:
            bnum = GC_b
            freq1, sample_num, mean1, std = self._stat_GC(even_samples, bnum)
            freq2, sample_num, mean2, std = self._stat_GC(bias_samples, bnum)
            for i in range(bnum):
                name = 'GC' + str(i)
                freq[name] = mean2[i]/mean1[i]
        if Harmonic:
            None

            # To do
        self.freq = freq
        print >> sys.stderr, "Done"
        return None

    def spectrum (self,
                  MM_orders,
                  Kmer_k_b,
                  PolyA_b,
                  GC_b,
                  Harmonic,
                  gnum=5,
                  reverse=False,
                  sym=True,
                  norm=True):
        
        score_idx = [ (self.score_list[i], i) for i in range(len(self.score_list)) ]
        score_idx = sorted(score_idx, reverse=reverse)

        total_seq_list = []
        for score, idx in score_idx:
            total_seq_list.append(self.seq_list[idx])
            if sym:
                total_seq_list.append(rev_comp(self.seq_list[idx]))
         
        group_seq_list = []
        
        size = (len(total_seq_list)) / gnum
        for i in range(gnum):
            st = i*size
            ed = st + size
            if i == gnum-1:
                ed = len(total_seq_list)
            group_seq_list.append(total_seq_list[st:ed])

        if sym:
            for i in range(gnum):
                new_seq_list = []
                for seq in group_seq_list[i]:
                    new_seq_list.append(rev_comp(seq))
                group_seq_list[i] += new_seq_list
                
        group_freq = []
        for i in range(gnum):
            freq = {}
            seq_list = group_seq_list[i]
            if MM_orders:
                for order in MM_orders:
                    name = 'MM' + str(order)
                    freq1, sample_num, mean, std, stdz_freq = self._stat_Markov(seq_list, order) 
                    if norm:
                        freq2, sample_num, mean, std, stdz_freq = self._stat_Markov(total_seq_list, order) 
                        freq_fold = [{} for i in range(len(freq1))]
                        for i in range(len(freq1)):
                            for nt in freq1[i]:
                                freq_fold[i][nt] = freq1[i][nt] / freq2[i][nt]
                        freq[name] = freq_fold        
                    else:
                        freq[name] = freq1
                        
            if Kmer_k_b:
                knum, bnum = Kmer_k_b
                freq1, sample_num, mean, std, stdz_freq = self._stat_Kmer(seq_list, knum, bnum)
                if norm:
                    freq2, sample_num, mean, std, stdz_freq = self._stat_Kmer(total_seq_list, knum, bnum)
                for i in range(bnum):
                    name = 'Kmer' + str(i)
                    if norm:
                        freq_fold = {}
                        for nt in freq1[i]:
                            freq_fold[nt] = freq1[i][nt] / freq2[i][nt]
                        freq[name] = freq_fold
                    else:
                        freq[name] = freq1[i]

            if PolyA_b:
                bnum = PolyA_b
                freq1, sample_num, mean1, std = self._stat_PolyA(seq_list, bnum)
                if norm:
                    freq2, sample_num, mean2, std = self._stat_PolyA(total_seq_list, bnum)
                for i in range(bnum):
                    name = 'PolyA' + str(i)
                    if norm:
                        freq[name] = mean1[i]/mean2[i]
                    else:
                        freq[name] = mean1[i]

            if GC_b:
                bnum = GC_b
                freq1, sample_num, mean1, std = self._stat_GC(seq_list, bnum)
                if norm:
                    freq2, sample_num, mean2, std = self._stat_GC(total_seq_list, bnum)
                for i in range(bnum):
                    name = 'GC' + str(i)
                    if norm:
                        freq[name] = mean1[i]/mean2[i]
                    else:
                        freq[name] = mean1[i]

            if Harmonic:
                None
            group_freq.append(freq)
            
        return group_freq
    

    def _var_Markov (self, seq_list, order, sym):
        nt_pos = {}
        nts = all_path(order+1, 'ATCG')
        for i in range(len(nts)):
            nt = nts[i]
            nt_pos[nt] = i
        if order % 2 == 0:
            palnum = 0
        else:
            palnum = 4**((order+1)/2)
        var_list = []
        for seq in seq_list:
            left = seq
            right = rev_comp(seq)
            row = []
            if sym:
                for i in range((len(seq) - order)/2):
                    nt1 = left[i:i+order+1]
                    nt2 = right[i:i+order+1]
                    temp = [0.0] * len(nts)
                    temp[nt_pos[nt1]] += 1
                    temp[nt_pos[nt2]] += 1
                    row += temp
                if (len(seq) - order) % 2 != 0:
                    #assert order % 2 == 0
                    i = (len(seq) - order)/2
                    #i = len(seq)/2 - (order+1)/2
                    nt = left[i:i+order+1]
                    temp = [0.0] * ((len(nts)+palnum)/2)
                    pos = min(nt_pos[nt], nt_pos[rev_comp(nt)])
                    temp[pos] += 1
                    row += temp
            else:
                for i in range(len(seq) - order):
                    nt = left[i:i+order+1]
                    temp = [0.0]* len(nts)
                    temp[nt_pos[nt]] += 1
                    row += temp
            var_list.append(row)
        return var_list

    def _var_Kmer (self, seq_list, knum, bnum, sym):        
        seqlen = self.NCPlen / bnum
        assert seqlen >= knum        
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2
        nt_pos = {}
        nts = all_path(knum, 'ATCG')
        for i in range(len(nts)):
            nt = nts[i]
            nt_pos[nt] = i
        if knum % 2 != 0:
            palnum = 0
        else:
            palnum = 4**(knum/2)
        var_list = []
        for seq in seq_list:
            left = seq
            right = rev_comp(seq)
            row = []
            if sym:
                for k in range(bnum/2):
                    temp = [0]*len(nts)
                    if k < bnum/2:
                        st = boundoff + k*seqlen
                    if k >= bnum/2:
                        st = boundoff + centeroff + k*seqlen
                    bseq1 = left[st:st+seqlen]
                    bseq2 = right[st:st+seqlen]
                    for i in range(seqlen-knum+1):
                        nt1 = bseq1[i:i+knum]
                        nt2 = bseq2[i:i+knum]
                        temp[nt_pos[nt1]] += 1
                        temp[nt_pos[nt2]] += 1
                    row += temp
                if bnum % 2 != 0:
                    #assert seqlen % 2 != 0
                    #st = len(seq)/2 - seqlen/2
                    st = boundoff + centeroff + (bnum/2)*seqlen
                    bseq = left[st:st+seqlen]
                    temp = [0] * ((len(nts)+palnum)/2)
                    for i in range(seqlen-knum+1):
                        nt = bseq[i:i+knum]
                        pos = min(nt_pos[nt], nt_pos[rev_comp(nt)])
                        temp[pos] += 1
                    row += temp
            else:
                for k in range(bnum):
                    temp = [0]*len(nts)
                    if k < bnum/2:
                        st = boundoff + k*seqlen
                    if k >= bnum/2:
                        st = boundoff + centeroff + k*seqlen
                    bseq = left[st:st+seqlen]
                    for i in range(seqlen-knum+1):
                        nt = bseq[i:i+knum]
                        temp[nt_pos[nt]] += 1
                    row += temp
            var_list.append(row)
        return var_list

    def _var_PolyA (self, seq_list, bnum, sym):        
        seqlen = self.NCPlen / bnum
        assert seqlen >= 1
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2        
        var_list = []
        for seq in seq_list:
            row = []
            for k in range(bnum):
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq = seq[st:st+seqlen]
                num_pos = Amer_len(bseq)
                count = 0
                for num, pos in num_pos.items():
                    if num >=5:
                        count += len(pos)
                row.append(count)
            if sym:
                sym_row = [row[i] + row[::-1][i] for i in range(bnum/2)]
                if bnum % 2 != 0:
                    sym_row += [row[bnum/2]]
                row = sym_row
            var_list.append(row)
        return var_list

    def _var_GC (self, seq_list, bnum, sym):
        seqlen = self.NCPlen / bnum
        assert seqlen >= 1
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2        
        var_list = []
        for seq in seq_list:
            row = []
            for k in range((bnum+1)/2):
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq = seq[st:st+seqlen]
                GC = GC_content(bseq)
                row.append(GC)
            if sym:
                sym_row = [row[i] + row[::-1][i] for i in range(bnum/2)]
                if bnum % 2 != 0:
                    sym_row += [row[bnum/2]]
                row = sym_row
            var_list.append(sym_row)
        return var_list

    def _var_Harmonic (self, pos_list):
        offset = self.templatelen / 2
        var_list = []
        for pos in pos_list:
            row = []
            row.append(0.5*((pos-offset)**2))
            var_list.append(row)
        return var_list

    def _var_dPolyA (self, seq_list, sym, lmin=3, lmax=20):
        var_list = []
        for seq in seq_list:
            row = []
            num_pos = Amer_len(seq)
            for i in range(lmin, lmax+1):
                temp = [0.0]*len(seq)
                try:
                    pos_list = num_pos[i]
                except:
                    if sym:
                        row += temp[:(self.NCPlen + 1)/2]
                    else:
                        row += temp
                    continue
                for pos in pos_list:
                    for j in range(i):
                        temp[pos+j] += 1
                if sym:
                    sym_temp = [temp[i] + temp[::-1][i] for i in range(self.NCPlen/2)]
                    if self.NCPlen % 2 != 0:
                        sym_temp += [temp[self.NCPlen/2]]
                    temp = sym_temp
                row += temp
            var_list.append(row)
        return var_list
    
    def train (self,
               MM_orders,
               Kmer_k_b,
               PolyA_b,
               GC_b,
               Harmonic,
               ref_key=None,
               dPolyA=False,
               alpha=0.5,
               k_fold=10,
               sym=True,
               graph=False):

        self.MM_orders, self.Kmer_k_b, = MM_orders, Kmer_k_b
        self.PolyA_b, self.GC_b, self.Harmonic = PolyA_b, GC_b, Harmonic
        score_list, seq_list, pos_list = self.score_list, self.seq_list, self.pos_list
        
        # reading variables
        print >> sys.stderr, "reading data values"
        if MM_orders:
            MM_vars = [ [] for i in range(len(seq_list)) ] 
            for order in sorted(MM_orders):
                temp = self._var_Markov(seq_list, order, sym=sym)
                for i in range(len(temp)):
                    MM_vars[i] += temp[i]
        if Kmer_k_b:
            knum, bnum = Kmer_k_b
            Kmer_vars = self._var_Kmer(seq_list, knum, bnum, sym=sym)
        if PolyA_b:
            PolyA_vars = self._var_PolyA(seq_list, PolyA_b, sym=sym)
        if GC_b:
            GC_vars = self._var_GC(seq_list, GC_b, sym=sym)
        if Harmonic:
            Harmonic_vars = self._var_Harmonic(pos_list)
        if dPolyA:
            lmin, lmax = dPolyA
            dPolyA_vars = self._var_dPolyA(seq_list, lmin=lmin, lmax=lmax, sym=sym)

        var_list = [ [] for i in range(len(seq_list))]
        for i in range(len(seq_list)):
            if MM_orders:
                var_list[i] += MM_vars[i]
            if Kmer_k_b:
                var_list[i] += Kmer_vars[i]
            if PolyA_b:
                var_list[i] += PolyA_vars[i]
            if GC_b:
                var_list[i] += GC_vars[i]
            if Harmonic:
                var_list[i] += Harmonic_vars[i]
            if dPolyA:
                var_list[i] += dPolyA_vars[i]

        """
        # adjust variables by reference key
        if ref_key:
            new_var_list = []
            new_nlogprob_list = []
            knum = self.key_list.index(ref_key)
            gnum = self.templatelen - self.NCPlen + 1 - 2*self.bound
            st, ed = knum*gnum, (knum+1)*gnum
            ref_var = np.asarray(var_list[st:ed])
            ref_nlogprob = np.asarray(nlogprob_list[st:ed])
            ref_prob = [np.exp(-value) for value in ref_nlogprob]
            #fig = plt.figure()
            #plt.plot(ref_prob)
            
            for i in range(len(self.key_list)):
                if i == knum:
                    continue
                var_group = np.asarray(var_list[i*gnum:(i+1)*gnum])
                nlogprob_group = np.asarray(nlogprob_list[i*gnum:(i+1)*gnum])
                new_var_list += list(var_group - ref_var)
                new_nlogprob_list += list((nlogprob_group - ref_nlogprob)/ref_nlogprob)
                #print self.key_list[i]
                #plt.title(self.key_list[i])
                #prob = [np.exp(-value) for value in nlogprob_group]
                #plt.plot(ref_prob)
                #plt.plot(prob)
                #plt.show()
            var_list = new_var_list
            nlogprob_list = new_nlogprob_list

        
        # adjust variables by reference point
        gnum = self.templatelen - self.NCPlen + 1 - 2*self.bound
        if gnum > 1:
            feature_list, target_list = [], []
            i = 0
            count = 0
            while i < len(var_list):
                ref_var = copy.deepcopy(var_list[i])
                ref_nlogprob = copy.deepcopy(nlogprob_list[i])
                j = 1
                count += 1
                while j < gnum:
                    row = np.asarray(var_list[i+j]) - np.asarray(ref_var)
                    feature_list.append(row)
                    target_list.append([nlogprob_list[i+j] - ref_nlogprob])
                    j +=1
                    count += 1
                i += j
            assert count == len(var_list)
        else:
            feature_list = var_list
            target_list = [[value] for value in nlogprob_list]
        """

        feature_list = var_list
        target_list = [[value] for value in score_list]

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
        pt = 0
        self.coeff = {}
        if MM_orders:
            for order in MM_orders:
                name = 'MM' + str(order)
                self.coeff[name] = [{} for i in range(self.NCPlen-order)]
                nts = all_path(order+1, 'ATCG')
                if order % 2 == 0:
                    palnum = 0
                else:
                    palnum = 4**((order+1)/2)
                if sym:
                    repeat = (self.NCPlen - order + 1)/2
                else:
                    repeat = self.NCPlen - order
                for k in range(repeat):
                    count = 0
                    #if order % 2 == 0 and k == (self.NCPlen - order + 1)/2 - 1:
                    if sym and (self.NCPlen - order) % 2 != 0 and k == (self.NCPlen - order + 1)/2 - 1:
                        total = (len(nts)+palnum)/2
                    else:
                        total = len(nts)
                    while count < total:
                        nt = nts[count]
                        self.coeff[name][k][nt] = reg.coef_[0][pt]
                        if sym:
                            self.coeff[name][self.NCPlen-order-1-k][rev_comp(nt)] = reg.coef_[0][pt]
                        count += 1
                        pt += 1
        if Kmer_k_b:
            knum, bnum = Kmer_k_b
            if knum % 2 != 0:
                palnum = 0
            else:
                palnum = 4**(knum/2)
            nts = all_path(knum, 'ATCG')
            if sym:
                repeat = (bnum+1)/2
            else:
                repeat = bnum
            for i in range(repeat):
                name1 = 'Kmer' + str(i)
                self.coeff[name1] = {}
                if sym:
                    name2 = 'Kmer' + str(bnum-1-i)
                    self.coeff[name2] = {}
                count = 0
                if sym and bnum % 2 != 0 and i == (bnum+1)/2 - 1:
                    total = (len(nts)+palnum)/2
                else:
                    total = len(nts)
                while count < total:
                    nt = nts[count]
                    self.coeff[name1][nt] = reg.coef_[0][pt]
                    if sym:
                        self.coeff[name2][rev_comp(nt)] = reg.coef_[0][pt]
                    count += 1
                    pt += 1
        if PolyA_b:
            bnum = PolyA_b
            if sym:
                repeat = (bnum+1)/2
            else:
                repeat = bnum
            for i in range(repeat):
                name1 = 'PolyA' + str(i)
                self.coeff[name1] = reg.coef_[0][pt]
                if sym:
                    name2 = 'PolyA' + str(bnum-1-i)
                    self.coeff[name2] = reg.coef_[0][pt]
                pt += 1
        if GC_b:
            bnum = GC_b
            if sym:
                repeat = (bnum+1)/2
            else:
                repeat = bnum
            for i in range(repeat):
                name1 = 'GC' + str(i)
                self.coeff[name1] = reg.coef_[0][pt]
                if sym:
                    name2 = 'GC' + str(bnum-1-i)
                    self.coeff[name2] = reg.coef_[0][pt]
                pt += 1
        if Harmonic:
            name = 'Harmonic'
            self.coeff[name] = reg.coef_[0][pt]
            pt += 1
        if dPolyA:
            name = 'dPolyA'
            self.coeff[name] = {}
            for i in range(lmin, lmax+1):
                self.coeff[name][i] = {}
                if sym:
                    repeat = (self.NCPlen + 1)/2
                else:
                    repeat = self.NCPlen
                for j in range(repeat):
                    self.coeff[name][i][j] = reg.coef_[0][pt]
                    pt +=1
            
        assert len(reg.coef_[0]) == pt
        print >> sys.stderr, "Done"
        #print self.coeff
        return None

    def _predict(self, target_seq, bound, sym):
        seq_list = []
        pos_list = []
        for i in range(bound, len(target_seq)-self.NCPlen-bound):
            NCPseq = target_seq[i:i+self.NCPlen]
            seq_list.append(NCPseq)
            pos_list.append(i+self.NCPlen/2)
        
        data_dict = self.coeff
        MM_orders = []
        Kmer_k_b = []
        PolyA_b = []
        GC_b = []
        Harmonic = None
        dPolyA = []
        for key in data_dict:
            if key.startswith('MM'):
                MM_orders.append(int(key[2:]))
            if key.startswith('Kmer'):
                knum = len(data_dict[key].keys()[0])
                Kmer_k_b.append(int(key[4:]))
            if key.startswith('PolyA'):
                PolyA_b.append(int(key[5:]))
            if key.startswith('GC'):
                GC_b.append(int(key[2:]))
            if key.startswith('Harmonic'):
                Harmonic = True
            if key.startswith('dPolyA'):
                lmin = min(data_dict['dPolyA'].keys())
                lmax = max(data_dict['dPolyA'].keys())
                dPolyA = [lmin, lmax]
        if MM_orders:
            MM_orders = sorted(MM_orders)
        if Kmer_k_b:
            bnum = max(Kmer_k_b) + 1
            Kmer_k_b = ([knum, bnum])
        if PolyA_b:
            PolyA_b = max(PolyA_b) + 1
        if GC_b:
            GC_b = max(GC_b) + 1

        if MM_orders:
            MM_vars = [ [] for i in range(len(seq_list)) ] 
            for order in sorted(MM_orders):
                temp = self._var_Markov(seq_list, order, sym=sym)
                for i in range(len(temp)):
                    MM_vars[i] += temp[i]
        if Kmer_k_b:
            knum, bnum = Kmer_k_b
            Kmer_vars = self._var_Kmer(seq_list, knum, bnum, sym=sym)
        if PolyA_b:
            PolyA_vars = self._var_PolyA(seq_list, PolyA_b, sym=sym)
        if GC_b:
            GC_vars = self._var_GC(seq_list, GC_b, sym=sym)
        if Harmonic:
            Harmonic_vars = self._var_Harmonic(pos_list)
        if dPolyA:
            lmin, lmax = dPolyA
            dPolyA_vars = self._var_dPolyA(seq_list, lmin=lmin, lmax=lmax, sym=sym)

        var_list = [ [] for i in range(len(seq_list))]
        for i in range(len(seq_list)):
            if MM_orders:
                var_list[i] += MM_vars[i]
            if Kmer_k_b:
                var_list[i] += Kmer_vars[i]
            if PolyA_b:
                var_list[i] += PolyA_vars[i]
            if GC_b:
                var_list[i] += GC_vars[i]
            if Harmonic:
                var_list[i] += Harmonic_vars[i]
            if dPolyA:
                var_list[i] += dPolyA_vars[i]

        Ypred = self.reg.predict(var_list)
        Ypred = [ value[0] for value in Ypred]
        #pred_prob = [np.exp(-value) for value in Ypred]
        #pred_prob = norm(pred_prob)
        pred_score = Ypred
        return pred_score

    def predict(self, target_seq_list, bound=0, sym=True):
        pred_score_list = []
        for target_seq in target_seq_list:
            pred_score = self._predict(target_seq, bound=bound, sym=sym)
            pred_score_list.append(pred_score)
        return pred_score_list
        

    def display(self,
                data_dict,
                vmin=None,
                vmax=None):
        
        MM_orders = []
        Kmer_k_b = []
        PolyA_b = []
        GC_b = []
        Harmonic = None
        dPolyA = []
        
        for key in data_dict:
            if key.startswith('MM'):
                MM_orders.append(int(key[2:]))
            if key.startswith('Kmer'):
                knum = len(data_dict[key].keys()[0])
                Kmer_k_b.append(int(key[4:]))
            if key.startswith('PolyA'):
                PolyA_b.append(int(key[5:]))
            if key.startswith('GC'):
                GC_b.append(int(key[2:]))
            if key.startswith('Harmonic'):
                Harmonic = True
            if key.startswith('dPolyA'):
                lmin = min(data_dict['dPolyA'].keys())
                lmax = max(data_dict['dPolyA'].keys())
                dPolyA = [lmin, lmax]

        if MM_orders:
            MM_orders = sorted(MM_orders)
        if Kmer_k_b:
            bnum = max(Kmer_k_b) + 1
            Kmer_k_b = ([knum, bnum])
        if PolyA_b:
            PolyA_b = max(PolyA_b) + 1
        if GC_b:
            GC_b = max(GC_b) + 1

        if MM_orders:
            for order in sorted(MM_orders):
                name = 'MM' + str(order)
                img = []
                nts = all_path(order+1, 'ATCG')
                labels = []
                for nt in nts:
                    labels.append(nt)
                    row = []
                    for i in range(self.NCPlen - order):
                        row.append(data_dict[name][i][nt])
                    img.append(row)
                fig = plt.figure()
                plt.imshow(img, interpolation='none', aspect='auto', vmin=vmin, vmax=vmax)
                plt.colorbar()
                plt.show()
                plt.close()

                fig = plt.figure()
                for i in range(len(img)):
                    row = img[i]
                    label = labels[i]
                    x_axis = [ i - len(row)*0.5 + 0.5  for i in range(len(row))]
                    plt.plot(x_axis, row, label=label)
                plt.legend()
                plt.show()
                plt.close()
                
        if 1 in MM_orders:
            name = 'MM1'
            AT = ['AA', 'AT', 'TA', 'TT']
            GC = ['GG', 'GC', 'CG', 'CC']
            Y1, Y2 = np.zeros(self.NCPlen-1), np.zeros(self.NCPlen-1)
            for i in range(self.NCPlen - 1):
                for key in AT:
                    Y1[i] += data_dict[name][i][key] / len(AT)
                for key in GC:
                    Y2[i] += data_dict[name][i][key] / len(GC)
            fig = plt.figure()
            x_axis = [ i - len(Y1)*0.5 + 0.5  for i in range(len(Y1))]
            plt.plot(x_axis, Y1, label='AA/AT/TA/TT')
            plt.plot(x_axis, Y2, label='GG/GC/CG/CC')
            plt.legend()
            plt.show()
            plt.close()

        if Kmer_k_b:
            knum, bnum = Kmer_k_b
            nts = all_path(knum, 'ATCG')                

            #for nt in nts:
            #    row = []
            #    for i in range(bnum):
            #        name = 'Kmer' + str(i)
            #        row.append(data_dict[name][nt])
            #    img.append(row)
            #print img
            #fig = plt.figure()
            #plt.imshow(img, interpolation='none', aspect='auto')
            #plt.colorbar()
            #plt.show()
            #plt.close()

            nt_value_list = []
            for i in range(bnum):
                nt_value = []
                for nt in nts:
                    name = 'Kmer' + str(i)
                    nt_value.append((nt,data_dict[name][nt]))
                nt_value_list.append(nt_value)

            for i in range(bnum):
                nt_value = nt_value_list[i]
                X, Y = [], []
                for j in range(len(nt_value)):
                    X.append(nt_value[j][0])
                    Y.append(nt_value[j][1])
                fig = plt.figure()
                plt.plot(range(len(Y)), Y, '.')
                plt.show()
                plt.close()

        if PolyA_b:
            bnum = PolyA_b
            Y = []
            for i in range(bnum):
                name = 'PolyA' + str(i)
                Y.append(data_dict[name])
            fig = plt.figure()
            plt.plot(Y)
            plt.xlabel('Bin')
            plt.show()
            plt.close()

        if GC_b:
            bnum = GC_b
            Y = []
            for i in range(bnum):
                name = 'GC' + str(i)
                Y.append(data_dict[name])
            fig = plt.figure()
            plt.plot(Y)
            plt.xlabel('Bin')
            plt.show()
            plt.close()

        if Harmonic:
            k = data_dict['Harmonic']
            X = np.asarray(range(self.templatelen))
            Y = [0.5*k*(x-self.templatelen/2)**2 for x in X]
            fig = plt.figure()
            plt.plot(X,Y)
            plt.show()
            plt.close()

        if dPolyA:
            name = 'dPolyA'
            lmin, lmax = dPolyA
            img = np.zeros((lmax-lmin+1, self.NCPlen))
            for i in range(lmin, lmax+1):
                for j in range(self.NCPlen/2 + 1):
                    img[i-lmin][j] = data_dict[name][i][j]
                    if j < self.NCPlen/2:
                        img[i-lmin][self.NCPlen-1-j] = data_dict[name][i][j]
            fig = plt.figure()
            plt.imshow(img, interpolation='none', aspect='auto', vmin=vmin, vmax=vmax)
            plt.colorbar()
            plt.show()
            plt.close()

            fig = plt.figure()
            for i in range(len(img)):
                row = img[i]
                plt.plot(row, label=str(i+lmin))
            plt.legend()
            plt.show()
            plt.close()
            
            self.img = img
            #fig = plt.figure()
            #plt.plot(row, label=str(i+lmin))
            #plt.legend()
            #plt.show()
            #plt.close()

        if self.shape:
            names = ["MGW", "HelT", "ProT", "Roll"]
            fig = plt.figure()
            for name in names:
                freq = data_dict[name]
                plt.plot(freq, label=name)
            plt.legend()
            plt.show()
            plt.close()

class LinSliderModel:
    def __init__(self, key_slider, NCPlen=147, templatelen=225, bound=0, shape=False):
        assert NCPlen % 2 != 0
        #self.key_slider = key_slider
        self.NCPlen = NCPlen
        self.templatelen = templatelen
        self.bound = bound
        self.nlogprob_list, self.seq_list, self.pos_list = [], [], []
        self.MGW_list, self.HelT_list, self.ProT_list, self.Roll_list = [], [], [], []
        self.key_list = sorted(key_slider.keys())
        for key in self.key_list:
            slider = key_slider[key]
            energy_profile = slider.energy_profile()
            #energy_profile = -np.log(slider.dyadmap)
            seq = slider.seq
            for i in range(self.NCPlen/2+bound, self.templatelen-self.NCPlen/2-bound):
                nlogprob = energy_profile[i]
                NCPseq = seq[i-self.NCPlen/2:i+self.NCPlen/2+1]
                self.nlogprob_list.append(nlogprob)
                self.seq_list.append(NCPseq)
                self.pos_list.append(i)
                if shape:
                    self.MGW_list.append(slider.MGW[i-self.NCPlen/2:i+self.NCPlen/2+1])
                    self.HelT_list.append(slider.HelT[i-self.NCPlen/2:i+self.NCPlen/2+1])
                    self.ProT_list.append(slider.ProT[i-self.NCPlen/2:i+self.NCPlen/2+1])
                    self.Roll_list.append(slider.Roll[i-self.NCPlen/2:i+self.NCPlen/2+1])
                
        self.reg = None
        self.coeff = None
        self.rsquare = None
        self.corr = None
        self.freq = None
        if shape:
            self.shape = True
        else:
            self.shape = False

    def _even_sampling(self, sym=True, shape=False):
        if not shape:
            seq_samples = []
            for i in range(len(self.seq_list)):
                seq = self.seq_list[i]
                seq_samples.append(seq)
                if sym:
                    seq_samples.append(rev_comp(seq))
            return seq_samples
        if shape:
            MGW_samples, HelT_samples, ProT_samples, Roll_samples = [], [], [], []
            for i in range(len(self.seq_list)):
                MGW_samples.append(self.MGW_list[i])
                HelT_samples.append(self.HelT_list[i])
                ProT_samples.append(self.ProT_list[i])
                Roll_samples.append(self.Roll_list[i])
                if sym:
                    MGW_samples.append(self.MGW_list[i][::-1])
                    HelT_samples.append(self.HelT_list[i][::-1])
                    ProT_samples.append(self.ProT_list[i][::-1])
                    Roll_samples.append(self.Roll_list[i][::-1])
            return [MGW_samples, HelT_samples, ProT_samples, Roll_samples]


    def _bias_sampling(self, scale=100, sym=True, shape=False):
        if not shape:
            seq_samples = []
            for i in range(len(self.seq_list)):
                seq = self.seq_list[i]
                nlogprob = self.nlogprob_list[i]
                prob = np.exp(-nlogprob)
                count = int(round(prob*scale))
                for i in range(count):
                    seq_samples.append(seq)
                    if sym:
                        seq_samples.append(rev_comp(seq))
            return seq_samples
        if shape:
            MGW_samples, HelT_samples, ProT_samples, Roll_samples = [], [], [], []
            for i in range(len(self.seq_list)):
                nlogprob = self.nlogprob_list[i]
                prob = np.exp(-nlogprob)
                count = int(round(prob*scale))
                for i in range(count):
                    MGW_samples.append(self.MGW_list[i])
                    HelT_samples.append(self.HelT_list[i])
                    ProT_samples.append(self.ProT_list[i])
                    Roll_samples.append(self.Roll_list[i])
                    if sym:
                        MGW_samples.append(self.MGW_list[i][::-1])
                        HelT_samples.append(self.HelT_list[i][::-1])
                        ProT_samples.append(self.ProT_list[i][::-1])
                        Roll_samples.append(self.Roll_list[i][::-1])

            return [MGW_samples, HelT_samples, ProT_samples, Roll_samples]

             
    def _stat_Markov(self, seq_list, order):
        ntdic = {}
        for nt in all_path(order+1, 'ATCG'):
            ntdic[nt] = 0.0

        sample_num = len(seq_list)
        
        freq = [ copy.deepcopy(ntdic) for i in range(self.NCPlen - order) ]
        for seq in seq_list:
            for i in range(len(seq) - order):
                nt = seq[i:i+1+order]
                freq[i][nt] += 1.0 / sample_num
        
        mean, std = [], []
        for ntdic in freq:
            mean.append(np.mean(ntdic.values()))
            std.append(np.std(ntdic.values()))

        stdz_freq = []
        for i in range(len(freq)):
            ntdic = freq[i]
            temp = {}
            for nt, value in ntdic.items():
                temp[nt] = (value - mean[i]) / std[i]
            stdz_freq.append(temp)

        return freq, sample_num, mean, std, stdz_freq
        
    def _stat_Kmer(self, seq_list, knum, bnum):
        seqlen = self.NCPlen / bnum
        assert seqlen >= knum
        
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2

        ntdic = {}
        for nt in all_path(knum, 'ATCG'):
            ntdic[nt] = 0.0
        freq = [ copy.deepcopy(ntdic) for i in range(bnum) ]

        sample_num = len(seq_list)*(seqlen-knum+1)
        
        for seq in seq_list:
            for k in range(bnum):
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq = seq[st:st+seqlen]                
                for i in range(seqlen - knum + 1):
                    nt = bseq[i:i+knum]
                    freq[k][nt] += 1.0 / sample_num

        mean, std = [], []
        for ntdic in freq:
            mean.append(np.mean(ntdic.values()))
            std.append(np.std(ntdic.values()))

        stdz_freq = []
        for i in range(len(freq)):
            ntdic = freq[i]
            temp = {}
            for nt, value in ntdic.items():
                temp[nt] = (value - mean[i]) / std[i]
            stdz_freq.append(temp)

        return freq, sample_num, mean, std, stdz_freq

    def _stat_PolyA(self, seq_list, bnum):
        seqlen = self.NCPlen / bnum
        assert seqlen >= 1
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2

        freq = [ [] for i in range(bnum) ]
        sample_num = len(seq_list)
        
        for seq in seq_list:
            for k in range(bnum):
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq = seq[st:st+seqlen]
                num_pos = Amer_len(bseq)
                score = 0.0
                for num, pos in num_pos.items():
                    if num >= 5:
                        score += num*len(pos)
                freq[k].append(score)
                
        mean, std = [], []
        for k in range(len(freq)):
            mean.append(np.mean(freq[k]))
            std.append(np.std(freq[k]))

        return freq, sample_num, mean, std

    def _stat_GC(self, seq_list, bnum):
        seqlen = self.NCPlen / bnum
        assert seqlen >= 1
        
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2

        freq = [ [] for i in range(bnum) ]

        sample_num = len(seq_list)
        
        for seq in seq_list:
            for k in range(bnum):
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq = seq[st:st+seqlen]
                GC = GC_content(bseq)
                freq[k].append(GC)

        mean, std = [], []
        for k in range(len(freq)):
            mean.append(np.mean(freq[k]))
            std.append(np.std(freq[k]))
        return freq, sample_num, mean, std

    def _stat_shape(self, shape_list):
        shape_list = np.asarray(shape_list)
        #print shape_list.shape
        freq = np.sum(shape_list, axis=0) / len(shape_list)
        return freq

    def report (self,
                MM_orders,
                Kmer_k_b,
                PolyA_b,
                GC_b,
                Harmonic,
                scale=100,
                shape=False,
                sym=True,
                graph=False):
        freq = {}
        if shape:
            print >> sys.stderr, "sampling data values"        
            even_samples_list = self._even_sampling(sym=sym, shape=True)
            bias_samples_list = self._bias_sampling(scale=scale, sym=sym, shape=True)
            names = ["MGW", "HelT", "ProT", "Roll"] 
            for i in range(len(even_samples_list)):
                even_samples, bias_samples = even_samples_list[i], bias_samples_list[i]
                freq[names[i]] = self._stat_shape(bias_samples) / self._stat_shape(even_samples)
                self.freq = freq
        else:
            # sampling data
            print >> sys.stderr, "sampling data values"        
            even_samples = self._even_sampling(sym=sym)
            bias_samples = self._bias_sampling(scale=scale, sym=sym)
        
        # frequency counts
        print >> sys.stderr, "counting sequence features"
        if MM_orders:
            for order in sorted(MM_orders):
                name = 'MM' + str(order)
                freq1, sample_num, mean, std, stdz_freq = self._stat_Markov(even_samples, order)
                freq2, sample_num, mean, std, stdz_freq = self._stat_Markov(bias_samples, order)
                freq_fold = [{} for i in range(self.NCPlen-order)]
                nts = all_path(order+1, 'ATCG')
                for i in range(self.NCPlen-order):
                    for nt in nts:
                        freq_fold[i][nt] = freq2[i][nt] / freq1[i][nt]
                freq[name] = freq_fold        
        if Kmer_k_b:
            knum, bnum = Kmer_k_b
            freq1, sample_num, mean, std, stdz_freq = self._stat_Kmer(even_samples, knum, bnum)
            freq2, sample_num, mean, std, stdz_freq = self._stat_Kmer(bias_samples, knum, bnum)
            nts = all_path(knum, 'ATCG')
            for i in range(bnum):
                name = 'Kmer' + str(i)
                freq_fold = {}
                for nt in nts:
                    freq_fold[nt] = freq2[i][nt] / freq1[i][nt]
                freq[name] = freq_fold
        if PolyA_b:
            bnum = PolyA_b
            freq1, sample_num, mean1, std = self._stat_PolyA(even_samples, bnum)
            freq2, sample_num, mean2, std = self._stat_PolyA(bias_samples, bnum)
            for i in range(bnum):
                name = 'PolyA' + str(i)
                freq[name] = mean2[i]/mean1[i]
        if GC_b:
            bnum = GC_b
            freq1, sample_num, mean1, std = self._stat_GC(even_samples, bnum)
            freq2, sample_num, mean2, std = self._stat_GC(bias_samples, bnum)
            for i in range(bnum):
                name = 'GC' + str(i)
                freq[name] = mean2[i]/mean1[i]
        if Harmonic:
            None

            # To do
        self.freq = freq
        print >> sys.stderr, "Done"
        return None
    

    def _var_Markov (self, seq_list, order):
        nt_pos = {}
        nts = all_path(order+1, 'ATCG')
        for i in range(len(nts)):
            nt = nts[i]
            nt_pos[nt] = i
        if order % 2 == 0:
            palnum = 0
        else:
            palnum = 4**((order+1)/2)
        var_list = []
        for seq in seq_list:
            left = seq
            right = rev_comp(seq)
            row = []
            for i in range((len(seq) - order)/2):
                nt1 = left[i:i+order+1]
                nt2 = right[i:i+order+1]
                temp = [0.0] * len(nts)
                temp[nt_pos[nt1]] += 1
                temp[nt_pos[nt2]] += 1
                row += temp
            if (len(seq) - order) % 2 != 0:
                assert order % 2 == 0
                i = len(seq)/2 - (order+1)/2
                nt = left[i:i+order+1]
                temp = [0.0] * ((len(nts)+palnum)/2)
                pos = min(nt_pos[nt], nt_pos[rev_comp(nt)])
                temp[pos] += 1
                row += temp
            var_list.append(row)
        return var_list

    def _var_Kmer (self, seq_list, knum, bnum):        
        seqlen = self.NCPlen / bnum
        assert seqlen >= knum        
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2
        nt_pos = {}
        nts = all_path(knum, 'ATCG')
        for i in range(len(nts)):
            nt = nts[i]
            nt_pos[nt] = i
        if knum % 2 != 0:
            palnum = 0
        else:
            palnum = 4**(knum/2)
        var_list = []
        for seq in seq_list:
            left = seq
            right = rev_comp(seq)
            row = []
            for k in range(bnum/2):
                temp = [0]*len(nts)
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq1 = left[st:st+seqlen]
                bseq2 = right[st:st+seqlen]
                for i in range(seqlen-knum+1):
                    nt1 = bseq1[i:i+knum]
                    nt2 = bseq2[i:i+knum]
                    temp[nt_pos[nt1]] += 1
                    temp[nt_pos[nt2]] += 1
                row += temp
            if bnum % 2 != 0:
                assert seqlen % 2 != 0
                st = len(seq)/2 - seqlen/2
                bseq = left[st:st+seqlen]
                temp = [0] * ((len(nts)+palnum)/2)
                for i in range(seqlen-knum+1):
                    nt = bseq[i:i+knum]
                    pos = min(nt_pos[nt], nt_pos[rev_comp(nt)])
                    temp[pos] += 1
                row += temp
            var_list.append(row)
        return var_list

    def _var_PolyA (self, seq_list, bnum):        
        seqlen = self.NCPlen / bnum
        assert seqlen >= 1
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2        
        var_list = []
        for seq in seq_list:
            row = []
            for k in range(bnum):
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq = seq[st:st+seqlen]
                num_pos = Amer_len(bseq)
                count = 0
                for num, pos in num_pos.items():
                    if num >=5:
                        count += len(pos)
                row.append(count)
            sym_row = [row[i] + row[::-1][i] for i in range(bnum/2)]
            if bnum % 2 != 0:
                sym_row += [row[bnum/2]]
            var_list.append(sym_row)
        return var_list

    def _var_GC (self, seq_list, bnum):
        seqlen = self.NCPlen / bnum
        assert seqlen >= 1
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2        
        var_list = []
        for seq in seq_list:
            row = []
            for k in range((bnum+1)/2):
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq = seq[st:st+seqlen]
                GC = GC_content(bseq)
                row.append(GC)
            sym_row = [row[i] + row[::-1][i] for i in range(bnum/2)]
            if bnum % 2 != 0:
                sym_row += [row[bnum/2]]
            var_list.append(sym_row)
        return var_list

    def _var_Harmonic (self, pos_list):
        offset = self.templatelen / 2
        var_list = []
        for pos in pos_list:
            row = []
            row.append(0.5*((pos-offset)**2))
            var_list.append(row)
        return var_list

    def _var_dPolyA (self, seq_list, lmin=3, lmax=20):
        var_list = []
        for seq in seq_list:
            row = []
            num_pos = Amer_len(seq)
            for i in range(lmin, lmax+1):
                temp = [0.0]*len(seq)
                try:
                    pos_list = num_pos[i]
                except:
                    row += temp[:self.NCPlen/2 + 1]
                    continue
                for pos in pos_list:
                    for j in range(i):
                        temp[pos+j] += 1
                sym_temp = [temp[i] + temp[::-1][i] for i in range(self.NCPlen/2)]
                sym_temp += [temp[self.NCPlen/2]]
                row += sym_temp
            var_list.append(row)
        return var_list
    
    def train (self,
               MM_orders,
               Kmer_k_b,
               PolyA_b,
               GC_b,
               Harmonic,
               ref_key=None,
               dPolyA=False,
               alpha=0.5,
               k_fold=10,
               graph=False):

        self.MM_orders, self.Kmer_k_b, = MM_orders, Kmer_k_b
        self.PolyA_b, self.GC_b, self.Harmonic = PolyA_b, GC_b, Harmonic
        nlogprob_list, seq_list, pos_list = self.nlogprob_list, self.seq_list, self.pos_list
        
        # reading variables
        print >> sys.stderr, "reading data values"
        if MM_orders:
            MM_vars = [ [] for i in range(len(seq_list)) ] 
            for order in sorted(MM_orders):
                temp = self._var_Markov(seq_list, order)
                for i in range(len(temp)):
                    MM_vars[i] += temp[i]
        if Kmer_k_b:
            knum, bnum = Kmer_k_b
            Kmer_vars = self._var_Kmer(seq_list, knum, bnum)
        if PolyA_b:
            PolyA_vars = self._var_PolyA(seq_list, PolyA_b)
        if GC_b:
            GC_vars = self._var_GC(seq_list, GC_b)
        if Harmonic:
            Harmonic_vars = self._var_Harmonic(pos_list)
        if dPolyA:
            lmin, lmax = dPolyA
            dPolyA_vars = self._var_dPolyA(seq_list, lmin=lmin, lmax=lmax)

        var_list = [ [] for i in range(len(seq_list))]
        for i in range(len(seq_list)):
            if MM_orders:
                var_list[i] += MM_vars[i]
            if Kmer_k_b:
                var_list[i] += Kmer_vars[i]
            if PolyA_b:
                var_list[i] += PolyA_vars[i]
            if GC_b:
                var_list[i] += GC_vars[i]
            if Harmonic:
                var_list[i] += Harmonic_vars[i]
            if dPolyA:
                var_list[i] += dPolyA_vars[i]


        # adjust variables by reference key
        if ref_key:
            new_var_list = []
            new_nlogprob_list = []
            knum = self.key_list.index(ref_key)
            gnum = self.templatelen - self.NCPlen + 1 - 2*self.bound
            st, ed = knum*gnum, (knum+1)*gnum
            ref_var = np.asarray(var_list[st:ed])
            ref_nlogprob = np.asarray(nlogprob_list[st:ed])
            ref_prob = [np.exp(-value) for value in ref_nlogprob]
            #fig = plt.figure()
            #plt.plot(ref_prob)
            
            for i in range(len(self.key_list)):
                if i == knum:
                    continue
                var_group = np.asarray(var_list[i*gnum:(i+1)*gnum])
                nlogprob_group = np.asarray(nlogprob_list[i*gnum:(i+1)*gnum])
                new_var_list += list(var_group - ref_var)
                new_nlogprob_list += list((nlogprob_group - ref_nlogprob)/ref_nlogprob)
                #print self.key_list[i]
                #plt.title(self.key_list[i])
                #prob = [np.exp(-value) for value in nlogprob_group]
                #plt.plot(ref_prob)
                #plt.plot(prob)
                #plt.show()
            var_list = new_var_list
            nlogprob_list = new_nlogprob_list

        
        # adjust variables by reference point
        gnum = self.templatelen - self.NCPlen + 1 - 2*self.bound
        if gnum > 1:
            feature_list, target_list = [], []
            i = 0
            count = 0
            while i < len(var_list):
                ref_var = copy.deepcopy(var_list[i])
                ref_nlogprob = copy.deepcopy(nlogprob_list[i])
                j = 1
                count += 1
                while j < gnum:
                    row = np.asarray(var_list[i+j]) - np.asarray(ref_var)
                    feature_list.append(row)
                    target_list.append([nlogprob_list[i+j] - ref_nlogprob])
                    j +=1
                    count += 1
                i += j
            assert count == len(var_list)
        else:
            feature_list = var_list
            target_list = [[value] for value in nlogprob_list]

        
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
        pt = 0
        self.coeff = {}
        if MM_orders:
            for order in MM_orders:
                name = 'MM' + str(order)
                self.coeff[name] = [{} for i in range(self.NCPlen-order)]
                nts = all_path(order+1, 'ATCG')
                if order % 2 == 0:
                    palnum = 0
                else:
                    palnum = 4**((order+1)/2)
                for k in range((self.NCPlen - order + 1)/2):
                    count = 0
                    if order % 2 == 0 and k == (self.NCPlen - order + 1)/2 - 1:
                        total = (len(nts)+palnum)/2
                    else:
                        total = len(nts)
                    while count < total:
                        nt = nts[count]
                        self.coeff[name][k][nt] = reg.coef_[0][pt]
                        self.coeff[name][self.NCPlen-order-1-k][rev_comp(nt)] = reg.coef_[0][pt]
                        count += 1
                        pt += 1
        if Kmer_k_b:
            knum, bnum = Kmer_k_b
            if knum % 2 != 0:
                palnum = 0
            else:
                palnum = 4**(knum/2)
            nts = all_path(knum, 'ATCG')
            for i in range((bnum+1)/2):
                name1 = 'Kmer' + str(i)
                name2 = 'Kmer' + str(bnum-1-i)
                self.coeff[name1] = {}
                self.coeff[name2] = {}
                count = 0
                if bnum % 2 != 0 and i == (bnum+1)/2 - 1:
                    total = (len(nts)+palnum)/2
                else:
                    total = len(nts)
                while count < total:
                    nt = nts[count]
                    self.coeff[name1][nt] = reg.coef_[0][pt]
                    self.coeff[name2][rev_comp(nt)] = reg.coef_[0][pt]
                    count += 1
                    pt += 1
        if PolyA_b:
            bnum = PolyA_b
            for i in range((bnum+1)/2):
                name1 = 'PolyA' + str(i)
                name2 = 'PolyA' + str(bnum-1-i)
                self.coeff[name1] = reg.coef_[0][pt]
                self.coeff[name2] = reg.coef_[0][pt]
                pt += 1
        if GC_b:
            bnum = GC_b
            for i in range((bnum+1)/2):
                name1 = 'GC' + str(i)
                name2 = 'GC' + str(bnum-1-i)
                self.coeff[name1] = reg.coef_[0][pt]
                self.coeff[name2] = reg.coef_[0][pt]
                pt += 1
        if Harmonic:
            name = 'Harmonic'
            self.coeff[name] = reg.coef_[0][pt]
            pt += 1
        if dPolyA:
            name = 'dPolyA'
            self.coeff[name] = {}
            for i in range(lmin, lmax+1):
                self.coeff[name][i] = {}
                for j in range(self.NCPlen/2 + 1):
                    self.coeff[name][i][j] = reg.coef_[0][pt]
                    pt +=1
            
        assert len(reg.coef_[0]) == pt
        print >> sys.stderr, "Done"
        #print self.coeff
        return None

    def _predict(self, seq, bound=0):
        seq_list = []
        pos_list = []
        for i in range(self.NCPlen/2+bound, self.templatelen-self.NCPlen/2-bound):
            NCPseq = seq[i-self.NCPlen/2:i+self.NCPlen/2+1]
            seq_list.append(NCPseq)
            pos_list.append(i)
        
        data_dict = self.coeff
        MM_orders = []
        Kmer_k_b = []
        PolyA_b = []
        GC_b = []
        Harmonic = None
        dPolyA = []
        for key in data_dict:
            if key.startswith('MM'):
                MM_orders.append(int(key[2:]))
            if key.startswith('Kmer'):
                knum = len(data_dict[key].keys()[0])
                Kmer_k_b.append(int(key[4:]))
            if key.startswith('PolyA'):
                PolyA_b.append(int(key[5:]))
            if key.startswith('GC'):
                GC_b.append(int(key[2:]))
            if key.startswith('Harmonic'):
                Harmonic = True
            if key.startswith('dPolyA'):
                lmin = min(data_dict['dPolyA'].keys())
                lmax = max(data_dict['dPolyA'].keys())
                dPolyA = [lmin, lmax]
        if MM_orders:
            MM_orders = sorted(MM_orders)
        if Kmer_k_b:
            bnum = max(Kmer_k_b) + 1
            Kmer_k_b = ([knum, bnum])
        if PolyA_b:
            PolyA_b = max(PolyA_b) + 1
        if GC_b:
            GC_b = max(GC_b) + 1

        if MM_orders:
            MM_vars = [ [] for i in range(len(seq_list)) ] 
            for order in sorted(MM_orders):
                temp = self._var_Markov(seq_list, order)
                for i in range(len(temp)):
                    MM_vars[i] += temp[i]
        if Kmer_k_b:
            knum, bnum = Kmer_k_b
            Kmer_vars = self._var_Kmer(seq_list, knum, bnum)
        if PolyA_b:
            PolyA_vars = self._var_PolyA(seq_list, PolyA_b)
        if GC_b:
            GC_vars = self._var_GC(seq_list, GC_b)
        if Harmonic:
            Harmonic_vars = self._var_Harmonic(pos_list)
        if dPolyA:
            lmin, lmax = dPolyA
            dPolyA_vars = self._var_dPolyA(seq_list, lmin=lmin, lmax=lmax)

        var_list = [ [] for i in range(len(seq_list))]
        for i in range(len(seq_list)):
            if MM_orders:
                var_list[i] += MM_vars[i]
            if Kmer_k_b:
                var_list[i] += Kmer_vars[i]
            if PolyA_b:
                var_list[i] += PolyA_vars[i]
            if GC_b:
                var_list[i] += GC_vars[i]
            if Harmonic:
                var_list[i] += Harmonic_vars[i]
            if dPolyA:
                var_list[i] += dPolyA_vars[i]

        Ypred = self.reg.predict(var_list)
        Ypred = [ value[0] for value in Ypred]
        pred_prob = [np.exp(-value) for value in Ypred]
        pred_prob = norm(pred_prob)
        return pred_prob

    def predict(self, key_slider, scale=2, bound=0):
        keys = random.sample(key_slider.keys(), 30)
        key_slider2 = {}
        for key in keys:
            #win,loc = key.split('-')
            #loc = int(loc)
            #st,ed = loc, loc+len(win)
            seq = key_slider[key].seq
            KDE = norm(key_slider[key].dyadmap)
            prob = list(np.asarray(self._predict(seq, bound=bound))**scale)
            prob = norm(prob)
            prob = ([0.0]*(self.NCPlen/2)) +  prob  + ([0.0]*(self.NCPlen/2))
            #fig = plt.figure()
            #plt.plot(KDE, label='exp')
            #plt.plot(prob, label='pre')
            #plt.axvspan(st, ed-1, alpha=0.5, color='red')
            #plt.legend()
            #plt.show()
            #plt.close()
            key_slider2[key] = Slider(key,'','','','',seq,prob,'','')
        sample_list = [keys]
        graph_edit.plot_map(key_slider, sample_list, norm_choice=True, obs_func = Slider.get_dyadmap, draw = 'polyA', slicing=0, note='Exp')
        graph_edit.plot_map(key_slider2, sample_list, norm_choice=True, obs_func = Slider.get_dyadmap, draw = 'polyA', slicing=0, note='Pred')

        #sample_list = [[key] for key in keys]
        #graph_edit.plot_signal(key_slider, sample_list, norm_choice=True, obs_func = Slider.get_dyadmap, draw = 'polyA', slicing=0, note='Exp')
        #graph_edit.plot_signal(key_slider2, sample_list, norm_choice=True, obs_func = Slider.get_dyadmap, draw = 'polyA', slicing=0, note='Pred')

        return None
        

    def display(self,
                data_dict,
                vmin=None,
                vmax=None):

        MM_orders = []
        Kmer_k_b = []
        PolyA_b = []
        GC_b = []
        Harmonic = None
        dPolyA = []
        

        for key in data_dict:
            if key.startswith('MM'):
                MM_orders.append(int(key[2:]))
            if key.startswith('Kmer'):
                knum = len(data_dict[key].keys()[0])
                Kmer_k_b.append(int(key[4:]))
            if key.startswith('PolyA'):
                PolyA_b.append(int(key[5:]))
            if key.startswith('GC'):
                GC_b.append(int(key[2:]))
            if key.startswith('Harmonic'):
                Harmonic = True
            if key.startswith('dPolyA'):
                lmin = min(data_dict['dPolyA'].keys())
                lmax = max(data_dict['dPolyA'].keys())
                dPolyA = [lmin, lmax]

        if MM_orders:
            MM_orders = sorted(MM_orders)
        if Kmer_k_b:
            bnum = max(Kmer_k_b) + 1
            Kmer_k_b = ([knum, bnum])
        if PolyA_b:
            PolyA_b = max(PolyA_b) + 1
        if GC_b:
            GC_b = max(GC_b) + 1

        if MM_orders:
            for order in sorted(MM_orders):
                name = 'MM' + str(order)
                img = []
                nts = all_path(order+1, 'ATCG')
                labels = []
                for nt in nts:
                    labels.append(nt)
                    row = []
                    for i in range(self.NCPlen - order):
                        row.append(data_dict[name][i][nt])
                    img.append(row)
                fig = plt.figure()
                plt.imshow(img, interpolation='none', aspect='auto', vmin=vmin, vmax=vmax)
                plt.colorbar()
                plt.show()
                plt.close()

                fig = plt.figure()
                for i in range(len(img)):
                    row = img[i]
                    label = labels[i]
                    x_axis = [ i - len(row)*0.5 + 0.5  for i in range(len(row))]
                    plt.plot(x_axis, row, label=label)
                plt.legend()
                plt.show()
                plt.close()
                
        if 1 in MM_orders:
            name = 'MM1'
            AT = ['AA', 'AT', 'TA', 'TT']
            GC = ['GG', 'GC', 'CG', 'CC']
            Y1, Y2 = np.zeros(self.NCPlen-1), np.zeros(self.NCPlen-1)
            for i in range(self.NCPlen - 1):
                for key in AT:
                    Y1[i] += data_dict[name][i][key] / len(AT)
                for key in GC:
                    Y2[i] += data_dict[name][i][key] / len(GC)
            fig = plt.figure()
            x_axis = [ i - len(Y1)*0.5 + 0.5  for i in range(len(Y1))]
            plt.plot(x_axis, Y1, label='AA/AT/TA/TT')
            plt.plot(x_axis, Y2, label='GG/GC/CG/CC')
            plt.legend()
            plt.show()
            plt.close()

        if Kmer_k_b:
            knum, bnum = Kmer_k_b
            nts = all_path(knum, 'ATCG')                

            #for nt in nts:
            #    row = []
            #    for i in range(bnum):
            #        name = 'Kmer' + str(i)
            #        row.append(data_dict[name][nt])
            #    img.append(row)
            #print img
            #fig = plt.figure()
            #plt.imshow(img, interpolation='none', aspect='auto')
            #plt.colorbar()
            #plt.show()
            #plt.close()

            nt_value_list = []
            for i in range(bnum):
                nt_value = []
                for nt in nts:
                    name = 'Kmer' + str(i)
                    nt_value.append((nt,data_dict[name][nt]))
                nt_value_list.append(nt_value)

            for i in range(bnum):
                nt_value = nt_value_list[i]
                X, Y = [], []
                for j in range(len(nt_value)):
                    X.append(nt_value[j][0])
                    Y.append(nt_value[j][1])
                fig = plt.figure()
                plt.plot(range(len(Y)), Y, '.')
                plt.show()
                plt.close()

        if PolyA_b:
            bnum = PolyA_b
            Y = []
            for i in range(bnum):
                name = 'PolyA' + str(i)
                Y.append(data_dict[name])
            fig = plt.figure()
            plt.plot(Y)
            plt.xlabel('Bin')
            plt.show()
            plt.close()

        if GC_b:
            bnum = GC_b
            Y = []
            for i in range(bnum):
                name = 'GC' + str(i)
                Y.append(data_dict[name])
            fig = plt.figure()
            plt.plot(Y)
            plt.xlabel('Bin')
            plt.show()
            plt.close()

        if Harmonic:
            k = data_dict['Harmonic']
            X = np.asarray(range(self.templatelen))
            Y = [0.5*k*(x-self.templatelen/2)**2 for x in X]
            fig = plt.figure()
            plt.plot(X,Y)
            plt.show()
            plt.close()

        if dPolyA:
            name = 'dPolyA'
            lmin, lmax = dPolyA
            img = np.zeros((lmax-lmin+1, self.NCPlen))
            for i in range(lmin, lmax+1):
                for j in range(self.NCPlen/2 + 1):
                    img[i-lmin][j] = data_dict[name][i][j]
                    if j < self.NCPlen/2:
                        img[i-lmin][self.NCPlen-1-j] = data_dict[name][i][j]
            fig = plt.figure()
            plt.imshow(img, interpolation='none', aspect='auto', vmin=vmin, vmax=vmax)
            plt.colorbar()
            plt.show()
            plt.close()

            fig = plt.figure()
            for i in range(len(img)):
                row = img[i]
                plt.plot(row, label=str(i+lmin))
            plt.legend()
            plt.show()
            plt.close()
            
            self.img = img
            #fig = plt.figure()
            #plt.plot(row, label=str(i+lmin))
            #plt.legend()
            #plt.show()
            #plt.close()

        if self.shape:
            names = ["MGW", "HelT", "ProT", "Roll"]
            fig = plt.figure()
            for name in names:
                freq = data_dict[name]
                plt.plot(freq, label=name)
            plt.legend()
            plt.show()
            plt.close()
