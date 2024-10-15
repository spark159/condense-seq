import sys
import numpy as np
import math
import copy
import matplotlib.pyplot as plt
from sklearn import linear_model
from scipy import signal
from scipy.fftpack import fft, ifft, ifftshift, fftshift
import sklearn.cluster
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import scipy.stats
import Interval_dict

# get histogram
def get_hist (data, binnum=1000, prob=False):
    hist={};
    if prob:
        deno=float(len(data))
    else:
        deno=1.0
    binwidth=float(max(data)-min(data))/binnum
    for value in data:
        bin=int((value-min(data))/binwidth)
        bincenter=min(data)+(bin+0.5)*binwidth
        if bincenter not in hist:
            hist[bincenter]=0
        hist[bincenter]+=1/deno
    return hist

# flatten the nested data into a list
def flatten (data,
             target_types=[list, set, tuple]):
    if type(data) not in target_types:
        return [data]
    new_data = []
    for element in data:
        new_data +=flatten(element)
    return new_data

# all possible k-mers with given nts
def all_kmer(k, nts='ATCG'):
    if k==1:
        return list(nts)
    output=[]
    for kmer in all_kmer(k-1):
        for nt in nts:
            output.append(kmer+nt)
    return output


# compute GC content of sequence
def GC_content(seq):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    return (num/float(len(seq)))*100


def print_Roman(num):
 
    # Storing roman values of digits from 0-9
    # when placed at different places
    m = ["", "M", "MM", "MMM"]
    c = ["", "C", "CC", "CCC", "CD", "D",
         "DC", "DCC", "DCCC", "CM "]
    x = ["", "X", "XX", "XXX", "XL", "L",
         "LX", "LXX", "LXXX", "XC"]
    i = ["", "I", "II", "III", "IV", "V",
         "VI", "VII", "VIII", "IX"]
 
    # Converting to roman
    thousands = m[num // 1000]
    hundreds = c[(num % 1000) // 100]
    tens = x[(num % 100) // 10]
    ones = i[num % 10]
 
    ans = (thousands + hundreds +
           tens + ones)
 
    return ans


# find a reverse complementary of sequence
def rev_comp (seq):
    rev_seq=seq[::-1]; new_seq=''
    for nt in rev_seq:
        nt = nt.upper()
        if nt == 'A':
            new_seq += 'T'
        elif nt == 'T':
            new_seq += 'A'
        elif nt == 'C':
            new_seq += 'G'
        elif nt == 'G':
            new_seq += 'C'
        else:
            new_seq += nt
    return new_seq


# check the sequence is palindromic
def is_pal (seq):
    if len(seq) % 2 != 0:
        return False
    for i in range(len(seq)/2):
        if nt[i] != rev_comp(nt[len(seq)-1-i]):
            return False
    return True


# find the longest length of poly-nt (pos=False)
# find all poly-nt locations and counts (pos=True)
def polynt_count (seq, nts, pos=False):
    len_pos = {}
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
            if count not in len_pos:
                len_pos[count] = []
            len_pos[count].append(i)
            i = j
        else:
            i +=1
    if pos:
        return len_pos    
    if len(len_pos) == 0:
        return 0
    return max(len_pos.keys())


# find the count of the given dinucleotide
# find the count of all existing dinucleotide (din=None)
# check both top and bottom strands (both=True)
def get_dincount(seq, din=None, both=False):
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
        if both:
            din = rev_comp(din)
            if din not in din_count:
                din_count[din] = 0
            din_count[din] += 1
    return din_count


def get_kmer_count(seq, klen, both=False):
    kmer_count = {}
    seq = seq.upper()
    seq_list = seq.split('N')
    for seq in seq_list:
        if len(seq) < klen:
            continue
        if both:
            rev_seq = rev_comp(seq)
        for i in range(len(seq)-klen-1):
            kmer = seq[i:i+klen]
            if kmer not in kmer_count:
                kmer_count[kmer] = 0
            kmer_count[kmer] += 1
            if both:
                kmer = rev_seq[i:i+klen]
                if kmer not in kmer_count:
                    kmer_count[kmer] = 0
                kmer_count[kmer] += 1
    return kmer_count


def get_fract_dict (ID_test,
                    ID_control,
                    IDs=None,
                    dummy_addon=0,
                    div_error=np.NaN):

    if IDs == None:
        IDs = set(ID_test.keys()) | set(ID_control.keys())

    ID_fract = {}
    for ID in IDs:
        try:
            test = float(ID_test[ID])
        except:
            test = 0.0
        try:
            control = float(ID_control[ID])
        except:
            control = 0.0
            
        test += dummy_addon
        control += dummy_addon

        if control <= 0:
            ID_fract[ID] = div_error
            continue

        fract = test / control
        ID_fract[ID] = fract
    return ID_fract

def get_fract_profile (ID_tests,
                       ID_controls,
                       IDs=None,
                       dummy_addon=0,
                       div_error=np.NaN):

    if IDs == None:
        IDs = set(ID_tests.keys()) | set(ID_controls.keys())
        IDs = sorted(list(IDs))
        
    data_len = len(ID_tests[IDs[0]])

    ID_fracts = {}
    for ID in IDs:
        fracts = []
        for i in range(data_len):
            try:
                test = ID_tests[ID][i]
            except:
                test = 0.0
            try:
                control = ID_controls[ID][i]
            except:
                control = 0.0

            test += dummy_addon
            control += dummy_addon

            if control <= 0:
                fracts.append(div_error)
                continue

            fract = test / control
            fracts.append(fract)

        ID_fracts[ID] = fracts
    return ID_fracts
    
                       
def stat_Markov(seq_list, NCPlen, order):
    ntdic = {}
    for nt in all_kmer(order+1, 'ATCG'):
        ntdic[nt] = 0.0

    sample_num = len(seq_list)

    freq = [ copy.deepcopy(ntdic) for i in range(NCPlen - order) ]
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


def get_ATGC_sig (freq_MM1):
    length = len(freq_MM1)
    nts = all_kmer(2, 'ATCG')
    AT_sig, GC_sig = np.zeros(length), np.zeros(length)
    for nt in nts:
        row = [ freq_MM1[i][nt] for i in range(length)]
        if nt in ['AA', 'AT', 'TA', 'TT']:
            AT_sig += np.asarray(row)
        if nt in ['GG', 'GC', 'CG', 'CC']:
            GC_sig += np.asarray(row)
    for i in range(length):
        AT_sig[i] = float(AT_sig[i]) / sum(freq_MM1[i].values())
        GC_sig[i] = float(GC_sig[i]) / sum(freq_MM1[i].values())
    return AT_sig, GC_sig


def stat_Kmer(seq_list, NCPlen, knum, bnum):
    seqlen = NCPlen / bnum
    assert seqlen >= knum

    extra = NCPlen % bnum
    boundoff = extra / 2
    centeroff = extra % 2

    ntdic = {}
    for nt in all_kmer(knum, 'ATCG'):
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


# rescale the data in old range (old_st, old_ed) into new range (new_st, new_ed)
def rescale (value,
             old_st,
             old_ed,
             new_st,
             new_ed,
             adjust_outsider=False):

    if value < old_st:
        if adjust_outsider:
            return new_st
        else:
            raise ValueError('value is outside of range')
    elif value > old_ed:
        if adjust_outsider:
            return new_ed
        else:
            raise ValueError('value is outside of range')
    else:
        return new_st + (new_ed - new_st)*float(value-old_st)/(old_ed-old_st)

# rescale the data list
def rescale_list (value_list,
                  old_st,
                  old_ed,
                  new_st,
                  new_ed,
                  adjust_outsider=False):
    output = []
    for value in value_list:
        new_value = rescale(value,
                            old_st,
                            old_ed,
                            new_st,
                            new_ed,
                            adjust_outsider=adjust_outsider)
        
        output.append(new_value)
    return output

def norm(L):
    total = sum(L)
    return [L[i]/float(total) for i in range(len(L))]

def standardize (data_list,
                 by_mean=True,
                 by_std=True):
    if by_mean:
        mean = np.nanmean(data_list)
    else:
        mean = 0.0
    if by_std:
        std = np.nanstd(data_list)
    else:
        std = 1.0
    return [float(data - mean)/std for data in data_list]

def standardize_dict (key_value,
                      keys=None,
                      by_mean=True,
                      by_std=True):
    if keys == None:
        keys = key_value.keys()
        
    if by_mean:
        mean = np.nanmean([key_value[key] for key in keys])
    else:
        mean = 0.0
    if by_std:
        std = np.nanstd([key_value[key] for key in keys])
    else:
        std = 1.0
        
    key_zscore = {}
    for key in keys:
        value = key_value[key]
        zscore = float(value - mean)/std
        key_zscore[key] = zscore
    return key_zscore

def get_corr(x, y):
    assert len(x) == len(y)
    x, y = np.asarray(x), np.asarray(y)
    selected = (~np.isnan(x))*(~np.isnan(y))
    x, y = x[selected], y[selected]
    n = len(x)
    if n <= 0:
        return np.nan
    #assert n > 0
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

def get_spearman_corr (X, Y):
    assert len(X) == len(Y)
    X, Y = np.asarray(X), np.asarray(Y)
    selected = (~np.isnan(X))*(~np.isnan(Y))
    X, Y = X[selected], Y[selected]
    if len(X) < 2:
        return np.nan
    return scipy.stats.spearmanr(X, Y)[0]


def acf(x):
    #xp = ifftshift((x - np.average(x))/np.std(x))
    xp = fftshift((x - np.average(x))/np.std(x))
    n, = xp.shape
    xp = np.r_[xp[:n//2], np.zeros_like(xp), xp[n//2:]] # zero-padding
    f = fft(xp)
    p = np.absolute(f)**2
    pi = ifft(p)
    return np.real(pi)[:n//2]/(np.arange(n//2)[::-1]+n//2)


def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size/2:]


def slow_acf(x):
    output = []
    mean = np.mean(x)
    std = np.std(x)
    for lag in range(len(x)-10):
        temp = 0.0
        for i in range(len(x)-lag):
            temp += (x[i]-mean)*(x[i+lag]-mean)
        temp = float(temp)/(len(x)-lag)
        temp = temp/(std**2)
        output.append(temp)
    return output


def test(x):
    output = []
    for lag in range(len(x)-10):
        temp = 0.0
        for i in range(len(x)-lag):
            temp += abs(x[i]-x[i+lag])
        temp = float(temp) / (len(x) -lag)
        output.append(temp)
    return output


def quantile_dict (ID_score, num, IDs=None, frac=None):
    def value_cmp(a, b):
        if a[1] <= b[1]:
            return -1
        else:
            return 1
        
    if IDs == None:
        IDs = ID_score.keys()

    IDscore = [[ID, ID_score[ID]] for ID in IDs]
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
    #ed = 0
    for i in range(num):
        size = size_list[i]
        st = i*size
        ed = min((i+1)*size,len(IDscore))
        #size = size_list[i]
        #st = ed             
        #ed = st + size
        temp = [IDscore[j][0] for j in range(st,ed)]
        output.append(temp)
    return output

# partition data into given windows of value
def partition (ID_value,
               p_wins):
    p_num = len(p_wins)
    p_IDs = [[] for i in range(p_num)]
    for ID in ID_value:
        value = ID_value[ID]
        for i in range(p_num):
            st, ed = p_wins[i]
            if value >= st and value < ed:
                break
        p_IDs[i].append(ID)
    return p_IDs
        

def NN_interpolate (raw_data_list):
    data_list = copy.deepcopy(raw_data_list)
    has_data = False
    # handle first one
    for i in range(len(data_list)):
        if not np.isnan(data_list[i]):
            data_list[:i] = [data_list[i]]*i
            has_data = True
            break
    # no data
    if not has_data:
        return data_list
    # handle extra
    while i < len(data_list)-1:
        j = i + 1
        while j < len(data_list)-1:
            if not np.isnan(data_list[j]):
                break
            j += 1
        mid = int(math.ceil((i+j)/2.0))
        #print i, mid, j
        data_list[i+1:mid] = [data_list[i]]*(mid-i-1)
        if not np.isnan(data_list[j]):
            data_list[mid:j] = [data_list[j]]*(j-mid)
        else:
            data_list[mid:j+1] = [data_list[i]]*(j-mid+1)
        i = j
    return data_list


def NN_interpolate_dict (ID_data):
    ID_newdata = {}
    for ID, data in ID_data.items():
        ID_newdata[ID] = NN_interpolate(data)
    return ID_newdata


def slow_moving_average (signal, win):
    assert win % 2 != 0
    new_sig = []
    for i in range(win/2):
        new_sig.append(0)
    for i in range(win/2, len(signal)-win/2):
        value = sum(signal[i-win/2:i+win/2+1])/float(win)
        new_sig.append(value)
    for i in range(win/2):
        new_sig.append(0)
    return new_sig

def slow_moving_average2 (signal, win):
    assert win % 2 != 0
    new_sig = []
    for i in range(len(signal)):
        neighbors = np.asarray(signal[max(0, i-win/2):min(i+win/2+1, len(signal))])
        neighbors = neighbors[~np.isnan(neighbors)]
        if len(neighbors) <= 0:
            new_sig.append(np.nan)
        else:
            new_sig.append(np.mean(neighbors))
    return new_sig


def moving_average (input_signal,
                    win,
                    interpolate=True):

    if np.isnan(input_signal).any():
        if interpolate:
            input_signal = NN_interpolate(input_signal)
    return signal.fftconvolve(input_signal, np.ones((win,))/win, mode='same')

def binning (data_list, bin_num):
    offset = min(data_list)
    bin_size = int(math.ceil(float(max(data_list)-min(data_list)+1) / bin_num))
    idx_list = []
    for i in range(len(data_list)):
        data = data_list[i]
        idx = int(round(float(data - offset)/bin_size))
        assert idx >=0 and idx < bin_num
        idx_list.append(idx)
    return idx_list

# categorize genomic positions into states falling in the window
def categorize (state_intervals,
                ID_pos,
                hash_func=None,
                max_pos=None,
                domain_size=None,
                silent=False):

    # if hash function not provided
    if hash_func == None:
        dID_interval = {}
        for state in state_intervals:        
            intervals = state_intervals[state]
            for i in range(len(intervals)):
                dID = state + ':' + str(i)
                assert dID not in dID_interval
                dID_interval[dID] = intervals[i]
        
        state_dict = Interval_dict.double_hash(dID_interval,
                                               max_pos=max_pos,
                                               domain_size=domain_size,
                                               silent=silent)
    else:
        state_dict = hash_func

    state_IDs = {state:[] for state in state_intervals}
    for ID in ID_pos:
        pos = ID_pos[ID]
        find_dIDs = state_dict.find(pos)
        for dID in find_dIDs:
            state = dID.split(':')[0]
            state_IDs[state].append(ID)

    return state_IDs


# categorize genomic bins into states with maximum overlap
def categorize_bin (state_intervals,
                    binID_interval,
                    hash_func=None,
                    max_pos=None,
                    domain_size=None,
                    silent=False):

    # if hash function not provided
    if hash_func == None:
        dID_interval = {}
        for state in state_intervals:
            intervals = state_intervals[state]
            for i in range(len(intervals)):
                dID = state + ':' + str(i)
                assert dID not in dID_interval
                dID_interval[dID] = intervals[i]

        state_dict = Interval_dict.double_hash(dID_interval,
                                               max_pos=max_pos,
                                               domain_size=domain_size,
                                               silent=silent)

    else:
        state_dict = hash_func

    state_binIDs = {state:[] for state in state_intervals}
    for binID in binID_interval:
        st, ed = binID_interval[binID]
        find_dIDs = state_dict.insert(st, ed, 1)
        dID_value = state_dict.get()

        weight_state = []
        for dID in find_dIDs:
            weight = dID_value[dID]
            state = dID.split(':')[0]
            weight_state.append((weight, state))
        weight_state = sorted(weight_state, reverse=True)

        best_weight = weight_state[0][0]
        for weight, state in weight_state:
            if weight < best_weight:
                break
            state_binIDs[state].append(binID)

        state_dict.clear()                

    return state_binIDs


# categorize regular genomic bins (binstep/binsize) into states with maximum overlap
def categorize_rbin (state_intervals,
                     bin_size,
                     bin_step,
                     hash_func=None,
                     max_pos=None,
                     binID_interval=None,
                     silent=False):

    # if hash function not provided
    if hash_func == None:
        bin_dict = Interval_dict.bin_hash(bin_size,
                                          bin_step,
                                          max_pos=max_pos,
                                          ID_interval=binID_interval,
                                          silent=silent)
    else:
        bin_dict = hash_func
        
    binID_WeightStates = {}
    for state in state_intervals:
        for interval in state_intervals[state]:
            st, ed = interval
            find_binIDs = bin_dict.insert_range(st, ed, 1)
            binID_value = bin_dict.get()
            for binID in find_binIDs:
                weight = binID_value[binID]
                if binID not in binID_WeightStates:
                    binID_WeightStates[binID] = []
                binID_WeightStates[binID].append((weight, state))
            bin_dict.clear()

    state_binIDs = {state:[] for state in state_intervals}
    for binID in binID_WeightStates:
        WeightStates = sorted(binID_WeightStates[binID], reverse=True)
        best_weight = WeightStates[0][0]

        for weight, state in WeightStates:
            if weight < best_weight:
                break
            state_binIDs[state].append(binID)

    return state_binIDs

### binning average the data in genomic bins
def bin_data_mean (binID_interval,
                   ID_loc,
                   ID_value,
                   hash_func=None,
                   max_pos=None,
                   domain_size=None,
                   min_sample_size=1,
                   skip_nan=False,
                   silent=False):

    # if hash function not provided
    if hash_func == None:
        Int_dict = Interval_dict.double_hash(binID_interval,
                                             domain_size=domain_size,
                                             max_pos=max_pos,
                                             silent=silent)
    else:
        Int_dict = hash_func

    Int_dict_sum = copy.deepcopy(Int_dict)
    Int_dict_count = copy.deepcopy(Int_dict)

    for ID in ID_value:
        value = ID_value[ID]
        loc = ID_loc[ID]

        try:
            st, ed = loc
            Int_dict_sum.insert_range(st, ed, value)
            Int_dict_count.insert_range(st, ed, 1)
        except:
            pos = loc
            Int_dict_sum.insert(pos, value)
            Int_dict_count.insert(pos, 1)

    binID_sum = Int_dict_sum.get()
    binID_count = Int_dict_count.get()

    binID_mean = {}
    for binID in binID_interval:
        try:
            total = binID_sum[binID]
            count = binID_count[binID]
            if count < min_sample_size:
                mean = np.nan
            else:
                mean = float(total) / count
        except:
            mean = np.nan        
        if skip_nan and np.isnan(mean):
            continue
        binID_mean[binID] = mean
        
    return binID_mean

### binning average the data in regular bins
def rbin_data_mean (bin_size,
                    bin_step,
                    ID_loc,
                    ID_value,
                    binID_interval=None,
                    hash_func=None,
                    max_pos=None,
                    min_sample_size=1,
                    skip_nan=False,
                    silent=False):
    
    # if hash function not provided
    if hash_func == None:
        Int_dict = Interval_dict.bin_hash(bin_size,
                                          bin_size,
                                          max_pos=max_pos,
                                          ID_interval=binID_interval,
                                          silent=silent)
    else:
        Int_dict = hash_func

    Int_dict_sum = copy.deepcopy(Int_dict)
    Int_dict_count = copy.deepcopy(Int_dict)

    for ID in ID_value:
        value = ID_value[ID]
        loc = ID_loc[ID]

        try:
            st, ed = loc
            Int_dict_sum.insert_range(st, ed, value)
            Int_dict_count.insert_range(st, ed, 1)
        except:
            pos = loc
            Int_dict_sum.insert(pos, value)
            Int_dict_count.insert(pos, 1)

    binID_sum = Int_dict_sum.get()
    binID_count = Int_dict_count.get()

    if binID_interval != None:
        binIDs = binID_interval.keys()
    else:
        binIDs = range(0, max_pos / bin_step + 1)

    binID_mean = {}
    for binID in binIDs:
        try:
            total = binID_sum[binID]
            count = binID_count[binID]
            if count < min_sample_size:
                mean = np.nan
            else:
                mean = float(total) / count
        except:
            mean = np.nan        
        if skip_nan and np.isnan(mean):
            continue
        binID_mean[binID] = mean        
    return binID_mean

### Spectral clustering
def Spectral_clustering (score_matrix,
                         cluster_num):
    idx_cID = sklearn.cluster.spectral_clustering(affinity=score_matrix,
                                                  n_clusters=cluster_num,
                                                  random_state=0)
    cID_idxs = {}
    for i in range(len(idx_cID)):
        cID = idx_cID[i]
        if cID not in cID_idxs:
            cID_idxs[cID] = []
        cID_idxs[cID].append(i)
    return idx_cID, cID_idxs

### Hierarchial clustering
def Hierarchial_clustering (dist_matrix):
    y = squareform(dist_matrix)
    Z = linkage(y, 'ward', optimal_ordering=True)
    #Z = linkage(y, 'average', optimal_ordering=False)

    #idx_cID = [cID-1 for cID in fcluster(Z, t=0.01, criterion='distance')]
    idx_cID = [cID-1 for cID in fcluster(Z, t=3, criterion='maxclust')]
    cID_idxs = {}
    for i in range(len(idx_cID)):
        cID = idx_cID[i]
        if cID not in cID_idxs:
            cID_idxs[cID] = []
        cID_idxs[cID].append(i)
    return Z, idx_cID, cID_idxs

### OPTICS clustering
def OPTICS_clustering (dist_matrix):
    idx_cID = sklearn.cluster.OPTICS(min_samples=2,
                                     metric='precomputed').fit_predict(dist_matrix)

    cID_idxs = {}
    for i in range(len(idx_cID)):
        cID = idx_cID[i]
        if cID not in cID_idxs:
            cID_idxs[cID] = []
        cID_idxs[cID].append(i)
    return idx_cID, cID_idxs

### decode Z (linkage) information
def decode_Z (Z, idx_key):
    node_children = {i:{} for i in range(len(idx_key))}
    node_dist = {i:None for i in range(len(idx_key))}
    node_keys = {i:{idx_key[i]} for i in range(len(idx_key))}
    for i in range(len(Z)):
        node1, node2 = int(Z[i][0]), int(Z[i][1])
        new_node = max(node_keys.keys()) + 1
        node_children[new_node] = set([node1, node2])
        node_dist[new_node] = float(Z[i][2])
        node_keys[new_node] = node_keys[node1] | node_keys[node2]
    return node_children, node_dist, node_keys

### encode Z (linkage)
def encode_Z (node_children, node_dist, node_keys):
    idx_node = []
    for node in node_children:
        idx_node.append(node)
    idx_node = sorted(idx_node)

    node_idx = {}
    for i in range(len(idx_node)):
        node = idx_node[i]
        node_idx[node] = i

    dist_row = []
    for node in node_children:
        if len(node_children[node]) <= 0:
            continue
        cnode1, cnode2 = sorted(list(node_children[node]))
        dist, size = node_dist[node], len(node_keys[node])
        row = [node_idx[cnode1], node_idx[cnode2], dist, size]
        dist_row.append((dist, row))
    dist_row = sorted(dist_row)
    Z = np.asarray([row for dist, row in dist_row])
    return idx_node, Z

### get Cohen's kappa
def get_kappa (bvec1, bvec2):
    assert len(bvec1) == len(bvec2)
    bpair_count = {(0,0):0, (1,0):0, (0,1):0, (1,1):0}
    for i in range(len(bvec1)):
        bpair = (bvec1[i], bvec2[i])
        bpair_count[bpair] +=1

    sum0_ = bpair_count[(0,0)] + bpair_count[(0,1)]
    sum1_ = bpair_count[(1,0)] + bpair_count[(1,1)]
    sum_0 = bpair_count[(0,0)] + bpair_count[(1,0)]
    sum_1 = bpair_count[(0,1)] + bpair_count[(1,1)]
    total = sum(bpair_count.values())
    obs = float(bpair_count[(1,1)]+bpair_count[(0,0)])/total
    exp = float(sum_1*sum1_+sum_0*sum0_)/(total*total)
    kappa = float(obs-exp)/(1-exp)
    assert kappa <= 1
    return kappa


# correlation analysis
def correlate (sig1,
               sig2,
               max_dist=sys.maxint,
               circular=False,
               clip=10):
    
    sig1 = sig1[clip:len(sig1)-clip]
    sig2 = sig2[clip:len(sig2)-clip]
    sig1 = np.asarray(sig1) - np.mean(sig1)
    sig2 = np.asarray(sig2) - np.mean(sig2)
    dist_products = {}
    for i in range(len(sig1)):
        for j in range(i, min(len(sig1), i+max_dist+1)):
            dist = j - i
            product = sig1[i]*sig2[j]
            if dist not in dist_products:
                dist_products[dist] = []
            dist_products[dist].append(product)
    corr_sig = [0.0]*(max(dist_products.keys())+1)
    for dist, products in dist_products.items():
        corr_sig[dist] = np.mean(products)
    return corr_sig


def FFT (sig, clip=10):
    sig = sig[clip:len(sig)-clip]
    N = len(sig)
    sig_ft = fft(sig)[1:N/2]
    periods = [float(N)/k for k in range(1, N/2)]
    amplts = np.abs(sig_ft)/float(N)
    phases = np.arctan2(sig_ft.imag, sig_ft.real) / np.pi # unit of pi
    shifts = np.asarray([(phases[k-1] * N) / (2*np.pi*k) for k in range(1, N/2)]) /np.pi # unit of pi
    return periods, amplts, phases

