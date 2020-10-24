import numpy as np
import math
import copy
import matplotlib.pyplot as plt
from sklearn import linear_model
from scipy import signal
from scipy.fftpack import fft, ifft, ifftshift, fftshift

def norm(L):
    total = sum(L)
    return [L[i]/float(total) for i in range(len(L))]

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
        new_sig.append(np.mean(neighbors))
    return new_sig


def moving_average (input_signal, win, interpolate=True):
    if np.isnan(input_signal).any():
        if interpolate:
            input_signal = NN_interpolate(input_signal)
    return signal.fftconvolve(input_signal, np.ones((win,))/win, mode='same')


def neutralize_score_by_target (ID_score, ID_target):
    def standardization (data_list):
        if len(data_list) <= 1:
            return data_list
        mean = np.mean(data_list)
        std = np.std(data_list)
        #return [ float(data - mean) for data in data_list]
        return [ float(data - mean)/std for data in data_list]
    target_scores = {}
    target_IDs = {}
    for ID in ID_score:
        target = ID_target[ID]
        score = ID_score[ID]
        if target not in target_scores:
            target_scores[target] = []
            target_IDs[target] = []
        target_scores[target].append(score)
        target_IDs[target].append(ID)
    new_ID_score = {}
    for target in target_scores:
        new_scores = standardization(target_scores[target])
        IDs = target_IDs[target]
        for i in range(len(IDs)):
            ID = IDs[i]
            new_score = new_scores[i]
            new_ID_score[ID] = new_score
    return new_ID_score

def partial_corr (X, Y, Zs):
    if len(Zs) <= 0:
        return get_corr(X, Y)
    Z = Zs.pop()
    return (partial_corr(X,Y,Zs) - partial_corr(X,Z,Zs)*partial_corr(Z,Y,Zs)) / (np.sqrt(1-partial_corr(X,Z,Zs)**2)*np.sqrt(1-partial_corr(Z,Y,Zs)**2))


def semi_partial_corr (X, Y, Zs):
    if len(Zs) <= 0:
        return get_corr(X, Y)
    Z = Zs.pop()
    return (semi_partial_corr(X,Y,Zs) - semi_partial_corr(X,Z,Zs)*semi_partial_corr(Z,Y,Zs)) / (np.sqrt(1-semi_partial_corr(X,Z,Zs)**2))


def neutralize (ID_score, ID_target_list, graph=False):
    ID_newscore = {}
    IDs = ID_score.keys()
    Y, X = [], []
    for ID in IDs:
        Y.append([ID_score[ID]])
        x = [ID_target[ID] for ID_target in ID_target_list]
        X.append(x)
    reg = linear_model.Ridge(alpha=0.5)
    reg.fit (X, Y)
    Ypred = reg.predict(X)

    if graph and len(ID_target_list) == 1:
        fig = plt.figure()
        plt.plot([ x[0] for x in X], [ y[0] for y in Y ], '.', alpha=0.1)
        plt.plot([ x[0] for x in X], [ y[0] for y in Ypred], '--')
        plt.show()
        plt.close()
        
    for i in range(len(IDs)):
        ID = IDs[i]
        newscore = Y[i][0] - Ypred[i][0]
        ID_newscore[ID] = newscore
    return ID_newscore

def standardize (data_list):
    if len(data_list) <= 1:
        return data_list
    mean = np.mean(data_list)
    std = np.std(data_list)
    return [ float(data - mean)/std for data in data_list]

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
        
