import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import seaborn as sns
import copy
import math
from scipy import signal

def tuple_cmp(a, b):
    if a[0] <= b[0]:
        return -1
    else:
        return 1

def norm(L):
    total = sum(L)
    return [L[i]/float(total) for i in range(len(L))]

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

def moving_average (signal, win):
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

def draw_along_genome (ID_pos, ID_sig_list, win, labels, ylabel, scatt=False, note=""):
    posID = []
    for ID in ID_pos:
        pos = ID_pos[ID]
        posID.append([pos, ID])
    posID = sorted(posID, cmp=tuple_cmp)
    X = []
    Y_list = [[] for i in range(len(ID_sig_list))]
    for pos, ID in posID:
        X.append(pos)
        for i in range(len(ID_sig_list)):
            ID_sig = ID_sig_list[i]
            sig = ID_sig[ID]
            Y_list[i].append(sig)
    fig = plt.figure(figsize=(10,5))
    if scatt:
        for i in range(len(Y_list)):
            Y = Y_list[i]
            plt.plot(X, Y, 'k.', alpha=0.2)
    for i in range(len(Y_list)):
        Y = signal.fftconvolve(Y_list[i], np.ones((win,))/win, mode='same')
        plt.plot(X, Y, label = labels[i])
    plt.xlabel("Genomic coordinate (bp)")
    plt.ylabel(ylabel)
    plt.legend()
    #plt.savefig("Gwide_" + note + ".png", bbox_inches='tight')
    plt.show()
    plt.close()

def draw_along_genome_pair (ID_pos, ID_sig1, ID_sig2,  win, ylabel1, ylabel2, xlabel="Genomic coordinate (bp)", note=""):
    posID = []
    for ID in ID_pos:
        pos = ID_pos[ID]
        posID.append([pos, ID])
    posID = sorted(posID, cmp=tuple_cmp)
    X = []
    Y1, Y2 = [], []
    for pos, ID in posID:
        X.append(pos)
        Y1.append(ID_sig1[ID])
        Y2.append(ID_sig2[ID])
    Y1 = signal.fftconvolve(Y1, np.ones((win,))/win, mode='same')
    Y2 = signal.fftconvolve(Y2, np.ones((win,))/win, mode='same')
    fig, ax1 = plt.subplots(figsize=(10,5))
    ax1.plot(X, Y1, 'b', alpha=0.5)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel1, color='b')
    ax1.tick_params('y', colors='b')
    ax2 = ax1.twinx()
    ax2.plot(X, Y2, 'r', alpha=0.5)
    ax2.set_ylabel(ylabel2, color='r')
    ax2.tick_params('y', colors='r')
    fig.tight_layout()
    #plt.savefig("Gwide_pair_" + note + ".png", bbox_inches='tight')
    plt.show()
    plt.close()


def Scatter_plot (ID_xvalue, ID_score, xlim=[0, 100], ylim=[-2, 2.5], markersize=None, xlabel="AT content (%)", ylabel = 'Condensability (A.U.)', note=""):
    X, Y = [], []
    xvalue_scores = {}
    for ID in ID_xvalue:
        xvalue, score = ID_xvalue[ID], ID_score[ID]
        X.append(xvalue)
        Y.append(score)
        if xvalue not in xvalue_scores:
            xvalue_scores[xvalue] = []
        xvalue_scores[xvalue].append(score)
    Xmean, Ymean, Yerr = [], [], []
    for xvalue in xvalue_scores:
        if len(xvalue_scores[xvalue]) <= 1:
            continue
        Xmean.append(xvalue)
        Ymean.append(np.mean(xvalue_scores[xvalue]))
        Yerr.append(np.std(xvalue_scores[xvalue]/np.sqrt(len(xvalue_scores[xvalue]))))
    fig = plt.figure()
    plt.plot(X, Y, 'k.', alpha=0.01, markersize=markersize)
    plt.errorbar(Xmean, Ymean, yerr=Yerr, fmt='.')
    plt.plot(Xmean, Ymean,'.')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #plt.title(chr_name)
    #plt.ylim([0,4])
    #plt.ylim([-2,2.5])
    #plt.xlim([0, 100])
    plt.xlim(xlim)
    plt.ylim(ylim)
    leg = plt.legend(loc='best', numpoints=1, prop={'size': 6})
    for lh in leg.legendHandles:
        lh._legmarker.set_markersize(15)
        lh._legmarker.set_alpha(1)
    plt.savefig(xlabel + '_vs_' + ylabel + "_" + note + "_scatter.png")
    #plt.show()
    plt.close()

def PartitionBoxplot (ID_score, ID_target, frac, xlabel="", ylabel = 'Condensability - GC dependence (A.U.)', note=""):
    frac = sorted(norm(frac), reverse=True)
    group = quantile(ID_target, len(frac), frac=frac)
    boxdata = []
    for IDs in group:
        #print len(IDs)
        temp = []
        for ID in IDs:
            temp.append(ID_score[ID])
        boxdata.append(temp)
    fig = plt.figure()
    plt.xlabel('Partitions by ' + xlabel)
    plt.ylabel(ylabel)
    plt.boxplot(boxdata, 0, "")
    plt.savefig(xlabel + '_vs_' + ylabel + "_" + note + "_pbox.png")
    #plt.show()
    plt.close()

def PartitionScatterplot (ID_xvalue, ID_score, ID_target, frac, xlim=[0, 100], ylim=[-2,2.5], xlabel="AT content (%)", ylabel = 'Condensability (A.U.)', note=""):
    frac = sorted(norm(frac), reverse=True)
    group = quantile(ID_target, len(frac), frac=frac)
    X_list, Y_list = [], []
    label_list = []
    for i in range(len(group)):
        IDs = group[i]
        X = [ID_xvalue[ID] for ID in IDs]
        Y = [ID_score[ID] for ID in IDs]
        label = str(i+1) + 'th partition ' + '(' + str(len(IDs)) + ')'
        X_list.append(X)
        Y_list.append(Y)
        label_list.append(label)
    alpha_list = np.linspace(0.02, 1, num=len(group))
    color_list = np.linspace(0.01, 1, num=len(group))
    cmap = mpl.cm.get_cmap("jet")
    fig = plt.figure()
    for i in range(len(X_list)):
        X, Y = X_list[i], Y_list[i]
        plt.plot(X, Y, '.', color=cmap(color_list[i]), alpha=alpha_list[i], label=label_list[i])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #plt.ylim([-2,2.5])
    #plt.xlim([0, 100])
    plt.xlim(xlim)
    plt.ylim(ylim)
    leg = plt.legend(loc='best', numpoints=1, prop={'size': 6})
    for lh in leg.legendHandles:
        lh._legmarker.set_markersize(15)
        lh._legmarker.set_alpha(1)
    plt.savefig(xlabel + '_vs_' + ylabel + "_" + note + "_pscatter.png")
    #plt.show()
    plt.close()

def PartitionMeanplot (ID_xvalue, ID_score, ID_target, frac, xlim=[0,100], ylim=[-0.8, 0.6], xlabel="AT content (%)", ylabel = 'Condensability (A.U.)', note=""):
    frac = sorted(norm(frac), reverse=True)
    group = quantile(ID_target, len(frac), frac=frac)
    X_list, Y_list, Z_list = [], [], []
    label_list = []
    for i in range(len(group)):
        IDs = group[i]
        xvalue_scores = {}
        for ID in IDs:
            xvalue = ID_xvalue[ID]
            if xvalue not in xvalue_scores:
                xvalue_scores[xvalue] = []
            xvalue_scores[xvalue].append(ID_score[ID])
        X, Y, Z = [], [], []
        for xvalue in xvalue_scores:
            if len (xvalue_scores[xvalue]) <= 1:
                continue
            X.append(xvalue)
            Y.append(np.mean(xvalue_scores[xvalue]))
            Z.append(np.std(xvalue_scores[xvalue]/np.sqrt(len(xvalue_scores[xvalue]))))
        X_list.append(X)
        Y_list.append(Y)
        Z_list.append(Z)
        label = str(i+1) + 'th partition ' + '(' + str(len(IDs)) + ')'
        label_list.append(label)
    alpha_list = np.linspace(0.02, 1, num=len(group))
    color_list = np.linspace(0.01, 1, num=len(group))
    cmap = mpl.cm.get_cmap("jet")
    fig = plt.figure()
    for i in range(len(X_list)):
        X, Y, Z = X_list[i], Y_list[i], Z_list[i]
        plt.errorbar(X, Y, yerr=Z, fmt='.', color=cmap(color_list[i]))
        plt.plot(X, Y, '.', color=cmap(color_list[i]), label=label_list[i])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #plt.ylim([-0.8, 0.6])
    #plt.xlim([0, 100])
    plt.xlim(xlim)
    plt.ylim(ylim)
    leg = plt.legend(loc='best', numpoints=1, prop={'size': 6})
    for lh in leg.legendHandles:
        lh._legmarker.set_markersize(15)
        lh._legmarker.set_alpha(1)
    plt.savefig(xlabel + '_vs_' + ylabel + "_" + note + "_pmean.png")
    #plt.show()
    plt.close()
