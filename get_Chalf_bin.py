import glob
import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy
import numpy as np
import random
import pickle
from scipy.optimize import curve_fit
from sklearn import linear_model

def sigmoid(x, L ,x0, k):
    y = L / (1 + np.exp(k*(x-x0)))
    return (y)

def read_titration (fname):
    tnum_conc = {}
    tnum_frac = {}
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        conc, frac, tnum = cols[0], cols[7], cols[-1]
        try:
            tnum = int(tnum)
        except:
            continue
        conc = float(conc)
        frac = float(frac)
        tnum_conc[tnum] = conc
        tnum_frac[tnum] = frac
    return tnum_conc, tnum_frac

def read_bin_num (fname, skip_star=False, chr_choice=None):
    chr_binID_range = {}
    chr_binID_tnum_count = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        if First:
            tnums = [int(name.rsplit('.', 1)[0].split('-')[-1]) for name in cols[4:-1]]
            First = False
            continue
        if cols[-1] == '*':
            if skip_star:
                continue
            cols = cols[:-1]
        _, chr, st, ed = cols[:4]
        if chr_choice and chr not in chr_choice:
            continue
        st, ed = int(st), int(ed)
        counts = [int(count) for count in cols[4:]]
        if sum(counts) <=0:
            continue        
        if chr not in chr_binID_range:
            chr_binID_range[chr] = []
        chr_binID_range[chr].append((st, ed))
        if chr not in chr_binID_tnum_count:
            chr_binID_tnum_count[chr] = []
        tnum_count = {tnum:count for tnum, count in zip(tnums, counts)}
        chr_binID_tnum_count[chr].append(tnum_count)

    return chr_binID_range, chr_binID_tnum_count


# parameters
path = "/home/spark159/../../media/spark159/sw/"

# sample information
#cell = 'H1'
#cell = 'GM'
#cell = 'mCD8T'
#sample = 'NCP'
#sample = 'WT-NCP'
#sample = 'inht-NCP'
#sample = 'KO-NCP'
#agent = 'sp'

# experiment list (cell, sample, agent)
exp_list = [('H1', 'NCP', 'sp'),
            ('H1', 'NCP', 'spd'),
            ('H1', 'NCP', 'CoH'),
            ('H1', 'NCP', 'PEG'),
            ('H1', 'NCP', 'Ca'),
            ('H1', 'NCP', 'Mg'),
            ('H1', 'NCP', 'HP1a'),
            ('H1', 'NCP', 'HP1bSUV'),
            ('H1', 'NCP', 'LKH'),
            ('H1', 'NCP', 'Ki67'),
            ('H1', 'NCP', 'FUS')]

bin_size = 10000
draw_graph = False

for cell, sample, agent in exp_list:

    # set species and gender
    if cell in ['H1', 'GM']:
        species = 'human'
    elif cell in ['mCD8T']:
        species = 'mouse'

    if cell in ['H1']:
        gender = 'male'
    elif cell in ['GM', 'mCD8T']:
        gender = 'female'

    # set chromosome list
    if species == 'human':
        chr_list = ['chr' + str(i) for i in range(1, 23)]
    elif species == 'mouse':
        chr_list = ['chr' + str(i) for i in range(1, 20)]
    chr_list += ['chrX']

    if gender == 'male':
        chr_list += ['chrY']

    chr_list = ['chr1'] # temporal


    # file name
    nfname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000))+'kb']) +'_num.cn'
    tfname = '_'.join([cell, sample, agent, 'titration']) + '.csv'

    print "working on %s" % (nfname)

    # load data
    tnum_conc, tnum_frac = read_titration(tfname)
    chr_binID_range, chr_binID_tnum_count = read_bin_num(nfname, chr_choice=chr_list)

    # fitting titration curve and extract C-half
    rsq_list = []
    chr_binID_Chalf = {}
    total_count, nan_count = 0, 0
    for chr in chr_list:
        for binID in range(len(chr_binID_tnum_count[chr])):
            tnum_count = chr_binID_tnum_count[chr][binID]
            control = tnum_count[0]
            if control <= 0:
                continue
            X, Y = [], []
            for tnum in sorted(tnum_count):
                if tnum == 0:
                    #control = tnum_count[0]
                    X.append(0.0)
                    Y.append(1.0)
                    continue
                conc = tnum_conc[tnum]
                frac = float(tnum_count[tnum])/control
                X.append(conc)
                Y.append(frac)

            try:
                p0 = [max(Y), np.median(X), 1]
                bounds = ([0.0, 0.0, 0.0], [max(Y)+max(Y)*0.1, np.inf, np.inf])
                popt, pcov = curve_fit(sigmoid, X, Y, p0, bounds = bounds,  method='dogbox')
                residuals = np.asarray(Y)- sigmoid(X, *popt)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((np.asarray(Y)-np.mean(Y))**2)
                r_squared = 1 - (ss_res / ss_tot)
                rsq_list.append(r_squared)
                Chalf = popt[1]
            except:
                Chalf = np.nan
                nan_count +=1

            if r_squared < 0.5:
                Chalf = np.nan
                nan_count +=1

            if chr not in chr_binID_Chalf:
                chr_binID_Chalf[chr] = []
            chr_binID_Chalf[chr].append(Chalf)
            total_count +=1

            if draw_graph:
                fig = plt.figure()
                plt.plot(X, Y, '.', markersize=10, alpha=0.2)
                if np.isnan(Chalf):
                    plt.show()
                    plt.close()
                    continue
                pred_X = np.linspace(min(X), max(X), 1000)
                pred_Y = sigmoid(pred_X, *popt)
                plt.plot(pred_X, pred_Y, 'k-', alpha=0.2)
                plt.xlabel("Concentration")
                plt.ylabel("Soluble fraction")
                if agent in ['HP1a']:
                    plt.xscale('log', basex=2)
                elif agent in ['sp', 'spd', 'CoH']:
                    plt.xscale('log', basex=10)
                plt.title("%s binID:%d" % (chr, binID))
                #plt.ylim([0, 1.2])
                plt.show()
                plt.close()

    print "fitting failure %d/%d" % (nan_count, total_count)
    print "fitting failure %f percentage" % (100*float(nan_count)/total_count)


    # writing Chalf file
    fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb']) + '_Chalf.cn'
    f = open(fname, 'w')
    s = 'BinID\tChromosome\tStart\tEnd\tChalf'
    print >> f, s

    for chr in chr_list:
        for binID in range(len(chr_binID_Chalf[chr])):
            Chalf = chr_binID_Chalf[chr][binID]
            if np.isnan(Chalf):
                continue
            st, ed = chr_binID_range[chr][binID]
            print >> f, '%d\t%s\t%d\t%d\t%f' % (binID, chr, st, ed, Chalf)

    f.close()
    print
