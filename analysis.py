from termcolor import colored
import os, sys, subprocess, re
import copy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn.cluster
import sklearn.covariance
#import graph
import SliderClass
import math
import random


ref_length = 225

global dyad_axis
dyad_axis = (1+ref_length)/2 -1


def GC_content(seq):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    return (num/float(len(seq)))*100

def Amer_len_detail(seq, pos=True):
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


def random_map (data, N):
    output = [0]*len(data)
    for i in range(N):
        index = int(random.choice(data))
        output[index] +=1
    return output

def color_A (seq):
    text = ''
    for nt in seq:
        if nt == 'A':
            text += colored(nt, 'red')
        else:
            text += nt
    return text

def key_cmp(a, b):
    if a[0] <= b[0]:
        return -1
    else:
        return 1

def value_cmp(a, b):
    if a[1] <= b[1]:
        return -1
    else:
        return 1

def all_path(N, states='AC'):
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output

def create_map(cuts, N=ref_length, norm_choice=False):
    map=[0.0]*N
    total = float(len(cuts))
    if total <= 0:
        return map
    for loc in cuts:
        map[loc] +=1
    if norm_choice:
        map = [e / total for e in map]
    return map

def norm(L):
    total = float(sum(L))
    newL = []
    for e in L:
        newL.append(e/total)
    return newL

def sub_background (map, frac=0.1):
    thres= min(map) + frac*(max(map)-min(map))
    #new = [0 for i in range(len(map))]
    #for i in range(len(map)):
    #    if map[i] > thres:
    #        new[i] = map[i]
    #return new
    return [max(value-thres, 0.0) for value in map]

def find_peaks(map, back= False, num=1):
    nmap = map[:]
    if back:
        nmap = sub_background(map)

    peak_sig={}
    for i in range(1, len(nmap)-1):
        if nmap[i] > nmap[i-1] and nmap[i] > nmap[i+1]:
            peak_sig[i] = nmap[i]

    peak_sig=[[peak, sig]  for peak, sig in peak_sig.items()]
    peak_sig=sorted(peak_sig, cmp=value_cmp, reverse=True)

    peaks=[]
    for i in range(min(len(peak_sig), num)):
        peaks.append(peak_sig[i])
    return peaks

def find_common_peaks(key_slider, peak_num=7, sample_list=None):
    if sample_list == None:
        sample_list = [key_slider.keys()]
        
    clcpeaks_list, crcpeaks_list, cdpeaks_list = [], [], []
    for i in range(len(sample_list)):
        av_lcutmap, av_rcutmap, av_dyadmap = None, None, None
        key_list = sample_list[i]
        for j in range(len(key_list)):
            key = key_list[j]
            slider = key_slider[key]
            lcut_map = slider.left_cutmap
            rcut_map = slider.right_cutmap
            dyad_map = slider.dyadmap

            if not av_lcutmap:
                av_lcutmap = lcut_map
            else:
                av_lcutmap = [av_lcutmap[k]+lcut_map[k] for k in range(len(lcut_map))]
            if not av_rcutmap:
                av_rcutmap = rcut_map
            else:
                av_rcutmap = [av_rcutmap[k]+rcut_map[k] for k in range(len(rcut_map))]
            if not av_dyadmap:
                av_dyadmap = dyad_map
            else:
                av_dyadmap = [av_dyadmap[k]+dyad_map[k] for k in range(len(dyad_map))]

        av_lcutmap = [av_lcutmap[l]/len(key_list) for l in range(len(av_lcutmap))]
        av_rcutmap = [av_rcutmap[l]/len(key_list) for l in range(len(av_rcutmap))]
        av_dyadmap = [av_dyadmap[l]/len(key_list) for l in range(len(av_dyadmap))]

        av_lcutpeaks = find_peaks(av_lcutmap, num=peak_num)
        av_rcutpeaks = find_peaks(av_rcutmap, num=peak_num)
        av_dyadpeaks = find_peaks(av_dyadmap, num=peak_num)

        clcpeaks = [peak for peak, sig in av_lcutpeaks]
        crcpeaks = [peak for peak, sig in av_rcutpeaks]
        cdpeaks = [peak for peak, sig in av_dyadpeaks]

        clcpeaks_list.append(clcpeaks)
        crcpeaks_list.append(crcpeaks)
        cdpeaks_list.append(cdpeaks)
        
    if len(sample_list) <= 1:
        return clcpeaks_list[0], crcpeaks_list[0], cdpeaks_list[0]

    return clcpeaks_list, crcpeaks_list, cdpeaks_list

def map_reduction (key_slider, clcpeaks, crcpeaks, cdpeaks):
    rkey_slider = {}
    for key in key_slider:
        slider = key_slider[key]
        lcutmap = slider.left_cutmap
        rcutmap = slider.right_cutmap
        dyadmap = slider.dyadmap
        rlcutmap, rrcutmap, rdyadmap = [], [], []
        assert len(lcutmap) == len(rcutmap) == len(dyadmap)
        for j in range(len(lcutmap)):
            if j in clcpeaks:
                rlcutmap.append(lcutmap[j])
            if j not in clcpeaks:
                rlcutmap.append(0)
            if j in crcpeaks:
                rrcutmap.append(rcutmap[j])
            if j not in crcpeaks:
                rrcutmap.append(0)
            if j in cdpeaks:
                rdyadmap.append(dyadmap[j])
            if j not in cdpeaks:
                rdyadmap.append(0)
        if sum(rdyadmap) < 50:
            continue
        #rlcutmap, rrcutmap, rdyadmap = norm(rlcutmap), norm(rrcutmap), norm(rdyadmap)
        #rdyadmap = norm(rdyadmap)
        assert key not in rkey_slider
        rslider = copy.deepcopy(slider)
        rslider.left_cutmap = rlcutmap
        rslider.right_cutmap = rrcutmap
        rslider.dyadmap = rdyadmap
        rkey_slider[key] = rslider
    return rkey_slider

# get cut/dyad peak list
def get_peaks(key_slider, left_cut_peak_num, right_cut_peak_num, dyad_peak_num, signal=False, sample_list=None):    
    if sample_list == None:
        sample_list = [key_slider.keys()]

    cutpeaks_list, dyadpeaks_list = [], []
    for i in range(len(sample_list)): 
        key_list = sample_list[i]
        cutpeaks, dyadpeaks = {}, {}
        for j in range(len(key_list)):
            key = key_list[j]
            slider = key_slider[key]
            lcut_peak = find_peaks(slider.left_cutmap, num = left_cut_peak_num)
            rcut_peak = find_peaks(slider.right_cutmap, num = right_cut_peak_num)
            dyad_peak = find_peaks(slider.dyadmap, num = dyad_peak_num)

            if not signal:
                lcut_peak = [peak for peak, sig in lcut_peak]
                rcut_peak = [peak for peak, sig in rcut_peak]
                dyad_peak = [peak for peak, sig in dyad_peak]
                
            assert key not in cutpeaks
            assert key not in dyadpeaks
            cutpeaks[key] = {}
            cutpeaks[key]['L'] = lcut_peak
            cutpeaks[key]['R'] = rcut_peak
            dyadpeaks[key] = dyad_peak

        cutpeaks_list.append(cutpeaks)
        dyadpeaks_list.append(dyadpeaks)

    if len(sample_list) <= 1:
        return cutpeaks_list[0], dyadpeaks_list[0]
    return cutpeaks_list, dyadpeaks_list

def Kmeans(dic, cluster_num, sample_list=None, type_targets=[None, []]):
    if sample_list == None:
        sample_list = [dic.keys()]
    key_cdx_list, cdx_key_list = [], []
    type, targets = type_targets
    for u in range(len(sample_list)):
        key_list = sample_list[u]
        X = []
        for key in key_list:
            if type == "slicing":
                temp = []
                for st, ed in sorted(targets):
                    temp += dic[key][st:ed]
                X.append(temp)
            if type == "exclude":
                for k in sorted(targets, reversed=True):
                    del dic[key][k]
                X.append(dic[key])
            elif type == "include":
                X.append([dic[key][k] for k in targets])
            else:
                X.append(dic[key])
        X = []
        for key in key_list:
            X.append(dic[key])
        y = sklearn.cluster.KMeans(init='k-means++', n_init=10, n_clusters = cluster_num, max_iter = 10000).fit_predict(np.asarray(X))
        #y = sklearn.cluster.DBSCAN().fit_predict(np.asarray(X))
        #print y
        #y = sklearn.cluster.SpectralClustering(n_clusters = cluster_num).fit_predict(np.asarray(X))
        #model = sklearn.cluster.AgglomerativeClustering(n_clusters=cluster_num, linkage='ward').fit(np.asarray(X))
        #y = model.labels_
        key_cdx, cdx_key = {}, {}
        for i in range(len(y)):
            key, cdx = key_list[i], y[i]
            assert key not in key_cdx
            key_cdx[key] = cdx
            if cdx not in cdx_key:
                cdx_key[cdx] = []
            cdx_key[cdx].append(key)
        key_cdx_list.append(key_cdx)
        cdx_key_list.append(cdx_key)

    if len(sample_list) <= 1:
        return key_cdx_list[0], cdx_key_list[0]
    return key_cdx_list, cdx_key_list

def PCA(dic, comp_num, sample_list=None, norm_choice=False):
    if sample_list == None:
        sample_list = [dic.keys()]

    key_comp_list, in_comp_list, out_comp_list = [], [], []
    for u in range(len(sample_list)):
        key_list = sample_list[u]
        X = []
        for key in key_list:
            temp = dic[key]
            if norm_choice:
                temp = norm(temp)
            X.append(temp)
            
        pca = sklearn.decomposition.PCA(n_components = comp_num)
        #pca = sklearn.decomposition.IncrementalPCA(n_components = comp_num)
        pca.fit(X)
        Y = pca.transform(X)
        clf = sklearn.covariance.EllipticEnvelope(contamination=0.01)
        clf.fit(Y)
        Y_pred = clf.predict(Y)

        key_comp = {}
        in_comp, out_comp = {}, {}
        for i in range(len(key_list)):
            key, comp, pred = key_list[i], Y[i], Y_pred[i]
            assert key not in key_comp
            key_comp[key] = comp
            if pred > 0:
                in_comp[key] = comp
            else:
                assert pred < 0
                out_comp[key] = comp
        key_comp_list.append(key_comp)
        in_comp_list.append(in_comp)
        out_comp_list.append(out_comp)

        for key in out_comp:
            print color_A(key)
        
        fig = plt.figure()
        A = np.asarray(in_comp.values())
        B = np.asarray(out_comp.values())
        plt.scatter(A[:,0], A[:,1], s = 2, c='k', edgecolor='k')
        sns.kdeplot(A[:,0], A[:,1], shade=False, shade_lowest=False)
        plt.scatter(B[:,0], B[:,1], s = 2, c='r', edgecolor='r')
        plt.xlabel('PC 1')
        plt.ylabel('PC 2')
        plt.savefig('PCA.png')
        plt.close()

        graph.plot_weblogo(out_comp.keys(), note='PCA_outliers')
        
    if len(sample_list) <= 1:
        return key_comp_list[0], in_comp_list[0], out_comp_list[0]
    return key_comp_list, in_comp_list, out_comp_list

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

def KL_div (dist1, dist2):
    coeff = 0.0
    for i in range(len(dist1)):
        if  dist1[i] == 0 or dist2[i] == 0:
            continue
        coeff += dist1[i] * np.log2(float(dist1[i])/dist2[i])
    return coeff

def JS_dist (dist1, dist2):
    assert len(dist1) == len(dist2)
    def JS_div (dist1, dist2):
        def KL_div (dist1, dist2):
            coeff = 0.0
            for i in range(len(dist1)):
                if  dist1[i] == 0 or dist2[i] == 0:
                    continue
                coeff += dist1[i] * np.log2(float(dist1[i])/dist2[i])
            return coeff
        m = [(dist1[i] + dist2[i])*0.5 for i in range(len(dist1))]
        return 0.5*KL_div(dist1, m) + 0.5*KL_div(dist2, m)
    if sum(dist1) != 1:
        dist1 = norm(dist1)
    if sum(dist2) != 1:
        dist2 = norm(dist2)
    return np.sqrt(JS_div(dist1, dist2))

def distance_matrix (key_slider):
    pair_dist = {}
    keys = key_slider.keys()
    for i in range(len(keys)-1):
        for j in range(i+1, len(keys)):
            key1, key2 = keys[i], keys[j]
            dist = JS_dist(key_slider[key1].dyadmap, key_slider[key2].dydmap)
            assert (key1, key2) not in pair_dist
            pair_dist[(key1, key2)] = dist
    return pair_dist        

def compare_map (key_slider1, key_slider2):
    def Bhatta_coeff (dist1, dist2):
        assert len(dist1) == len(dist2)
        coeff = 0.0
        for i in range(len(dist1)):
            coeff += math.sqrt(dist1[i] * dist2[i])
        return coeff

    def JS_div (dist1, dist2):
        assert len(dist1) == len(dist2)
        def KL_div (dist1, dist2):
            coeff = 0.0
            for i in range(len(dist1)):
                if  dist1[i] == 0 or dist2[i] == 0:
                    continue
                coeff += dist1[i] * np.log2(float(dist1[i])/dist2[i])
            return coeff
        m = [(dist1[i] + dist2[i])*0.5 for i in range(len(dist1))]
        return 0.5*KL_div(dist1, m) + 0.5*KL_div(dist2, m)
    
    keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))
    key_coeff = {}
    for key in keys:
        size, loc = key.split('-')
        size, loc = int(size), int(loc)
        if size > 13:
            continue
        dyadmap1 = key_slider1[key].dyadmap
        dyadmap2 = key_slider2[key].dyadmap
        #dyadmap3 = random_map(range(len(dyadmap1)), int(sum(dyadmap1)))
        #dyadmap1, dyadmap2, dyadmap3 = norm(dyadmap1), norm(dyadmap2), norm(dyadmap3)
        dyadmap1, dyadmap2  = norm(dyadmap1), norm(dyadmap2)
        #coeff = Bhatta_coeff(dyadmap1, dyadmap2)
        #limit = np.sqrt(JS_div(dyadmap1, dyadmap3))
        coeff = np.sqrt(JS_div(dyadmap1, dyadmap2))
        #if coeff > limit:
        #    coeff = limit
        #coeff = coeff / limit
        assert key not in key_coeff
        #key_coeff[key] = np.sqrt(1-coeff)
        key_coeff[key] = coeff

    #fig = plt.figure()
    #n, bins, patches = plt.hist(key_coeff.values(), bins=50, range=(0,1), facecolor='green', alpha=0.75)
    #plt.xlabel("Hellinger distance")
    #plt.xlabel("Bhattacharyya coefficient")
    #plt.ylabel("Counts")
    #plt.savefig("t.png")
    #plt.show()

    return key_coeff




    


