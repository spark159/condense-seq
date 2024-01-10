import matplotlib.pyplot as plt
import numpy as np
import math
from sklearn import linear_model
from mpl_toolkits.axes_grid import AxesGrid
import seaborn as sns

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

def read_ncount(fname):
    seq_ncount = {}
    for line in open(fname):
        if line.strip():
            id, seq, ncount = line.strip().split()
            assert id not in seq_ncount
            seq_ncount[seq] = float(ncount)
    return seq_ncount

def read_DNAshape(fname, id_seq):
    names = ['MGW', 'HelT', 'ProT', 'Roll']
    dic_list = [{} for i in range(len(names))]
    for i in range(len(names)):
        data = []
        for line in open(fname+"."+names[i]):
            line = line.strip()
            if line.startswith('>'):
                if data:
                    assert seq not in dic_list[i]
                    dic_list[i][seq] = data
                id = int(line[1:].strip().split('_')[1])
                seq = id_seq[id]
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
        assert seq not in dic_list[i]
        dic_list[i][seq] = data
    return dic_list

def seqTox (seq, st=None, en=None):
    if st == None and en == None:
        st, en = 0, len(seq)
    x = []
    for nt in seq[st:en]:
        temp = [0 for i in range(4)]
        if nt == 'A':
            temp[0] = 1
        elif nt == 'T':
            temp[1] = 1
        elif nt == 'G':
            temp[2] = 1
        elif nt == 'C':
            temp[3] = 1
        x += temp
    return x

def shapeTox (shape_profile):
    x = []
    for shape in shape_profile:
        for value in shape:
            if type(value) == str:
                continue
            x.append(value)
    return x
    
def removeNA (x):
    nx = []
    for e in x:
        if type(e) != str:
            nx.append(e)
    return nx
    
id_seq = read_ref('Blib.ref')
seq_MGW, seq_HelT, seq_ProT, seq_Roll = read_DNAshape('php6POcc7', id_seq)

fname_list = ['ncount_sp' + str(i+1) + '.txt' for i in range(7)]

f = open("seq_linearfit_error", 'w')

for k in range(len(fname_list)):
    fname = fname_list[k]
    seq_ncount = read_ncount(fname)

    # linear regreession seq to ncount
    X, Y = [], []
    for seq in seq_ncount:
        x = seqTox (seq, st=23, en=23+101)
        y = seq_ncount[seq]
        X.append(x)
        Y.append(y)
    #reg = linear_model.LinearRegression()
    reg = linear_model.Ridge(alpha=0.5)
    #reg = linear_model.RANSACRegressor(linear_model.LinearRegression())
    reg.fit(X, Y)
    print >> f, "sp%d Mean squared error: %.2f" % (k+1, np.mean((reg.predict(X) - np.asarray(Y)) ** 2))
    print >> f, "sp%d variance score: %.2f" % (k+1, reg.score(X,Y))

    pred_Y = reg.predict(X)
    fig = plt.figure()
    plt.scatter(pred_Y, Y, s=1, color='k')
    plt.xlabel('Prediction')
    plt.ylabel('Experiment')
    plt.savefig("NtpredVSexp_sp" + str(k+1) + '.png')
    plt.close()
    

    A_co, T_co, G_co, C_co = [0]*23, [0]*23, [0]*23, [0]*23
    for i in range(len(reg.coef_)):
        #coef = reg.estimator_.coef_
        coef = reg.coef_[i]
        mode = i % 4
        if mode  == 0:
            A_co.append(coef)
        elif mode == 1:
            T_co.append(coef)
        elif mode == 2:
            G_co.append(coef)
        elif mode == 3:
            C_co.append(coef)
    #print len(A_co), len(T_co), len(C_co)
    #print A_co
    A_co, T_co, G_co, C_co = A_co+[0]*23, T_co+[0]*23, G_co+[0]*23, C_co+[0]*23
    co_list = [A_co, T_co, G_co, C_co]

    fig = plt.figure()
    label_list = ['A','T','G','C']
    for j in range(len(co_list)):
        co = co_list[j]
        plt.plot(range(1,len(co)+1), co, label=label_list[j])
    plt.xlabel('Position (bp)')
    plt.ylabel('Coefficient')
    plt.ylim([-0.15,0.15])
    plt.legend(loc='best')
    plt.savefig("sp" + str(k+1) + "_linear_fit_signal.png", bbox_inches='tight')
    plt.close()
    
    img_list = []
    for i in range(4):
        img = []
        temp = []
        for coeff in co_list[i]:
            temp += [coeff] * 2
        for i in range(40):
            img.append(temp)
        img_list.append(img)
    
    fig = plt.figure()
    grid = AxesGrid(fig,
                    111,
                    nrows_ncols=(4, 1),
                    axes_pad=0.25,
                    label_mode="L",
                    cbar_location="right",
                    cbar_mode="single",
                    cbar_size="1.5%")

    labels = ['A','T','G','C']
    for u in range(4):
        ax = grid[u]
        nt = labels[u]
        img = img_list[u]
        ax.set_yticklabels(['0', nt])
        ax.set_xticklabels(['0', '0', '25', '50', '75', '100', '125'])
        im = ax.imshow(img, cmap="seismic", vmin=-0.2, vmax=0.2)

    grid.cbar_axes[0].colorbar(im)

    for cax in grid.cbar_axes:
        cax.toggle_label(True)

    plt.savefig("sp" + str(k+1) + "_linear_fit.png", bbox_inches='tight')
    plt.close()
f.close()

f = open("shape_linearfit_error", 'w')

for k in range(len(fname_list)):
    fname = fname_list[k]
    seq_ncount = read_ncount(fname)

    # linear regression shape to ncount
    X, Y = [], []
    for seq in seq_ncount:
        x1, x2, x3, x4 = seq_MGW[seq][23:23+101], seq_HelT[seq][23:23+100], seq_ProT[seq][23:23+101], seq_Roll[seq][23:23+100]
        x = shapeTox([x1,x2,x3,x4])
        y = seq_ncount[seq]
        X.append(x)
        Y.append(y)
    #reg = linear_model.LinearRegression()
    reg = linear_model.Ridge(alpha=0.5)
    reg.fit(X, Y)
    print >> f, "sp%d Mean squared error: %.2f" % (k+1, np.mean((reg.predict(X) - np.asarray(Y)) ** 2))
    print >> f, "sp%d variance score: %.2f" % (k+1, reg.score(X,Y))

    pred_Y = reg.predict(X)
    fig = plt.figure()
    plt.scatter(pred_Y, Y, s=1, color='k')
    plt.xlabel('Prediction')
    plt.ylabel('Experiment')
    plt.savefig("ShapepredVSexp_sp" + str(k+1) + '.png')
    plt.close()

    #MGW_co, HelT_co, ProT_co, Roll_co = [0,0], [0], [0,0], [0]
    MGW_co, HelT_co, ProT_co, Roll_co = [0]*23, [0]*24, [0]*23, [0]*24
    for i in range(len(reg.coef_)):
        coef = reg.coef_[i]
        if i < 101:
            MGW_co.append(coef)
        elif i >= 101 and i < 101+100:
            HelT_co.append(coef)
        elif i >= 101+100 and i < 101+100+101:
            ProT_co.append(coef)
        elif i >= 101+100+101 and i < 101+100+101+100:
            Roll_co.append(coef)
    MGW_co, HelT_co, ProT_co, Roll_co = MGW_co+[0]*23, HelT_co+[0]*23, ProT_co+[0]*23, Roll_co+[0]*23
    #MGW_co += [0,0]
    #HelT_co += [0,0]
    #ProT_co += [0,0]
    #Roll_co += [0,0]
    #print len(MGW_co), len(HelT_co), len(ProT_co), len(Roll_co)
    #print A_co
    co_list = [MGW_co, HelT_co, ProT_co, Roll_co]

    fig = plt.figure()
    label_list = ['MGW','HelT','ProT','Roll']
    for j in range(len(co_list)):
        co = co_list[j]
        plt.plot(range(1,len(co)+1), co, label=label_list[j])
    plt.xlabel('Position (bp)')
    plt.ylabel('Coefficient')
    plt.ylim([-0.15,0.15])
    plt.legend(loc='best')
    plt.savefig("sp" + str(k+1) + "_shape_linear_fit_signal.png", bbox_inches='tight')
    plt.close()


    img_list = []
    for i in range(4):
        temp = []
        img = []
        for coeff in co_list[i]:
            temp += [coeff] * 2
        for i in range(40):
            img.append(temp)
        img_list.append(img)

    fig = plt.figure()
    grid = AxesGrid(fig,
                    111,
                    nrows_ncols=(4, 1),
                    axes_pad=0.25,
                    label_mode="L",
                    cbar_location="right",
                    cbar_mode="single",
                    cbar_size="1.5%")

    labels = ['MGW','HelT','ProT','Roll']
    for u in range(4):
        ax = grid[u]
        nt = labels[u]
        img = img_list[u]
        ax.set_yticklabels(['0', nt])
        ax.set_xticklabels(['0', '0', '25', '50', '75', '100', '125'])
        im = ax.imshow(img, cmap="seismic", vmin=-0.2, vmax=0.2)

    grid.cbar_axes[0].colorbar(im)

    for cax in grid.cbar_axes:
        cax.toggle_label(True)

    plt.savefig("sp" + str(k+1) + "_shape_linear_fit.png", bbox_inches='tight')
    plt.close()
f.close()
