import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import math

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

def read_titration (fname):
    conc_list = []
    mean_list = []
    std_list = []
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        if not line:
            continue
        cols = line.strip().split()
        conc, mean, std = float(cols[0]), float(cols[1]), float(cols[2])
        conc_list.append(conc)
        mean_list.append(mean)
        std_list.append(std)
    return conc_list, mean_list, std_list

def GC_content(seq):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    return (num/float(len(seq)))*100

def average(l):
    sum = 0.0
    count =0
    for e in l:
        if type(e) != str:
            count +=1
            sum += e
    return sum/count

id_seq = read_ref('Blib.ref')
seq_MGW, seq_HelT, seq_ProT, seq_Roll = read_DNAshape('php6POcc7', id_seq)
_, mean_list, _ = read_titration("AakashDNAlib_spermine.csv")
fidx_to_tnum = [None, 4, 7, 8, [9,10], 11, 12]


fname_list = ['ncount_sp' + str(i+1) + '.txt' for i in range(7)]
color_list = ['tab:red', 'tab:green', 'tab:blue', 'tab:purple', 'tab:brown']

for i in range(len(fname_list)):
    fname = fname_list[i]
    seq_ncount = read_ncount(fname)
    X0 = []
    X1, X2, X3, X4, Y = [], [], [], [], []
    tnum = fidx_to_tnum[i]
    if tnum == None:
        mean_fract = 1.0
    elif type(tnum) == list:
        mean_fract = np.mean([mean_list[tnum[0]], mean_list[tnum[1]]])
    else:
        mean_fract = mean_list[tnum]
    
    for seq in seq_ncount:
        if seq_ncount[seq] <= 0:
            continue
        X0.append(100.0 - GC_content(seq))
        X1.append(average(seq_MGW[seq]))
        X2.append(average(seq_HelT[seq]))
        X3.append(average(seq_ProT[seq]))
        X4.append(average(seq_Roll[seq]))
        #Y.append(seq_ncount[seq])
        Y.append(-np.log2(seq_ncount[seq])-np.log2(mean_fract))

    names = ['AT content (%)', 'MGW', 'HelT', 'ProT', 'Roll']
    Xlist = [X0, X1,X2,X3,X4]
    f = open("corr_sp" + str(i+1) + ".txt", 'w')
    for j in range(len(Xlist)):
        X = Xlist[j]
        #fig = plt.figure(figsize=(2, 1.6))
        fig = plt.figure(figsize=(1.4, 1.2))
        #plt.title('sp%d' % i, fontsize=10)
        plt.plot(X, Y, '.', color=color_list[j], markersize=1, alpha=0.2)
        #plt.plot(X, Y, 'k.', markersize=1, alpha=0.25)
        plt.xlabel(names[j], fontsize=8)
        #plt.ylabel('Norm. counts')
        #plt.ylabel('Condensability (A.U.)', fontsize=6)
        #plt.savefig(names[j] + 'VSncount_sp' +str(i+1) + ".png")
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=5)
        ax.tick_params(axis='both', which='minor', labelsize=5)
        plt.savefig(names[j] + 'VScondensability_sp' +str(i) + ".png", dpi=500, bbox_inches='tight')
        plt.close()
        corr = scipy.stats.spearmanr(X,Y)[0]
        print >> f, '%s\t%f' % (names[j], corr)
    f.close()


# plot the bar graph of correlation
def read_corr (fname):
    name_corr = {}
    for line in open(fname):
        name, corr = line.strip().split('\t')
        if name.startswith('AT'):
            name = 'AT content'
        corr = float(corr)
        name_corr[name] = corr
    return name_corr

fname_list = ['corr_sp' + str(i+1) + '.txt' for i in range(1, 7)]
names = ['AT content', 'MGW', 'HelT', 'ProT', 'Roll']
color_list = ['tab:red', 'tab:green', 'tab:blue', 'tab:purple', 'tab:brown']
fig, axes = plt.subplots(nrows=1, ncols=len(fname_list), figsize=(6.6, 2.6), sharey=True)
for i in range(len(fname_list)):
    fname = fname_list[i]
    name_corr = read_corr(fname)
    corrs = [name_corr[name] for name in names]
    bars = axes[i].bar(range(len(corrs)), corrs, color=color_list[:len(names)], width=0.7)
    axes[i].set_xticks(range(len(corrs)))
    axes[i].set_xticklabels(names, fontsize=7, rotation=45, ha="right", rotation_mode="anchor")
    axes[i].set_title("sp%d" % (i+1), fontsize=10)
    if i <= 0:
        axes[i].set_ylabel("Spearman correlation", fontsize=8)
    axes[i].axhline(y=0, color='k', linestyle='-', lw=0.5)
plt.savefig("Blib_bar.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()
    
