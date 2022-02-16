import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import LinModel
import matplotlib.cm as cm

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

def get_ATGC_sig (freq):
    freq_MM1 = freq['MM1']
    length = len(freq_MM1)
    nts = LinModel.all_path(2, 'ATCG')
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

def get_nt_sig (freq):
    freq_MM1 = freq['MM0']
    length = len(freq_MM0)
    nts = 'ATCG'
    A_sig, T_sig = np.zeros(length), np.zeros(length)
    C_sig, G_sig = np.zeros(length), np.zeros(length)
    for nt in nts:
        row = [ freq_MM0[i][nt] for i in range(length)]
        if nt ==' A':
            A_sig += np.asarray(row)
    
            AT_sig += np.asarray(row)
        if nt in ['GG', 'GC', 'CG', 'CC']:
            GC_sig += np.asarray(row)
    for i in range(length):
        AT_sig[i] = float(AT_sig[i]) / sum(freq_MM1[i].values())
        GC_sig[i] = float(GC_sig[i]) / sum(freq_MM1[i].values())
    return AT_sig, GC_sig

def get_zscore(p_data):
    linear_data = []
    p_num = []
    for data in p_data:
        p_num.append(len(data))
        for value in data:
            linear_data.append(value)
    linear_data = stats.zscore(linear_data)
    new_data = []
    st = 0
    for num in p_num:
        new_data.append(linear_data[st:st+num])
        st += num
    assert st == len(linear_data)
    return new_data


id_seq = read_ref('Blib.ref')
seq_MGW, seq_HelT, seq_ProT, seq_Roll = read_DNAshape('php6POcc7', id_seq)
_, mean_list, _ = read_titration("AakashDNAlib_spermine.csv")
fidx_to_tnum = [None, 4, 7, 8, [9,10], 11, 12]

fname_list = ['ncount_sp' + str(i+1) + '.txt' for i in range(7)]
Models = []
for k in range(len(fname_list)):
    if k not in [3, 5]:
        continue
    fname = fname_list[k]
    seq_ncount = read_ncount(fname)

    tnum = fidx_to_tnum[k]
    if tnum == None:
        mean_fract = 1.0
    elif type(tnum) == list:
        mean_fract = np.mean([mean_list[tnum[0]], mean_list[tnum[1]]])
    else:
        mean_fract = mean_list[tnum]

    seq_list, score_list = [], []
    count_list = []
    for seq in seq_ncount:
        if seq_ncount[seq] <= 0:
            continue
        score = -np.log2(seq_ncount[seq]) - np.log2(mean_fract)
        seq = seq[23:23+101]
        seq_list.append(seq)
        score_list.append(score)
        count_list.append(1)

        ## add reverse complement
        #seq_list.append(rev_comp(seq))
        #score_list.append(score)
        #count_list.append(1)


    M = LinModel.SeqLinearModel(seq_list, score_list, count_list)

    
    # MM 0-th order 
    M.train(MM_orders=[0], Kmer_k_b=None, PolyA_b=None, GC_b=None, Harmonic=None, sym=True, graph=False)
    fig = plt.figure(figsize=(2.6, 2))
    nts = 'ATCG'
    color_list = ['tab:red', 'tab:orange', 'tab:green', 'tab:blue']
    for i in range(len(nts)):
        nt = nts[i]
        coeff = M.coeff['MM0']
        length = len(coeff)
        sig = [coeff[u][nt] for u in range(length)]
        plt.plot(sig, color=color_list[i], lw=1.5, alpha=1, label=nt)
        plt.gca().tick_params('y', labelsize=8)
        plt.ylabel("Coefficient", fontsize=8)
        plt.xlabel("Position (bp)", fontsize=8)
        line_list = [101/2 - i*20 for i in range(3)] + [101/2 + i*20 for i in range(1, 3)]
        for line in line_list:
            plt.axvline(x=line, color='k', linestyle='--', alpha=0.25)
        plt.xticks([101/2 + 10*i for i in range(-5, 6, 2)], [str(10*i) for i in range(-5, 6, 2)], fontsize=5)
        plt.gca().tick_params('x', labelsize=6)
        plt.gca().tick_params('y', labelsize=6)

    plt.legend(loc='upper right', fontsize=8)
    plt.savefig('Blib_ntCoeff_' + str(k) + '.svg', format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

    # MM n-th order
    order = 1
    M.train(MM_orders=[order], Kmer_k_b=None, PolyA_b=None, GC_b=None, Harmonic=None, sym=True, graph=False)
    dins = [din for _, din in sorted([(GC_content(din), din) for din in LinModel.all_path(order+1, 'ATCG')])]
    img = []
    for din in dins:
        coeff = M.coeff['MM' + str(order)]
        length = len(coeff)
        sig = [coeff[u][din] for u in range(length)]
        img.append(sig)

    limit = min(abs(np.min(np.asarray(img))), abs(np.max(np.asarray(img))))

    fig = plt.figure(figsize=(2.6, 2))
    #fig = plt.figure()
    im = plt.imshow(img, aspect='auto', cmap='seismic', vmin=-limit, vmax=limit)
    plt.yticks(range(len(dins)), dins, fontsize=8)
    plt.xlabel("Position (bp)", fontsize=8)
    line_list = [101/2 - i*20 for i in range(3)] + [101/2 + i*20 for i in range(1, 3)]
    for line in line_list:
        plt.axvline(x=line, color='k', linestyle='--', alpha=0.25)
    plt.xticks([101/2 + 10*i for i in range(-5, 6, 2)], [str(10*i) for i in range(-5, 6, 2)], fontsize=5)
    plt.gca().tick_params('x', labelsize=6)
    #plt.gca().tick_params('y', labelsize=8)
    #cbar = plt.colorbar(fraction=0.02)
    #cbar.set_ticks([-limit,0, limit])
    #cbar.set_ticklabels(['<' + str(round(-limit,1)),str(0), str(round(limit,1)) + '>'])
    #cbar.ax.tick_params(labelsize=6)
    #plt.legend(loc='upper right', fontsize=8)
    plt.savefig('Blib_dinCoeff_' + str(k) + '.svg', format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

    fig = plt.figure(figsize=(1.1,1))
    plt.subplot(1,2,1)
    cbar = plt.colorbar(im, cax=plt.gca())
    #cbar.ax.set_yticklabels([str(vmin), str(vmax)], fontsize=8)
    cbar.set_ticks([-limit,0, limit])
    cbar.set_ticklabels(['<' + str(round(-limit,1)),str(0), str(round(limit,1)) + '>'])
    cbar.ax.tick_params(labelsize=5)
    cbar.ax.set_ylabel('Coefficient', rotation=-90, va="bottom", fontsize=6)
    plt.tight_layout()
    plt.savefig('dincoeff_' + str(k)  + '_cbar.svg', format='svg', bbox_inches='tight')
    plt.close()


    # kmer
    kmer = 2
    M.train(MM_orders=[], Kmer_k_b=[kmer,1], PolyA_b=None, GC_b=None, Harmonic=None, sym=True, graph=False)
    kmers = [kmer for _, kmer in sorted([(GC_content(kmer), kmer) for kmer in LinModel.all_path(kmer, 'ATCG')])]

    fig = plt.figure()
    X, Y = [], []
    for i in range(len(kmers)):
        X.append(i)
        kmer = kmers[i]
        Y.append(M.coeff['Kmer0'][kmer])
    plt.plot(X, Y, '.')
    plt.xticks(X, kmers, rotation=45, fontsize=5)
    plt.show()
    plt.close()
    

    """
    group_freq = M.spectrum(MM_orders=[0, 1], Kmer_k_b=[4,1], gnum=7, PolyA_b=None, GC_b=None, Harmonic=None)

    color_list = np.linspace(0.2, 1, num=len(group_freq))

    # plot ATCG freq along position

    nts = 'ATCG'
    cmap_list = ['Reds', 'Purples','Greens', 'Blues']
    for i in range(len(nts)):
        nt = nts[i]
        cmap = cm.get_cmap(cmap_list[i])
        #fig = plt.figure(figsize=(2.6, 2))
        fig = plt.figure(figsize=(1.4, 1.2))
        for j in range(len(group_freq)):
            freq = group_freq[j]
            freq_MM0 = freq['MM0']
            length = len(freq_MM0)
            sig = [(freq_MM0[u][nt])/sum(freq_MM0[u].values()) for u in range(length)]
            plt.plot(sig, color=cmap(color_list[j]), lw=1.5, alpha=1)
        plt.title(nt, fontsize=10)
        #plt.ylabel("Frequency", fontsize=8)
        #plt.xlabel("Super Helical Location", fontsize=8)
        line_list = [101/2 - i*20 for i in range(3)] + [101/2 + i*20 for i in range(1, 3)]
        for line in line_list:
            plt.axvline(x=line, color='k', linestyle='--', alpha=0.25)
        plt.xticks([101/2 + 10*i for i in range(-5, 6, 2)], [str(10*i) for i in range(-5, 6, 2)], fontsize=5)
        plt.gca().tick_params('x', labelsize=6)
        plt.gca().tick_params('y', labelsize=6)

        plt.savefig('Blib_ntfreq_' + nt + str(k) + '.svg', format='svg', bbox_inches='tight')
        #plt.show()
        plt.close()

    # plot A/T C/G freq along position
    
    nts = 'AC'
    cmap_list = ['OrRd', 'GnBu']
    tick_colors = ['r', 'b']
    for i in range(len(nts)):
        nt = nts[i]
        cmap = cm.get_cmap(cmap_list[i])
        fig = plt.figure(figsize=(2.6, 2))
        for j in range(len(group_freq)):
            freq = group_freq[j]
            freq_MM0 = freq['MM0']
            length = len(freq_MM0)
            sig = [(freq_MM0[u][nt]+freq_MM0[u][rev_comp(nt)])/sum(freq_MM0[u].values()) for u in range(length)]
            plt.plot(sig, color=cmap(color_list[j]), lw=1.5, alpha=1)
        #plt.title(nt)
        plt.gca().tick_params('y', colors=tick_colors[i], labelsize=8)
        plt.ylabel("%s/%s frequency" % (nt, rev_comp(nt)), fontsize=8, color=tick_colors[i])
        if i % 2 != 0:
            plt.gca().yaxis.tick_right()
            plt.gca().yaxis.set_label_position('right')
        plt.xlabel("Super Helical Location", fontsize=8)
        line_list = [101/2 - i*20 for i in range(3)] + [101/2 + i*20 for i in range(1, 3)]
        for line in line_list:
            plt.axvline(x=line, color='k', linestyle='--', alpha=0.25)
        plt.xticks([101/2 + 10*i for i in range(-5, 6, 2)], [str(10*i) for i in range(-5, 6, 2)], fontsize=5)
        plt.savefig('Blib_ntfreq_' + nt + rev_comp(nt) + str(k) + '.svg', format='svg', bbox_inches='tight')
        #plt.show()
        plt.close()


    # plot all din freq along position
    order = 1
    dins = [din for _, din in sorted([(GC_content(din), din) for din in LinModel.all_path(order+1, 'ATCG')])]
    cmap_list = [ 'Greys' ] * len(dins)
    for i in range(len(dins)):
        din = dins[i]
        cmap = cm.get_cmap(cmap_list[i])
        #fig = plt.figure(figsize=(2.6, 2))
        fig = plt.figure(figsize=(1.4, 1.2))
        for j in range(len(group_freq)):
            freq = group_freq[j]
            freq_MM1 = freq['MM1']
            length = len(freq_MM1)
            sig = [(freq_MM1[u][din])/sum(freq_MM1[u].values()) for u in range(length)]
            plt.plot(sig, color=cmap(color_list[j]), lw=1.5, alpha=1)
        plt.title(din, fontsize=10)
        #plt.gca().tick_params('y', labelsize=8)
        #plt.ylabel("%s frequency" % (din), fontsize=8)
        #plt.xlabel("Super Helical Location", fontsize=8)
        line_list = [101/2 - i*20 for i in range(3)] + [101/2 + i*20 for i in range(1, 3)]
        for line in line_list:
            plt.axvline(x=line, color='k', linestyle='--', alpha=0.25)
        plt.xticks([101/2 + 10*i for i in range(-5, 6, 2)], [str(10*i) for i in range(-5, 6, 2)], fontsize=5)
        plt.gca().tick_params('x', labelsize=6)
        plt.gca().tick_params('y', labelsize=6)
        plt.savefig('Blib_dingroup_' + din + str(k) + '.svg', format='svg', bbox_inches='tight')
        #plt.show()
        plt.close()

    # plot AT-rich, GC-rich din freq along position
    
    cmap1, cmap2 = cm.get_cmap("OrRd"), cm.get_cmap("GnBu")
    cmap = cm.get_cmap("rainbow")

    fig = plt.figure(figsize=(2.6, 2))
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    for i in range(len(group_freq)):
        freq = group_freq[i]
        AT_sig, GC_sig = get_ATGC_sig(freq)
        #ax1.plot(AT_sig, color=cmap1(color_list[i]), lw=1.5, alpha=1)
        ax2.plot(GC_sig, color=cmap2(color_list[i]), lw=1.5, alpha=1)
    #ax1.set_ylabel('AA/AT/TA/TT freqeuncy', color='r', fontsize=8)
    ax2.set_ylabel('CC/CG/GC/GG freqeuncy', color='b', fontsize=8)
    ax1.set_xlabel("Super Helical Location", fontsize=8)

    #ax1.tick_params('y', colors='r', labelsize=8)
    ax2.tick_params('y', colors='b', labelsize=8)

    ax1.set_yticks([])
    ax1.set_yticklabels([])

    #ax2.set_yticks([])
    #ax2.set_yticklabels([])
    
    line_list = [101/2 - i*20 for i in range(3)] + [101/2 + i*20 for i in range(1, 3)]
    for line in line_list:
        plt.axvline(x=line, color='k', linestyle='--', alpha=0.25)
    plt.xticks([101/2 + 10*i for i in range(-5, 6, 2)], [str(10*i) for i in range(-5, 6, 2)], fontsize=5)
    #plt.ylabel("Relative frequency")
    #plt.legend()
    #plt.ylim([0.22, 0.28])
    #plt.savefig('ATGCperiod_' + "slide" + '.png')
    plt.savefig('Blib_ATGCperiod_' + "slide" + str(k) + '.svg', format='svg', bbox_inches='tight')
    plt.title("Dinucleotide periodicity", fontsize=8)
    #plt.show()
    plt.close()
    
    """

    
    

fig = plt.figure(figsize=(1,1))
plt.subplot(1,2,1)
img = [ [color, np.nan] for color in color_list ]
plt.imshow(img, cmap='OrRd', aspect=3)
img = [ [np.nan, color] for color in color_list ]
plt.imshow(img, cmap='GnBu', aspect=3)
plt.yticks(range(len(group_freq)), range(1, len(group_freq)+1), fontsize=8)
ax = plt.gca()
ax.yaxis.tick_right()
#plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
plt.xticks([])
#plt.savefig("cbar.png", bbox_inches='tight')
plt.savefig('cbar.svg', format='svg', bbox_inches='tight')
#plt.show()
plt.close()
