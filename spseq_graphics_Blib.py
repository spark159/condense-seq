
from argparse import ArgumentParser, FileType
import matplotlib.pyplot as plt
import numpy as np
import random  
import math

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

def display_graph (norm_fname_list,
                   test_fname_list,
                   ref_fname,
                   count_cutoff):
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
        
    def quicksort(sort_fname, id_seq):
        type_count = {}
        chosen=[]
        for line in open(sort_fname):
            if line.strip():
                read_id, type, seq =line.strip().split()
                if type != "mutant":
                    ref_id = int(type[3:])
                    chosen.append(id_seq[ref_id])
                    type = "valid"
                if type not in type_count:
                    type_count[type] = 0
                type_count[type] +=1
        return type_count, chosen

    def read_sortfile(sort_fname):
        seq_list={}
        for line in open(sort_fname):
            if line.strip():
                #print line
                read_id, type, insert, cutloc, seq =line.strip().split()
                seq_list[read_id]=[type,insert]
        return seq_list

    def sort(seq_list, choice="freeDNA"):
        type_count={}; chosen=[]
        for id in seq_list:
            type, insert=seq_list[id]
            if type not in type_count:
                type_count[type]=0
            type_count[type] +=1
            if type==choice:
                chosen.append(insert)
        return type_count, chosen 

    def seq_counting(seq_list):
        output={}
        for seq in seq_list:
            if seq not in output:
                output[seq] =0
            output[seq] +=1
        return output
    
    def GC_content(seq):
        num=0.0
        for nt in seq:
            if nt in 'GC':
                num+=1
        return (num/float(len(seq)))*100

    def GC_counting(seq_list):
        GC_count={}
        for seq in seq_list:
            GC=GC_content(seq)
            if GC not in GC_count:
                GC_count[GC] = 0
            GC_count[GC] +=1
        return GC_count
    
    def normalize (prob):
        total=0.0
        for key in prob:
            total += prob[key]
        for key in prob:
            prob[key]=prob[key]/float(total)
        return

    def get_hist (data, binnum=1000, prob=True):
        hist={};
        if prob:
            deno=float(len(data))
        else:
            deno=1
        binwidth=float(max(data)-min(data))/binnum
        for value in data:
            bin=int((value-min(data))/binwidth)
            bincenter=min(data)+(bin+0.5)*binwidth
            if bincenter not in hist:
                hist[bincenter]=0
            hist[bincenter]+=1/deno
        return hist
    
    def all_path(N, states='AC'):
        if N==1:
            return list(states)
        output=[]
        for path in all_path(N-1):
            for state in states:
                output.append(path+state)
        return output

    def random_pick (data, N, bias=0):
        output=[]; i=0
        while i < N:
            pick=data[random.randint(0,len(data)-1)]
            GC=GC_content(pick)/100.0
            if random.random() < math.exp(-bias*GC):
                output.append(pick)
                i +=1
        return output

    def nt_freq (seq_list):
        freq=[ {} for i in range(len(seq_list[0])) ]
        N=float(len(seq_list))
        for seq in seq_list:
            for i in range(len(seq)):
                nt=seq[i]
                if nt not in freq[i]:
                    freq[i][nt]=0.0
                freq[i][nt]+=1.0/N
        return freq

    def GC_profile(freq):
        AT=[]; GC=[]
        for pos in freq:
            count1=0.0; count2=0.0
            for nt in pos.keys():
                if nt in 'AT':
                    count1 +=pos[nt]
                else:
                    count2 +=pos[nt]
            AT.append(count1)
            GC.append(count2)
        return AT, GC

    def count_cmp(a, b):
        if a[1] >= b[1]:
            return -1
        else:
            return 1

    def pick_best(seq_count, N=1000):
        seqcount=[[key,value] for key, value in seq_count.items()]
        seqcount=sorted(seqcount, cmp=count_cmp)
        return seqcount[:N]

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
        return diffprod / math.sqrt(xdiff2 * ydiff2)

    id_seq = read_ref(ref_fname)
    
    # combine all norm seq files into one data set
    norm_seqs=[]
    for norm_fname in norm_fname_list:
        type_count, good_seqs=quicksort(norm_fname, id_seq)
        norm_seqs += good_seqs

    # each test seq files organized into separte data set
    test_seqs_list=[]
    for test_fname in test_fname_list:
        type_count, good_seqs=quicksort(test_fname, id_seq)
        test_seqs_list.append(good_seqs)

    """
    # find valid target inserts
    norm_insert_count = seq_counting(norm_seqs)
    test_insert_count = seq_counting(test_seqs)
    
    common_insert = list(set(norm_insert_count.keys()) & set(test_insert_count.keys()))
    target_insert=[]
    for insert in common_insert:
        if norm_insert_count[insert] > count_cutoff and test_insert_count[insert] > count_cutoff:
            target_insert.append(insert)

    temp1=[]; temp2=[]
    for insert in norm_seqs:
        if insert in target_insert:
            temp1.append(insert)
    for insert in test_seqs:
        if insert in target_insert:
            temp2.append(insert)
    norm_seqs=temp1
    test_seqs=temp2
    
    norm_insert_count = seq_counting(norm_seqs)
    test_insert_count = seq_counting(test_seqs)
    """
    
    # make simulated insertion list
    all_insert=id_seq.values()
    N = len(norm_seqs)
    ideal_seqs=random_pick(all_insert,N)
    #bias_seqs=random_pick(all_insert,N,4)
    
    # plot GC content vs. frequency
    fig = plt.figure()
    
    norm_GC_count=GC_counting(norm_seqs); normalize(norm_GC_count)
    plt.plot(norm_GC_count.keys(), norm_GC_count.values(), color='black', marker='D', markersize=5, linestyle= 'None', label='Input')
    idealGC_count=GC_counting(ideal_seqs); normalize(idealGC_count)
    plt.plot(idealGC_count.keys(), idealGC_count.values(),color='green',marker='x',  markersize=5, linestyle= 'None', label='Ideal')
    for i in range(len(test_seqs_list)):
        test_seqs = test_seqs_list[i]
        test_GC_count=GC_counting(test_seqs); normalize(test_GC_count)
        plt.plot(test_GC_count.keys(), test_GC_count.values(), marker="o", markersize=5, linestyle= 'None', label='sp'+str(i+1) )
    #biasGC_count=GC_counting(bias_seqs); normalize(biasGC_count)
    #plt.scatter(biasGC_count.keys(), biasGC_count.values(),color='red',marker='*', label='Biased')
    plt.xlabel('GC contents (%)')
    plt.ylabel('Frequency')
    plt.legend(loc='best', numpoints=1)
    plt.savefig("GC_freq.png")
    plt.close()
    #plt.show()

    # plot seqID vs. coverage
    fig = plt.figure()
    norm_insert_count = seq_counting(norm_seqs)
    plt.plot(range(len(norm_insert_count)), sorted(norm_insert_count.values()),color='black', label='Input')
    plt.axvline(len(norm_insert_count), color='k', linestyle='dotted')
    idealinsert_count=seq_counting(ideal_seqs)
    plt.plot(range(len(idealinsert_count)), sorted(idealinsert_count.values()),color='green', label='Ideal')
    plt.axvline(len(idealinsert_count), color='g', linestyle='dotted')
    test_insert_count_list = []
    for i in range(len(test_seqs_list)):
        test_seqs = test_seqs_list[i]
        test_insert_count = seq_counting(test_seqs)
        test_insert_count_list.append(test_insert_count)
        plt.plot(range(len(test_insert_count)), sorted(test_insert_count.values()), label='sp'+str(i+1))
        plt.axvline(len(test_insert_count), color='b', linestyle='dotted')
    #biasinsert_count=seq_counting(bias_seqs)
    #plt.plot(range(len(biasinsert_count)), sorted(biasinsert_count.values()),color='red', label='Biased')
    #plt.axvline(len(biasinsert_count), color='r', linestyle='dotted')
    plt.ylim([0,500])
    plt.yticks([0,500,50])
    plt.xlabel('Seq ID')
    plt.ylabel('Counts')
    plt.legend(loc='best')
    plt.savefig("coverage.png")
    plt.close()
    #plt.show()



    # GC content per pool
    norm_pool=''
    for insert, count in norm_insert_count.items():
        norm_pool += insert*count
    GCnorm = GC_content(norm_pool) 
    print "norm", str(GCnorm)
    for i in range(len(test_insert_count_list)):
        test_insert_count = test_insert_count_list[i]
        test_pool=''
        for insert, count in test_insert_count.items():
            test_pool += insert*count
        GCtest = GC_content(test_pool)
        print 'sp' + str(i+1), str(GCtest), str(float(GCtest)/GCnorm)

    """
    # find valid target inserts
    common_insert = list(set(norm_insert_count.keys()) & set(test_insert_count.keys()))
    target_insert=[]
    for insert in common_insert:
        if norm_insert_count[insert] > count_cutoff and test_insert_count[insert] > count_cutoff:
            target_insert.append(insert)
    """
    target_insert = all_insert
   
    # GC contents vs norm. counts per seq
    normalize(norm_insert_count)
    for i in range(len(test_insert_count_list)):
        normalize(test_insert_count_list[i])
        #print test_insert_count_list[i]

    for i in range(len(test_insert_count_list)):
        test_insert_count = test_insert_count_list[i]

        insert_ncount={}; insert_GC = {}
        GC_ncount = {}

        X=[]; Y=[]
        X2=[]; Y2=[]
        for insert in target_insert:
            if insert not in norm_insert_count:
                continue
            if insert not in test_insert_count:
                test_count = 0.0
            else:
                test_count = test_insert_count[insert]
            X2.append(norm_insert_count[insert])
            Y2.append(test_count)
            GC = GC_content(insert)
            ncount = float(test_count) / norm_insert_count[insert]
            insert_ncount[insert] = ncount
            insert_GC[insert] = GC
            X.append(GC)
            Y.append(ncount)
            if GC not in GC_ncount:
                GC_ncount[GC] = []
            GC_ncount[GC].append(ncount)

        f = open("GCcorr_sp"+str(i+1) +".txt", 'w')
        corr = get_corr(X,Y)
        print >> f, "GC", corr
        f.close()
            
        fig = plt.figure()
        plt.scatter(X2,Y2)
        plt.plot([0,0.004],[0,0.004], 'k--')
        plt.text(0.004,0.004, 'y=x', fontsize=20)
        plt.xlim([-0.001,0.005])
        plt.ylim([-0.001,0.005])
        plt.xlabel('Norm. Input Counts')
        plt.ylabel('Norm. Data Counts')
        plt.savefig("inputVStest_scatter" + str(i+1) + ".png")
        plt.close()
        #plt.show()

        """
        par=np.polyfit(X,Y,1)
        slope=par[0]; intercept=par[1]
        xl = range(100)
        yl = [slope*x + intercept  for x in range(100)]
        """
        fig = plt.figure()
        plt.scatter(X,Y)
        #plt.plot(xl,yl,'k--')
        #plt.text(xl[-1],yl[-1], 'Slope:' + str(slope), fontsize=10)
        plt.xlabel('GC contents (%)')
        plt.ylabel('Norm. counts')
        plt.savefig("GCVScount_scatter" + str(i+1) + ".png")
        plt.close()
        #plt.show()

        x, y, z = [], [], []
        for GC in GC_ncount:
            x.append(GC)
            y.append(np.mean(GC_ncount[GC]))
            z.append(np.std(GC_ncount[GC]))
        fig = plt.figure()
        plt.plot(x,y,'.')
        plt.errorbar(x,y,yerr=z,fmt='o')
        plt.xlabel('GC contents')
        plt.ylabel('Norm. counts')
        #$plt.show()
        plt.savefig('GCVScount.png'+ str(i+1) + ".png")
        plt.close()
                                        
        f = open("ncount_sp"+str(i+1)+'.txt', 'w')

        for i in range(len(insert_ncount)):
            insert = insert_ncount.keys()[i]
            ncount = insert_ncount[insert]
            print >> f, '%d\t%s\t%f' % (i, insert, ncount)

        f.close()

        
    """
    plt.plot(range(len(insert_ncount)), insert_ncount.values())
    plt.show()
    
    
    best_insert_count=pick_best(seq_counting(good_seqs), N=2000)
    tmp=[]
    for key, item in best_insert_count:
        tmp += [key]*item
    good_seqs=tmp
    
    
    # pie chart for classify
    types=type_count.keys(); counts=type_count.values()
    explode=[]
    for key in types:
        if key=='freeDNA':
            explode.append(0.05)
        else:
            explode.append(0)
    colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue', 'lightpink', 'lime']
    plt.pie(counts,colors=colors, labels=types, shadow=True, startangle=90, autopct='%1.1f%%', explode=explode)
    plt.axis('equal')
    plt.show()
    
    # GC/AT contents profile
    freq=nt_freq(good_seqs)
    ATprofile, GCprofile=GC_profile(freq)
    ideal,_=GC_profile(nt_freq(ideal_seqs))
    bias_AT,bias_GC=GC_profile(nt_freq(bias_seqs))
    plt.plot(range(len(ATprofile)),ATprofile,'r', label='AT frequency')
    plt.plot(range(len(GCprofile)),GCprofile,'b', label='GC frequency')   
    plt.plot(range(len(ideal)),ideal,'k--')
    plt.plot(range(len(bias_AT)),bias_AT,'r:')
    plt.plot(range(len(bias_GC)),bias_GC,'b:')
    plt.xlabel('Position (bp)')
    plt.ylabel('frequency')
    plt.legend(loc='upper right')
    plt.show()

    # print 100 best abundance sequences
    insertcount=pick_best(insert_count)
    for i in range(100):
        print "%s\t%d" % (insertcount[i][0], insertcount[i][1])
    
    
    # read coverage histogram
    hist=get_hist(insert_count.values())
    X=[];Y=[]
    for x, y in sorted(hist.items()):
        X.append(x); Y.append(y)
    plt.plot(X,Y, label='data')
    hist2=get_hist(idealinsert_count.values())
    X2=[];Y2=[]
    for x2, y2 in sorted(hist2.items()):
        X2.append(x2); Y2.append(y2)
    plt.plot(X2,Y2, label='ideal')
    plt.xlabel('Read Coverage')
    plt.ylabel('Probability')
    plt.legend(loc='best')
    plt.show()
    """
    
if __name__ == '__main__':
    parser = ArgumentParser(description='graphical analysis of normalized representation')
    parser.add_argument('-n',
                        dest='norm_fname_list',
                        nargs = '+',
                        type=str,
                        help='normalization filenames')
    parser.add_argument('-f',
                        dest='test_fname_list',
                        nargs = '+',
                        type=str,
                        help='test filenames')
    parser.add_argument('-x',
                        dest = 'ref_fname',
                        type=str,
                        help='input filename')
    parser.add_argument('--count--cutoff',
                        dest='count_cutoff',
                        type=int,
                        default=50,
                        help='count cutoff')
    args = parser.parse_args()
                        
    display_graph (args.norm_fname_list,
                   args.test_fname_list,
                   args.ref_fname,
                   args.count_cutoff)
