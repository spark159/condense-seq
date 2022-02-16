from DNAextinct import get_fullEC
from argparse import ArgumentParser, FileType
import matplotlib.pyplot as plt
import random  
import math

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

def display_graph (sort_fname,
                   ref_fname):

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

    id_seq = read_ref(ref_fname)
    type_count, good_seqs = quicksort(sort_fname, id_seq)

    print type_count

    """
    best_insert_count=pick_best(seq_counting(good_seqs), N=2000)
    tmp=[]
    for key, item in best_insert_count:
        tmp += [key]*item
    good_seqs=tmp
    """
    
    # pie chart for classify
    types=type_count.keys(); counts=type_count.values()
    explode=[]
    for key in types:
        if key=='valid':
            explode.append(0.05)
        else:
            explode.append(0)
    colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue', 'lightpink', 'lime']
    plt.pie(counts,colors=colors, labels=types, shadow=True, startangle=90, autopct='%1.1f%%', explode=explode)
    plt.axis('equal')
    plt.show()

    # make simulated insertion list
    #all_insert=all_path(16)
    all_insert = id_seq.values()
    assert len(good_seqs)==type_count['valid']
    N=type_count['valid']
    ideal_seqs=random_pick(all_insert,N)
    #bias_seqs=random_pick(all_insert,N,4)
    
    # plot GC content vs. frequency
    GC_count=GC_counting(good_seqs); normalize(GC_count)
    plt.scatter(GC_count.keys(), GC_count.values(), label='Data')
    idealGC_count=GC_counting(ideal_seqs); normalize(idealGC_count)
    plt.scatter(idealGC_count.keys(), idealGC_count.values(),color='green',marker='x', label='Ideal')
    #biasGC_count=GC_counting(bias_seqs); normalize(biasGC_count)
    #plt.scatter(biasGC_count.keys(), biasGC_count.values(),color='red',marker='*', label='Biased')
    plt.xlabel('GC contents (%)')
    plt.ylabel('Frequency')
    plt.legend(loc='best')
    plt.show()

    # plot seqID vs. coverage
    insert_count=seq_counting(good_seqs)
    plt.plot(range(len(insert_count)), sorted(insert_count.values()),color='blue', label='Data')
    plt.axvline(len(insert_count), color='b', linestyle='dotted')
    idealinsert_count=seq_counting(ideal_seqs)
    plt.plot(range(len(idealinsert_count)), sorted(idealinsert_count.values()),color='green', label='Ideal')
    plt.axvline(len(idealinsert_count), color='g', linestyle='dotted')
    #biasinsert_count=seq_counting(bias_seqs)
    #plt.plot(range(len(biasinsert_count)), sorted(biasinsert_count.values()),color='red', label='Biased')
    #plt.axvline(len(biasinsert_count), color='r', linestyle='dotted')
    plt.ylim([0,500])
    plt.yticks(range(0,500,50))
    plt.xlabel('Seq ID')
    plt.ylabel('Counts')
    plt.legend(loc='best')
    plt.show()

    """
    # GC/AT contents profile
    freq=nt_freq(good_seqs)
    ATprofile, GCprofile=GC_profile(freq)
    ideal,_=GC_profile(nt_freq(ideal_seqs))
    #bias_AT,bias_GC=GC_profile(nt_freq(bias_seqs))
    plt.plot(range(len(ATprofile)),ATprofile,'r', label='AT frequency')
    plt.plot(range(len(GCprofile)),GCprofile,'b', label='GC frequency')   
    plt.plot(range(len(ideal)),ideal,'k--')
    #plt.plot(range(len(bias_AT)),bias_AT,'r:')
    #plt.plot(range(len(bias_GC)),bias_GC,'b:')
    plt.xlabel('Position (bp)')
    plt.ylabel('frequency')
    plt.legend(loc='upper right')
    plt.show()
    """

    # print 100 best abundance sequences
    insertcount=pick_best(insert_count)
    for i in range(100):
        print "%s\t%d" % (insertcount[i][0], insertcount[i][1])
    
    """
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
    parser = ArgumentParser(description='graphical analysis of insert')
    parser.add_argument('-f',
                        dest = 'sort_fname',
                        type=str,
                        help='input filename')
    parser.add_argument('-x',
                        dest = 'ref_fname',
                        type=str,
                        help='input filename')

    args = parser.parse_args()
                        
    display_graph (args.sort_fname,
                   args.ref_fname)
