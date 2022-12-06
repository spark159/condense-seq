import numpy as np

def read_eigenfile (fname, bin_size=1000000):
    eigen_list = []
    interval_list = []
    i = 0
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        try:
            value = float(line)
        except:
            value = np.nan
        st = i*bin_size
        ed = (i+1)*bin_size
        eigen_list.append(value)
        interval_list.append((st,ed))
        i +=1
    return eigen_list, interval_list

def make_bedgraph (fname, binID_value, bin_size, header=None):
    f = open(fname, 'w')
    if header == None:
        header = 'track type=bedGraph name="BedGraph Format" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20'
        print >> f, header
    for i in range(len(binID_value)):
        st, ed = i*bin_size, (i+1)*bin_size
        value = float(binID_value[i])
        try:
            value = float(value)
        except:
            continue
        if np.isnan(value):
            continue
        print >> f, "chr1\t" + "%d\t%d\t%f" % (st, ed, value)
        #print >> f, "chr1\t" + "%d\t%d\t%f" % (st, ed, -value)
    f.close()

#bin_size = 50000
bin_size = 100000
#bin_size = 1000000

path = ""
#path = "./data/"
#fnames = ['eigen_WT_50kb.txt', 'eigen_CohesinKO_50kb.txt']
#fnames = ["eigen_GM12878_50kb.txt"]
#fnames = ["eigen_mouseCD8Tcell_1Mb.txt"]
fnames = ['eigen_mouseCD8Tcell_100kb.txt']
for i in range(len(fnames)):
    fname = fnames[i]
    eigen_list, interval_list = read_eigenfile(path+fname, bin_size=bin_size)
    new_fname = fname.split('.')[0] + '.bedgraph'
    make_bedgraph(path + new_fname, eigen_list, bin_size=bin_size)
