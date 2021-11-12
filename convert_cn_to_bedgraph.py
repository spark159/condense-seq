import numpy as np


# file name
path = ""
#infname = path + "H1_NCP_sp_chr1_167win25step_anot.cn"
#outfname = path + "H1_NCP_sp_chr1_167win25step_score2.bedgraph"

infname = path + "H1_NCP_sp_chr1_167win25step_rlen_Bsig.cn"
outfname = path + "H1_NCP_sp_chr1_167win25step_rlen_Bsig.bedgraph"


# set parameters
bin_size = 167
step_size = 25

# data selection
#column = 3 # score1
#column = 4 # score2
column = 5 # control

# bedgraph file paremters
#header = 'track type=bedGraph name="score2" description="score" visibility=full doWiggle=on'
header = 'track type=bedGraph name="readlen" description="readlen" visibility=full doWiggle=on'
#header = None


# reading and writing file
f = open(outfname, 'w')
if header:
    print >> f, header

First = True
for line in open(infname):
    if First:
        features = line.strip().split()
        First = False
        continue
    cols = line.strip().split()
    ID = int(cols[0])
    chr, pos = cols[1], int(cols[2])
    #st = pos - bin_size/2
    #ed = st + bin_size
    st, ed = pos, pos+1  
    try:
        value = float(cols[column])
    except:
        continue
    print >> f, "chr1\t" + "%d\t%d\t%f" % (st, ed, value)
f.close()
