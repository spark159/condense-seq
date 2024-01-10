import subprocess

# parameters
#path = "/home/spark159/../../media/spark159/sw/"
path = "/home/spark159/../../storage/"

# experiment list (cell, sample, agent)
exp_list = [('H1', 'NCP', 'sp', 8),
            ('H1', 'NCP', 'HP1a', 3),
            ('H1', 'NCP', 'LKH', 3),
            ('H1', 'NCP', 'Ki67', 4)]

#exp_list = [('H1', 'NCP', 'sp', 8)]

# set parameters
#bin_size = 171
#bin_size = 3
bin_size = 1000

# data type
dtype = 'score'

# header choice
header = False

# write bedgraph file
if True:
    for cell, sample, agent, tnum in exp_list:

        print "Converting %s %s %s" % (cell, sample, agent)

        infname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', dtype]) + '.cn'

        outfname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', dtype]) + '.bedGraph'
        fieldname = '-'.join([cell, sample, agent, str(tnum)]) + '.bam'

        # reading and writing file
        # bedgraph file paremters
        label = '-'.join([cell, sample, agent])

        if header:
            header = 'track type=bedGraph name="%s" description="%s" visibility=full doWiggle=on' % (label, label)

        f = open(outfname, 'w')

        if header:
            print >> f, header

        First = True
        col_pick = None
        for line in open(infname):
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            if First:
                for k in range(len(cols)):
                    if cols[k] == fieldname:
                        col_pick = k
                        break
                First = False
                continue
            ID = int(cols[0])
            chr, st, ed = cols[1], int(cols[2]), int(cols[3])
            value = float(cols[col_pick])
            print >> f, "%s\t%d\t%d\t%f" % (chr, st, ed, value)
        f.close()

# gzip file
if False:
    for cell, sample, agent, tnum in exp_list:
        print "Sorting %s %s %s" % (cell, sample, agent)
 
        infname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', dtype]) + '.bedGraph'

        cmd = ['gzip', infname]
        subprocess.call(cmd)


    
# sort file
if True:
    for cell, sample, agent, tnum in exp_list:
        print "Sorting %s %s %s" % (cell, sample, agent)

        # sort 
        infname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', dtype]) + '.bedGraph'
        outfname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', dtype]) + '.sorted.bedGraph'

        cmd = 'sort -k1,1 -k2,2n %s > %s' % (infname, outfname)
        subprocess.call(cmd, shell=True)

# make bigwig
if True:
    for cell, sample, agent, tnum in exp_list:

        print "Generate bigwig %s %s %s" % (cell, sample, agent)
        
        # sort 
        infname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', dtype]) + '.sorted.bedGraph'
        outfname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', dtype]) + '.bw'
        sizefname = path + "hg38.chrom.sizes.txt"

        cmd = [path+'bedGraphToBigWig', infname, sizefname, outfname]
        subprocess.call(cmd)
