import subprocess
import glob
import sys, re
import copy
import csv, codecs
import gzip
import math

### read scale table
def read_scale_table (fname):
    exp_scale = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split()

        if First:
            mark = cols[0]
            labels = cols[1:]
            First = False
            continue

        condition = cols[0]
        scales = cols[1:]
        for i in range(len(scales)):
            label = labels[i]
            scale = int(scales[i])
            exp = (mark, condition, label)
            exp_scale[exp] = scale
    return exp_scale    

### ref information
org_ref = {'human':"4D_hg38", 'mouse':"4D_mm10"}
cell_org = {'H1':'human', 'GM':'human', 'mCD8T:WT':'mouse', 'mCD8T:DFMO':'mouse', 'mCD8T:ODCKO':'mouse'}
cell_chrnames = {'H1':['chr%s' % (i) for i in range(1, 23)] + ['chrX', 'chrY'],
                 'GM':['chr%s' % (i) for i in range(1, 23)] + ['chrX'],
                 'mCD8T:WT':['chr%s' % (i) for i in range(1, 20)] + ['chrX'],
                 'mCD8T:DFMO':['chr%s' % (i) for i in range(1, 20)] + ['chrX'],
                 'mCD8T:ODCKO':['chr%s' % (i) for i in range(1, 20)] + ['chrX']}

### path file information
ref_path = '/home/spark159/data/2024_01_05_GEO/ref_files'
bam_path = '/home/spark159/data/MouseEpigeneticData/ODC_mouse_HistoneChipseq'
reps = [1, 2, 3]

### other parmeters
cell = 'mCD8T:WT'
refname = ref_path + '/' + org_ref[cell_org[cell]]
chr_names = cell_chrnames[cell]

### read scale table
exp_scale = {}
table_fnames = ['H3K27ac_spike.tsv', 'H3K27me3_spike.tsv']
for table_fname in table_fnames:
    exp_spike = read_scale_table(bam_path + '/' + table_fname)
    min_spike = min(exp_spike.values())
    for exp, spike in exp_spike.items():
        scale = min_spike/float(spike)
        exp_scale[exp] = scale

## batch readcount submit
for rep in reps:
    bamfiles = glob.glob(bam_path + '/' + 'Rep-' + str(rep) + '/02_bams_mm10/*.bam')
    for bamfile in bamfiles:
        cols = bamfile.split('/')[-1].rsplit('.', 1)[0].split('_')
        print cols
        if rep == 1:
            condition, mark = cols[1], cols[2]
        else:
            mark, condition = cols[0], cols[1]

        if mark == 'H3k27ac':
            mark = 'H3K27ac'
        elif mark == 'H3k27me3':
            mark = 'H3K27me3'

        exp = (mark, condition, 'Rep' + str(rep))
        print exp
        scale = exp_scale[exp]

        out_path = bam_path + '/' + 'gtab_files/'
        outfname = '_'.join([mark, condition, str(rep) + 'rep'])

        subprocess.call(["sbatch",
                         "readcount-submit",
                         "-f", bamfile,
                         "-x", refname+'.fa',
                         "-c", ','.join(chr_names),
                         "-s", str(scale),
                         "-o", out_path + outfname])
