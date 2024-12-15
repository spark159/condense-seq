import sys
import glob
import subprocess
from datetime import date

# cell information (organism, cell line, genotype, treatment)
cell_info = {'H1':('Homo sapiens', 'H1-hESC', 'WT', 'None'),
             'GM':('Homo sapiens', 'GM12878', 'WT', 'None'),
             'mCD8T:WT':('Mus musculus', 'CD8+ T cell', 'WT', 'None'),
             'mCD8T:DFMO':('Mus musculus', 'CD8+ T cell', 'WT', 'Difluoromethylornithine (DFMO) treated'),
             'mCD8T:ODCKO':('Mus musculus', 'CD8+ T cell', 'Ornithine decarboxylase (ODC) knockout', 'None'),
             'E14':('Mus musculus', 'E14 mESC', 'WT', 'None')}

# sample type fullnames
sample_fullname = {'NCP':'Native nucleosome',
                   'synNCP': 'Reconstituted nucleosome',
                   'DNA':'Nucleosomal DNA'}

# agent fullnames
agent_fullname = {'sp':'Spermine',
                  'spd':'Spermidine',
                  'CoH':'Cobalt Hexammine',
                  'PEG':'PEG 8000',
                  'Ca':'Calcium',
                  'HP1a':'HP1-alpha ',
                  'HP1bSUV':'HP1-beta/SUV39H1'}

# sort exp tuple (rep, cell, sample, agent, tnum)
def exp_cmp (exp1, exp2):
    cell_order = {'H1':0, 'GM':1, 'mCD8T:WT':2, 'mCD8T:DFMO':3, 'mCD8T:ODCKO':4, 'E14':5}
    sample_order = {'NCP':0, 'synNCP':2, 'DNA':3}
    agent_order = {'sp':0, 'spd':1, 'CoH':2, 'PEG':3, 'Ca':4, 'HP1a':5, 'HP1bSUV':6}

    rep1, cell1, sample1, agent1, tnum1 = exp1
    rep2, cell2, sample2, agent2, tnum2 = exp2
    if rep1 < rep2:
        return -1
    elif rep1 > rep2:
        return 1

    if cell_order[cell1] < cell_order[cell2]:
        return -1
    elif cell_order[cell1] > cell_order[cell2]:
        return 1

    if sample_order[sample1] < sample_order[sample2]:
        return -1
    elif sample_order[sample1] > sample_order[sample2]:
        return 1

    if agent_order[agent1] < agent_order[agent2]:
        return -1
    elif agent_order[agent1] > agent_order[agent2]:
        return 1

    if tnum1 < tnum2:
        return -1
    elif tnum1 > tnum2:
        return 1
    else:
        return 0
    
# read Condense-seq data table and get project ID and fastq file names
def read_table (fname):
    exp_IDs = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        if First:
            First = False
            cols = line.split(',')
            #print cols[-10:-2]
            continue
        cols = line.split(',')
        libname = cols[0]
        cell, sample, agent, tnum = libname.split('_')
        #print cols[-10:-2]
        ID_1rep, qn_1rep, ID_1rep_deep, qn_1rep_deep = cols[-10:-6]
        ID_2rep, qn_2rep, ID_2rep_deep, qn_2rep_deep = cols[-6:-2]

        exp_IDs[(1, cell, sample, agent, int(tnum))] = [ID_1rep, qn_1rep, ID_1rep_deep, qn_1rep_deep]
        exp_IDs[(2, cell, sample, agent, int(tnum))] = [ID_2rep, qn_2rep, ID_2rep_deep, qn_2rep_deep]
    return exp_IDs

# get all fastq files in the path
def get_fastq (path, header):
    cmd = ['rclone',
           'ls',
           '--include',
           header + '*.fastq.gz',
           path]
    
    rclone_ls = subprocess.Popen(cmd,
                                 stdout=subprocess.PIPE,
                                 stderr=open("/dev/null", 'w'))
    fnames = []
    for line in rclone_ls.stdout:
        if not line:
            continue
        _, fname = line.strip().split(' ')
        fnames.append(fname)
    return fnames

# segregate fname list to the tuple list of read pairs
def make_pairs (fnames):
    assert len(fnames) !=0
    assert len(fnames) % 2 ==0
    fname_pairs = []
    for i in range(len(fnames)/2):
        f1, f2 = fnames[2*i], fnames[2*i+1]
        head1, head2 = f1.rsplit('.', 2)[0], f2.rsplit('.', 2)[0]
        assert head1.endswith('R1_001')
        assert head2.endswith('R2_001')
        assert head1.split('R1_001')[0] == head2.split('R2_001')[0] 
        fname_pairs.append((f1, f2))
    return fname_pairs


### read Condense-seq data table
exp_IDs = read_table('../Condense-seq NGS data table revision.csv')
exps = sorted(exp_IDs.keys(), cmp=exp_cmp)

### list out fastq files for each exp
data_path = 'GEO_ftp:uploads/sparkly205@gmail.com_AhGDIowW/fastq_files'
exp_fastq = {}
for exp in exps:
    rep, cell, sample, agent, tnum = exp
    fhead = '_'.join([cell, sample, agent, str(tnum), str(rep)])

    fnames, fnames_deep = [], []
    for fname in get_fastq(data_path, fhead):
        if 'deep' in fname:
            fnames_deep.append(fname)
        else:
            fnames.append(fname)
    fnames = sorted(fnames)
    fnames_deep = sorted(fnames_deep)
    
    fname_pairs = make_pairs(sorted(fnames))
    if len(fnames_deep) > 0:
        fname_pairs_deep = make_pairs(sorted(fnames_deep))
    else:
        fname_pairs_deep = []

    exp_fastq[exp] = {}
    exp_fastq[exp]['shallow'] = fname_pairs
    exp_fastq[exp]['deep'] = fname_pairs_deep

# make GEO metadata table
fname = 'GEO_table1.txt'
f = open(fname, 'w')

cols = ['library name',
        'title',
        'organism',
        'cell line',
        'genotype',
        'treatment',
        'molecule',
        'single or pair-end',
        'instrument model',
        'description']
cols += ['processed data file' + str(i+1) for i in range(3)]
cols += ['raw file' + str(i+1) for i in range(12)]

print >> f, '\t'.join(cols)

data_path = ''
for exp in exps:
    rep, cell, sample, agent, tnum = exp
    
    libname = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep'])
    title = ','.join([cell, sample, agent, str(tnum), str(rep)+'rep'])
    org, cell_line, genotype, treatment = cell_info[cell]
    molecule = 'Genomic DNA'
    seq_type = 'pair-end'
    instrument = 'Illumina NovaSeq 6000'
    description = '%s used as condensing agent' % agent_fullname[agent]

    row = [libname, title, org, cell_line, genotype, treatment, molecule, seq_type, instrument, description]

    processed_files = []
    fhead = '_'.join([cell, sample, agent, str(rep)+'rep'])
    processed_files.append(fhead + '_10kb_score.gtab.gz')
    if len(exp_fastq[exp]['deep']) > 0:
        processed_files.append(fhead + '_deep_1kb_score.gtab.gz')
        processed_files.append(fhead + '_deep_score.tar')
    else:
        processed_files += ['']*2
    row += processed_files
    
    raw_files = []
    for fpair in exp_fastq[exp]['shallow']:
        raw_files += list(fpair)
    for fpair in exp_fastq[exp]['deep']:
        raw_files += list(fpair)
    raw_files += ['' for k in range(12 - len(raw_files))]
    row += raw_files

    print >> f, '\t'.join(row)
f.close()

f = open('GEO_table2.txt', 'w')
print >> f, '\t'.join(['R1', 'R2'])
for exp in exps:
    for fpair in exp_fastq[exp]['shallow']:
        print >> f, '\t'.join(list(fpair))
    for fpair in exp_fastq[exp]['deep']:
        print >> f, '\t'.join(list(fpair))
f.close()
