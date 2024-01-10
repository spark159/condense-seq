import sys
import glob
import subprocess
from datetime import date

# cell information (organism, cell line, genotype, treatment)
cell_info = {'H1':('Homo sapiens', 'H1-hESC', 'WT', 'None'),
             'GM':('Homo sapiens', 'GM12878', 'WT', 'None'),
             'mCD8T:WT':('Mus musculus', 'CD8+ T cell', 'WT', 'None'),
             'mCD8T:DFMO':('Mus musculus', 'CD8+ T cell', 'WT', 'Difluoromethylornithine (DFMO) treated'),
             'mCD8T:ODCKO':('Mus musculus', 'CD8+ T cell', 'Ornithine decarboxylase (ODC) knockout', 'None')}

# sample type fullnames
sample_fullname = {'NCP':'Native mono-nucleosome',
                   'DNA':'Nucleosomal DNA'}
# agent fullnames
agent_fullname = {'sp':'Spermine',
                  'spd':'Spermidine',
                  'CoH':'Cobalt Hexammine',
                  'PEG':'PEG 8000',
                  'Ca':'Calcium',
                  'HP1a':'HP1$\\alpha$',
                  'HP1bSUV':'HP1$\\beta$/SUV39H1'}

# sort exp tuple (rep, cell, sample, agent, tnum)
def exp_cmp (exp1, exp2):
    cell_order = {'H1':0, 'GM':1, 'mCD8T:WT':2, 'mCD8T:DFMO':3, 'mCD8T:ODCKO':4}
    sample_order = {'NCP':0, 'DNA':1}
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

exp_IDs = read_table('Condense-seq NGS data table.csv')
exps = sorted(exp_IDs.keys(), cmp=exp_cmp)


# list out fastq files for each exp
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

data_path = 'my_dropbox:My PhD works/Research projects/Condense-seq_project/data/fastq_files'
exp_fastq = {}
exp_newfastq = {}
for exp in exps:
    ID, qn, ID_deep, qn_deep = exp_IDs[exp]
    fnames = sorted(get_fastq(data_path + '/' + ID, qn))
    fname_pairs = make_pairs(fnames)

    if ID_deep != '':
        fnames_deep = sorted(get_fastq(data_path + '/' + ID_deep, qn_deep))
        fname_pairs_deep = make_pairs(fnames_deep)
    else:
        fnames_deep = []
        fname_pairs_deep = []

    exp_fastq[exp] = {}
    exp_fastq[exp]['shallow'] = fname_pairs
    exp_fastq[exp]['deep'] = fname_pairs_deep

    # make new file name
    assert exp not in exp_newfastq
    exp_newfastq[exp] = {}
    exp_newfastq[exp]['shallow'] = []
    exp_newfastq[exp]['deep'] = []
        
    rep, cell, sample, agent, tnum = exp
    new_fhead = [cell, sample, agent, str(tnum), str(rep)+'rep']
    for i in range(len(fname_pairs)):
        old_f1, old_f2 = fname_pairs[i]
        assert old_f1.startswith(qn)
        assert old_f2.startswith(qn)
        lnum = i+1 # lane number
        new_f1 = '_'.join(new_fhead + ['L' + str(lnum), 'R1_001']) + '.fastq.gz'
        new_f2 = '_'.join(new_fhead + ['L' + str(lnum), 'R2_001']) + '.fastq.gz'
        exp_newfastq[exp]['shallow'].append((new_f1, new_f2))

    new_fhead += ['deep']
    for i in range(len(fname_pairs_deep)):
        old_f1, old_f2 = fname_pairs_deep[i]
        assert old_f1.startswith(qn_deep)
        assert old_f2.startswith(qn_deep)
        lnum = i+1 # lane number
        new_f1 = '_'.join(new_fhead + ['L' + str(lnum), 'R1_001']) + '.fastq.gz'
        new_f2 = '_'.join(new_fhead + ['L' + str(lnum), 'R2_001']) + '.fastq.gz'
        exp_newfastq[exp]['deep'].append((new_f1, new_f2))

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
        'description',
        'processed data file']

cols += ['raw file ' + str(i+1) for i in range(12)]

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
    
    processed_file = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep']) +'.num.gz'
    
    raw_files = []
    for fpair in exp_newfastq[exp]['shallow']:
        raw_files += list(fpair)
    for fpair in exp_newfastq[exp]['deep']:
        raw_files += list(fpair)
    raw_files += ['' for k in range(12 - len(raw_files))]

    row = [libname, title, org, cell_line, genotype, treatment, molecule, seq_type, instrument, description, processed_file]
    row += raw_files

    print >> f, '\t'.join(row)
f.close()

f = open('GEO_table2.txt', 'w')
print >> f, '\t'.join(['R1', 'R2'])
for exp in exps:
    for fpair in exp_newfastq[exp]['shallow']:
        print >> f, '\t'.join(list(fpair))
    for fpair in exp_newfastq[exp]['deep']:
        print >> f, '\t'.join(list(fpair))
f.close()

#sys.exit(1)

# upload fastq files into GEO ftp server
# rclone_copyo
def rclone_copyto (from_path,
                   to_path,
                   old_fname,
                   new_fname):

    print 'Copy %s to %s' % (old_fname, new_fname)

    today = date.today().strftime("%b-%d-%Y")
    if len(glob.glob('./' + today)) == 0:
        subprocess.call(['mkdir', today])

    cmd = ['rclone',
           'copyto',
           '--ignore-existing',
           from_path + '/' + old_fname,
           to_path + '/' + new_fname]

    log_name = '-'.join([from_path.strip('/').split('/')[-1],old_fname]) 

    subprocess.call(cmd,
                    stdout=open("./" + today + "/" + log_name + "_out.txt", 'w'),
                    stderr=open("./" + today + "/" + log_name + "_err.txt", 'w'))
    
    return

from_path = 'my_dropbox:My PhD works/Research projects/Condense-seq_project/data/fastq_files'
to_path = 'GEO_ftp:uploads/sparkly205@gmail.com_AhGDIowW/fastq_files'

for exp in exps:
    ID, qn, ID_deep, qn_deep = exp_IDs[exp]
    fpairs = exp_fastq[exp]['shallow']
    fpairs_deep = exp_fastq[exp]['deep']
    new_fpairs = exp_newfastq[exp]['shallow']
    new_fpairs_deep = exp_newfastq[exp]['deep']
    rep, cell, sample, agent, tnum = exp
    new_fhead = [cell, sample, agent, str(tnum), str(rep)+'rep']
    print 'Data transfer: %s' % '_'.join(new_fhead)
    for i in range(len(fpairs)):
        old_f1, old_f2 = fpairs[i]
        new_f1, new_f2 = new_fpairs[i]
        rclone_copyto (from_path + '/' + ID,
                       to_path,
                       old_f1,
                       new_f1)
        rclone_copyto (from_path + '/' + ID,
                       to_path,
                       old_f2,
                       new_f2)
    for i in range(len(fpairs_deep)):
        old_f1, old_f2 = fpairs_deep[i]
        new_f1, new_f2 = new_fpairs_deep[i]
        rclone_copyto (from_path + '/' + ID_deep,
                       to_path,
                       old_f1,
                       new_f1)
        rclone_copyto (from_path + '/' + ID_deep,
                       to_path,
                       old_f2,
                       new_f2)
    print
