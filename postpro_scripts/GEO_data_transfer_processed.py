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
    cmd = ['./rclone',
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

# rclone_copyto
def rclone_copyto (from_path,
                   to_path,
                   fname):

    print 'Copying %s' % (fname)

    today = date.today().strftime("%b-%d-%Y")
    if len(glob.glob('./' + today)) == 0:
        subprocess.call(['mkdir', today])

    cmd = ['./rclone',
           'copyto',
           '--ignore-existing',
           from_path + '/' + fname,
           to_path + '/' + fname]

    log_name = '-'.join([from_path.strip('/').split('/')[-1],fname]) 

    subprocess.call(cmd,
                    stdout=open("./" + today + "/" + log_name + "_out.txt", 'w'),
                    stderr=open("./" + today + "/" + log_name + "_err.txt", 'w'))
    
    return


### read Condense-seq data table
exp_IDs = read_table('Condense-seq NGS data table revision.csv')
exps = sorted(exp_IDs.keys(), cmp=exp_cmp)

### transfer processed files into GEO
from_path = '/home/spark159/data/2024_01_05_GEO/processed_files'
to_path = 'GEO_ftp:uploads/sparkly205@gmail.com_AhGDIowW/processed_files_update'

for exp in exps:
    ID, qn, ID_deep, qn_deep = exp_IDs[exp]
    rep, cell, sample, agent, tnum = exp

    processed_files = []
    fhead = '_'.join([cell, sample, agent, str(rep)+'rep'])
    processed_files.append(fhead + '_10kb_score.gtab.gz')
    if ID_deep != '':
        processed_files.append(fhead + '_deep_1kb_score.gtab.gz')
        processed_files.append(fhead + '_deep_score.tar')

    for fname in processed_files:
        rclone_copyto (from_path,
                       to_path,
                       fname)
        
    print
