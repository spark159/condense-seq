import subprocess
import glob
import sys, re
import copy
import csv, codecs

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

# read titration file
def read_titration (fname):
    tnum_tfrac = {}
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        cols = line.strip().split()
        try:
            tnum = int(cols[-1])
        except:
            continue
        total_frac = float(cols[-3])
        assert tnum not in tnum_tfrac
        tnum_tfrac[tnum] = total_frac
    return tnum_tfrac


# read epigenetic information csv file
def read_csv_file (fname,
                   mode='row',
                   delim=',',
                   header=True,
                   rowID=True,
                   jump=None):
    if rowID:
        col_st = 1
    else:
        col_st = 0
        
    ID_field_value = {}
    First = True
    counter = -1

    #for cols in csv.reader(open(fname), delimiter=delim):
    for cols in csv.reader(codecs.EncodedFile(open(fname), 'utf-8', 'utf-8-sig'),
                           delimiter=delim):

        if First and header:
            field_names = cols[col_st:]
            First = False
            continue
        elif First and not header:
            field_names = range(len(cols[col_st:]))
            First = False
            pass

        counter += 1
        if jump and counter % jump != 0:
            continue

        if rowID:
            ID = cols[0]
        else:
            ID = counter

        if ID not in ID_field_value:
            ID_field_value[ID] = {}

        cols = cols[col_st:]
        #print cols
        for i in range(len(cols)):
            field = field_names[i]
            try:
                value = float(cols[i])
            except:
                value = cols[i]
            if field not in ID_field_value[ID]:
                ID_field_value[ID][field] = value
            else:
                if type(ID_field_value[ID][field]) != list:
                    ID_field_value[ID][field] = [ID_field_value[ID][field]]
                ID_field_value[ID][field].append(value)

    if mode == 'row':
        return ID_field_value

    if mode == 'col' or mode == 'both':
        field_ID_value = {}
        for ID in ID_field_value:
            field_value = ID_field_value[ID]
            for field in field_value:
                value = field_value[field]
                if field not in field_ID_value:
                    field_ID_value[field] = {}
                field_ID_value[field][ID] = value

    if mode == 'col':
        return field_ID_value

    if mode == 'both':
        return ID_field_value, field_ID_value


### ref information
org_ref = {'human':"4D_hg38", 'mouse':"4D_mm10"}
cell_org = {'H1':'human', 'GM':'human', 'mCD8T:WT':'mouse', 'mCD8T:DFMO':'mouse', 'mCD8T:ODCKO':'mouse'}
cell_chrnames = {'H1':['chr%s' % (i) for i in range(1, 23)] + ['chrX', 'chrY'],
                 'GM':['chr%s' % (i) for i in range(1, 23)] + ['chrX'],
                 'mCD8T:WT':['chr%s' % (i) for i in range(1, 20)] + ['chrX'],
                 'mCD8T:DFMO':['chr%s' % (i) for i in range(1, 20)] + ['chrX'],
                 'mCD8T:ODCKO':['chr%s' % (i) for i in range(1, 20)] + ['chrX']}

### path information
data_path = '/home/spark159/data'
ref_path = '/home/spark159/data/2024_01_05_GEO/ref_files'
bam_path = '/home/spark159/data/2024_01_05_GEO/bam_files'
titr_path = '/home/spark159/data/2024_01_05_GEO/titration_files'
output_path = '/home/spark159/data/2024_01_05_GEO/processed_files'
org_anot_path = {'human':'/home/spark159/data/HumanEpigeneticData',
                 'mouse':'/home/spark159/data/MouseEpigeneticData'}

### other parameters
trim_to = 50 # trim fastq reads to make 50 bp
bin_sizes = [1000, 10000] # genomic bin size in bp for counting reads (shallow)
Nlen = 171 # NCP length window for computing coverage area (deep)
sliding_windows = [(167, 25)] # bin size and step for sliding windows (deep)
# input molecule number scale for each agent
agent_mscale = {'sp':1, 'spd':1, 'CoH':1, 'PEG':1, 'Ca':5, 'HP1a':1, 'HP1bSUV':1}
# logistic fitting parameters


### read condense-seq data table
table_fname = 'Condense-seq NGS data table.csv'
exp_IDs = read_table(bam_path + '/' + table_fname)
exps = sorted(exp_IDs.keys(), cmp=exp_cmp)

### gather all fastq files
if False:
    exp_fastq = {}
    for exp in exps:
        ID, qn, ID_deep, qn_deep = exp_IDs[exp]
        fnames = glob.glob(data_path+'/'+ID+'/'+qn+'*.fastq.gz')
        fnames = sorted([fname.rsplit('/', 1)[-1] for fname in fnames])
        fname_pairs = make_pairs(fnames)

        if ID_deep != '':
            fnames = glob.glob(data_path+'/'+ID_deep+'/'+qn_deep+'*.fastq.gz')
            fnames_deep = sorted([fname.rsplit('/', 1)[-1] for fname in fnames])
            fname_pairs_deep = make_pairs(fnames_deep)      
        else:
            fnames_deep = []
            fname_pairs_deep = []

        exp_fastq[exp] = {}
        exp_fastq[exp]['shallow'] = fname_pairs
        exp_fastq[exp]['deep'] = fname_pairs_deep

    
### Condense-seq analysis for each exp
## Bowtie2 alignment
if False:
    for exp in exps:
        rep, cell, sample, agent, tnum = exp
        ID, qn, ID_deep, qn_deep = exp_IDs[exp]

        fheader = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep'])
        refname = ref_path + '/' + org_ref[cell_org[cell]]

        for depth in ['shallow', 'deep']:            
            fname_pairs = exp_fastq[exp][depth]

            if len(fname_pairs) == 0:
                continue

            if depth == 'shallow':
                folder = ID
            elif depth == 'deep':
                folder = ID_deep

            fnames1, fnames2 = [], []
            for fname1, fname2 in fname_pairs:
                fnames1.append(data_path + '/' + folder + '/' + fname1)
                fnames2.append(data_path + '/' + folder + '/' + fname2)

            outfname = bam_path + '/' + fheader
            if depth == 'deep':
                outfname += '_deep'

            subprocess.call(["sbatch",
                             "bowtie-submit_ver2",
                             "-f", ','.join(fnames1),
                             "-g", ','.join(fnames2),
                             "-t", str(trim_to),
                             "-x", refname,
                             "-o", outfname])

## get read length distribution
if False:
    for exp in exps:
        rep, cell, sample, agent, tnum = exp
        ID, qn, ID_deep, qn_deep = exp_IDs[exp]
        fheader = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep'])
        
        for depth in ['shallow', 'deep']:
            if depth == 'deep' and ID_deep == '':
                continue

            if depth == 'deep':
                fheader += '_deep'

            fname = bam_path + '/' + fheader + '.bam'
            outfname = output_path + '/' + fheader

            subprocess.call(["sbatch",
                             "lendist-submit",
                             "-f", fname,
                             "-o", outfname])


## get read counts for each genomic bins
if False:
    for bin_size in bin_sizes:
        for exp in exps:
            rep, cell, sample, agent, tnum = exp
            ID, qn, ID_deep, qn_deep = exp_IDs[exp]
            
            fheader = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep'])
            refname = ref_path + '/' + org_ref[cell_org[cell]]
            chr_names = cell_chrnames[cell]
            
            for depth in ['shallow', 'deep']:
                if depth == 'deep' and ID_deep == '':
                    continue

                if depth == 'deep':
                    fheader += '_deep'

                fname = bam_path + '/' + fheader + '.bam'
                outfname = output_path + '/' + fheader + '_' + str(int(bin_size/1000)) + 'kb'

                subprocess.call(["sbatch",
                                 "bincount-submit",
                                 "-f", fname,
                                 "-c", ','.join(chr_names),
                                 "-o", outfname,
                                 "-x", refname+'.fa',
                                 "-w", str(bin_size)])


## get coverage along the genome (deep only)
if False:
    for exp in exps:
        ID, qn, ID_deep, qn_deep = exp_IDs[exp]
        if ID_deep == '':
            continue

        rep, cell, sample, agent, tnum = exp
        fheader = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep', 'deep'])
        refname = ref_path + '/' + org_ref[cell_org[cell]]

        fname = bam_path + '/' + fheader + '.bam'
        chr_names = cell_chrnames[cell]

        for chr_name in chr_names:    
            outfname = output_path + '/' + fheader + '_' + chr_name

            subprocess.call(["sbatch",
                             "coverage-submit",
                             "-f", fname,
                             "-x", refname+'.fa',
                             "-c", chr_name,
                             "-o", outfname])


## NCP peak calling for input (deep only, input only)
if False:
    for exp in exps:
        ID, qn, ID_deep, qn_deep = exp_IDs[exp]
        if ID_deep == '':
            continue

        rep, cell, sample, agent, tnum = exp        
        if tnum != 0:
            continue
        
        fheader = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep', 'deep'])
        refname = ref_path + '/' + org_ref[cell_org[cell]]

        fname = bam_path + '/' + fheader + '.bam'
        chr_names = cell_chrnames[cell]

        for chr_name in chr_names:
            outfname = output_path + '/' + fheader + '_' + chr_name

            subprocess.call(["sbatch",
                             "NCPpeak-submit",
                             "-f", fname,
                             "-x", refname+'.fa',
                             "-c", chr_name,
                             "-o", outfname])


## compute coverage area under each input NCP peaks (deep only)
if False:
    for exp in exps:
        ID, qn, ID_deep, qn_deep = exp_IDs[exp]
        if ID_deep == '':
            continue

        rep, cell, sample, agent, tnum = exp

        chr_names = cell_chrnames[cell]
        for chr_name in chr_names:
            peak_fheader = '_'.join([cell, sample, agent, str(0), str(rep)+'rep', 'deep', chr_name])
            fheader = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep', 'deep', chr_name])

            peak_fname = output_path + '/' + peak_fheader + '_peak.cn'
            cov_fname = output_path + '/' + fheader + '_cov.cn'
            outfname = output_path + '/' + fheader

            subprocess.call(["sbatch",
                             "NCPcov-submit",
                             "-f", peak_fname,
                             "-g", cov_fname,
                             "-n", str(Nlen),
                             "-c", chr_name,
                             "-o", outfname])


## compute coverage area under each sliding windows (deep only)
if False:
    for exp in exps:
        ID, qn, ID_deep, qn_deep = exp_IDs[exp]
        if ID_deep == '':
            continue

        rep, cell, sample, agent, tnum = exp
        chr_names = cell_chrnames[cell]

        fheader = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep', 'deep'])
        for Bsize, Bstep in sliding_windows:
            for chr_name in chr_names:
                cov_fname = output_path + '/' + fheader + '_' + chr_name + '_' + '_cov.cn'
                outfname = output_path + '/' + fheader + '_' + chr_name
                outfname += '_%dwin%dstep' % (Bsize, Bstep)

                subprocess.call(["sbatch",
                                 "Bindata-submit",
                                 "-f", fname,
                                 "-x", refname+'.fa',
                                 "-c", chr_name,
                                 "-w", str(Bsize),
                                 "-s", str(Bstep),
                                 "-o", outfname])
                

## get condensability score (direct)
if False:
    for exp in exps:
        rep, cell, sample, agent, tnum = exp
        if tnum == 0:
            continue

        ID, qn, ID_deep, qn_deep = exp_IDs[exp]
        tfname = titr_path + '/' + '_'.join([cell, sample, agent, 'titration']) + '.csv'
        chr_names = cell_chrnames[cell]
        
        for bin_size in bin_sizes:
            fheader = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep'])
            input_fheader = '_'.join([cell, sample, agent, str(0), str(rep)+'rep'])
            fname = output_path + '/' + fheader  + '_' + str(int(bin_size/1000)) + 'kb' + '_bin.cn'
            input_fname = output_path + '/' + input_fheader  + '_' + str(int(bin_size/1000)) + 'kb' + '_bin.cn'
            
            subprocess.call(["sbatch",
                             "NCPscore-submit_ver2",
                             "-f", fname,
                             "-i", input_fname,
                             "-g", tfname,
                             "-t", str(tnum),
                             "-c", ','.join(chr_names)])
            
            if ID_deep != '':

                fheader = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep', 'deep'])
                input_fheader = '_'.join([cell, sample, agent, str(0), str(rep)+'rep', 'deep'])
                fname = output_path + '/' + fheader  + '_' + str(int(bin_size/1000)) + 'kb' + '_bin.cn'
                input_fname = output_path + '/' + input_fheader  + '_' + str(int(bin_size/1000)) + 'kb' + '_bin.cn'

                subprocess.call(["sbatch",
                                 "NCPscore-submit_ver2",
                                 "-f", fname,
                                 "-i", input_fname,
                                 "-g", tfname,
                                 "-t", str(tnum),
                                 "-c", ','.join(chr_names)])


        # deeply sequnced case
        if ID_deep != '':
            fheader = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep', 'deep'])
            input_fheader = '_'.join([cell, sample, agent, str(0), str(rep)+'rep', 'deep'])
            
            fnames = []
            input_fnames = []
            for chr_name in chr_names:
                fname = output_path + '/' + fheader + '_' + chr_name + '_Ncov.cn'
                fnames.append(fname)
                input_fname = output_path + '/' + input_fheader + '_' + chr_name + '_Ncov.cn'
                input_fnames.append(input_fname)

            subprocess.call(["sbatch",
                             "NCPscore-submit_ver2",
                             "-f", ','.join(fnames),
                             "-i", ','.join(input_fnames),
                             "-g", tfname,
                             "-t", str(tnum),
                             "-c", ','.join(chr_names)])


    


## estimate physical NCP number
if False:
    for exp in exps:
        ID, qn, ID_deep, qn_deep = exp_IDs[exp]
        rep, cell, sample, agent, tnum = exp
        mscale = agent_mscale[agent]
        tfname = titr_path + '/' + '_'.join([cell, sample, agent, 'titration']) + '.csv'
        chr_names = cell_chrnames[cell]

        for bin_size in bin_sizes:
            fheader = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep'])
            fname = output_path + '/' + fheader  + '_' + str(int(bin_size/1000)) + 'kb' + '_bin.cn'
            subprocess.call(["sbatch",
                             "NCPnum-submit",
                             "-f", fname,
                             "-g", tfname,
                             "-t", str(tnum),
                             "-c", ','.join(chr_names),
                             "-i", str(mscale)])

            if ID_deep != '':
                fheader = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep', 'deep'])
                fname = output_path + '/' + fheader  + '_' + str(int(bin_size/1000)) + 'kb' + '_bin.cn'
                subprocess.call(["sbatch",
                                 "NCPnum-submit",
                                 "-f", fname,
                                 "-g", tfname,
                                 "-t", str(tnum),
                                 "-c", ','.join(chr_names),
                                 "-i", str(mscale)])

        # deeply sequnced case
        if ID_deep != '':
            fheader = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep', 'deep'])
            fnames = []
            chr_names = cell_chrnames[cell]
            for chr_name in chr_names:
                fname = output_path + '/' + fheader + '_' + chr_name + '_Ncov.cn'
                fnames.append(fname)

            subprocess.call(["sbatch",
                             "NCPnum-submit",
                             "-f", ','.join(fnames),
                             "-g", tfname,
                             "-t", str(tnum),
                             "-c", ','.join(chr_names),
                             "-i", str(mscale)])


## concatenate num/score files over titrations
if False:
    # group exps over titrations
    dexp_tnums = {}
    dexp_tnums_deep = {}
    for exp in exps:
        rep, cell, sample, agent, tnum = exp
        dexp = (rep, cell, sample, agent)
        if dexp not in dexp_tnums:
            dexp_tnums[dexp] = []
        dexp_tnums[dexp].append(tnum)

        ID, qn, ID_deep, qn_deep = exp_IDs[exp]
        if ID_deep != '':
            if dexp not in dexp_tnums_deep:
                dexp_tnums_deep[dexp] = []
            dexp_tnums_deep[dexp].append(tnum)

    # concatenate num/score files (shallowly sequenced)
    for extension in ['_num.cn', '_score.cn']:
        for dexp, tnums in dexp_tnums.items():
            rep, cell, sample, agent = dexp
            tnums = copy.deepcopy(sorted(tnums))
            
            if extension == '_score.cn':
                assert tnums[0] == 0
                tnums = tnums[1:]

            for bin_size in bin_sizes:
                fnames = []
                colnums = []
                for i in range(len(tnums)):
                    tnum = tnums[i]
                    fheader = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep'])
                    fname = output_path + '/' + fheader  + '_' + str(int(bin_size/1000)) + 'kb' + extension
                    if i == 0:
                        fnames += [fname]*4
                        colnums += range(1,5)
                    fnames.append(fname)
                    colnums.append(5)

                fheader = '_'.join([cell, sample, agent, str(rep)+'rep'])
                outfname = output_path + '/' + fheader  + '_' + str(int(bin_size/1000)) + 'kb' + extension
                colnums = [str(colnum) for colnum in colnums]
                
                subprocess.call(["sbatch",
                                 "concat-files",
                                 "-f", ','.join(fnames),
                                 "-c", ','.join(colnums),
                                 "-o", outfname])

    # concatenate num/score files (deeply sequenced)
    for extension in ['_num.cn', '_score.cn']:
        for dexp, tnums in dexp_tnums_deep.items():
            rep, cell, sample, agent = dexp
            tnums = copy.deepcopy(sorted(tnums))

            if extension == '_score.cn':
                assert tnums[0] == 0
                tnums = tnums[1:]

            for bin_size in bin_sizes:
                fnames = []
                colnums = []
                for i in range(len(tnums)):
                    tnum = tnums[i]
                    fheader = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep', 'deep'])
                    fname = output_path + '/' + fheader  + '_' + str(int(bin_size/1000)) + 'kb' + extension
                    if i == 0:
                        fnames += [fname]*4
                        colnums += range(1,5)
                    fnames.append(fname)
                    colnums.append(5)

                fheader = '_'.join([cell, sample, agent, str(rep)+'rep', 'deep'])
                outfname = output_path + '/' + fheader  + '_' + str(int(bin_size/1000)) + 'kb' + extension
                colnums = [str(colnum) for colnum in colnums]
                
                subprocess.call(["sbatch",
                                 "concat-files",
                                 "-f", ','.join(fnames),
                                 "-c", ','.join(colnums),
                                 "-o", outfname])
            

            chr_names = cell_chrnames[cell]
            for chr_name in chr_names:
                fnames = []
                colnums = []
                for i in range(len(tnums)):
                    tnum = tnums[i]
                    fheader = '_'.join([cell, sample, agent, str(tnum), str(rep)+'rep','deep'])
                    fname = output_path + '/' + fheader + '_' + chr_name + extension

                    if i == 0:
                        fnames += [fname]*3
                        colnums += range(1,4)
                    fnames.append(fname)
                    colnums.append(4)

                fheader = '_'.join([cell, sample, agent, str(rep)+'rep', 'deep'])
                outfname = output_path + '/' + fheader  + '_' + chr_name + extension
                colnums = [str(colnum) for colnum in colnums]

                subprocess.call(["sbatch",
                                 "concat-files",
                                 "-f", ','.join(fnames),
                                 "-c", ','.join(colnums),
                                 "-o", outfname])


### gzip num/score files over chromosomes (deeply sequenced)
if False:
    # group exps over titrations
    dexp_tnums_deep = {}
    for exp in exps:
        ID, qn, ID_deep, qn_deep = exp_IDs[exp]
        if ID_deep == '':
            continue
        
        rep, cell, sample, agent, tnum = exp
        dexp = (rep, cell, sample, agent)
       
        if dexp not in dexp_tnums_deep:
            dexp_tnums_deep[dexp] = []
        dexp_tnums_deep[dexp].append(tnum)
    dexps_deep = sorted(dexp_tnums_deep)
    
    # gzip files
    for extension in ['_num.cn', '_score.cn']:
        for dexp in dexps_deep:
            rep, cell, sample, agent = dexp
            chr_names = cell_chrnames[cell]

            fnames = []
            fheader = '_'.join([cell, sample, agent, str(rep)+'rep', 'deep'])
            for chr_name in chr_names:
                fname = output_path + '/' + fheader  + '_' + chr_name + extension
                fnames.append(fname)
            outfname = output_path + '/' + fheader  + extension.rsplit('.', 1)[0] + '.tar.gz'

            tar_cmd = ['tar',
                       '-czvf']

            tar_cmd.append(outfname)
            tar_cmd += fnames

            subprocess.call(tar_cmd)

### logistic regression of concatenated num files over titrations
if False:
    # group exps over titrations
    dexp_tnums = {}
    for exp in exps:
        rep, cell, sample, agent, tnum = exp
        dexp = (rep, cell, sample, agent)
        if dexp not in dexp_tnums:
            dexp_tnums[dexp] = []
        dexp_tnums[dexp].append(tnum)

    # logistic regression for 10kb data set with >5 titration points
    bin_size = 10000
    min_tnum = 5
    min_rsq = 0.5
    models = ['sigmoid', 'hill']
    for dexp, tnums in dexp_tnums.items():
        if len(tnums) < min_tnum:
            continue
        rep, cell, sample, agent = dexp
        tnums = sorted(tnums)
        fheader = '_'.join([cell, sample, agent, str(rep)+'rep'])
        fname = output_path + '/' + fheader  + '_' + str(int(bin_size/1000)) + 'kb' + '_num.cn'
        tfname = titr_path + '/' + '_'.join([cell, sample, agent, 'titration']) + '.csv'

        for model in models:
            
            outfname = output_path + '/' + fheader  + '_' + str(int(bin_size/1000)) + 'kb' + '_' + model
            subprocess.call(["sbatch",
                             "logistic-submit",
                             "-f", ','.join([fname]),
                             "-g", tfname,
                             "-t", ','.join([str(tnum) for tnum in tnums]),
                             "-m", model,
                             "-r", str(min_rsq),
                             "-o", outfname])

### make annotation table
if True:
    # get histone PTM chip-seq/bisulfite information
    cell_mtype_mark_bedfiles = {}
    for org in ['human', 'mouse']:
        for mtype in ['HistoneChipseq', 'BSseq']:
            anot_path = org_anot_path[org] + '/' + mtype
            anot_fname = anot_path + '/' + '%s_%s_info.csv' % (org, mtype)

            try:
                ID_field_value = read_csv_file(anot_fname, rowID=False)
            except:
                continue

            for ID in ID_field_value:
                mark = ID_field_value[ID]['Name']
                cell = ID_field_value[ID]['Cell']

                # remap cell names
                if cell == 'GM12878':
                    cell = 'GM'
                elif cell == 'H1-hESC':
                    cell = 'H1'
                elif cell == 'Mouse CD8 T cell (invitro activated)':
                    cell = 'mCD8T:WT'
                else:
                    continue

                # find bed file code
                try:
                    file_codes = ID_field_value[ID]['Bed files'].split(',')
                except:
                    file_codes = ID_field_value[ID]['Peak files'].split(',')

                # map file code to bed file
                bed_files = []
                for file_code in file_codes:
                    search = glob.glob(anot_path + '/' + file_code + '.bed*')
                    assert len(search) == 1
                    bed_file = search[0]
                    bed_files.append(bed_file)

                if cell not in cell_mtype_mark_bedfiles:
                    cell_mtype_mark_bedfiles[cell] = {}
                if mtype not in cell_mtype_mark_bedfiles[cell]:
                    cell_mtype_mark_bedfiles[cell][mtype] = {}
                if mark not in cell_mtype_mark_bedfiles[cell][mtype]:
                    cell_mtype_mark_bedfiles[cell][mtype][mark] = []
                cell_mtype_mark_bedfiles[cell][mtype][mark] += bed_files
                  
    # group exps over titrations
    dexp_tnums = {}
    dexp_tnums_deep = {}
    for exp in exps:
        rep, cell, sample, agent, tnum = exp
        dexp = (rep, cell, sample, agent)
        if dexp not in dexp_tnums:
            dexp_tnums[dexp] = []
        dexp_tnums[dexp].append(tnum)

        ID, qn, ID_deep, qn_deep = exp_IDs[exp]
        if ID_deep != '':
            if dexp not in dexp_tnums_deep:
                dexp_tnums_deep[dexp] = []
            dexp_tnums_deep[dexp].append(tnum)

    # make annotation table
    for dexp in dexp_tnums:
        rep, cell, sample, agent = dexp
            
        chip_input = []
        try:
            for mark, bedfiles in cell_mtype_mark_bedfiles[cell]['HistoneChipseq'].items():
                chip_input.append(','.join([mark] + bedfiles))
        except:
            pass
        
        bs_input = []
        try:
            for mark, bedfiles in cell_mtype_mark_bedfiles[cell]['BSseq'].items():
                bs_input.append(','.join([mark] + bedfiles))
        except:
            pass

        if len(chip_input + bs_input) <=0:
            continue

        refname = ref_path + '/' + org_ref[cell_org[cell]]
        chr_names = cell_chrnames[cell]
        fheader = '_'.join([cell, sample, agent, str(rep)+'rep'])
        extension = '_score.cn'

        # shallow case
        for bin_size in bin_sizes:
            fname = output_path + '/' + fheader  + '_' + str(int(bin_size/1000)) + 'kb' + extension
            outfname = output_path + '/' + fheader  + '_' + str(int(bin_size/1000)) + 'kb' + extension.rsplit('.',1)[0]
            
            subprocess.call(["sbatch",
                             "maketable_submit",
                             "-f", fname,
                             "-x", refname+'.fa',
                             "-c", ','.join(chr_names),
                             "-o", outfname,
                             "-w", str(bin_size),
                             "-s", str(bin_size),
                             "-b", ':'.join(bs_input),
                             "-h", ':'.join(chip_input)])

        # deep case
        if dexp in dexp_tnums_deep:
            for chr_name in chr_names:
                fname = output_path + '/' + fheader  + '_deep_' + chr_name + extension
                outfname = output_path + '/' + fheader  + '_deep_' + chr_name + extension.rsplit('.',1)[0]
                
                subprocess.call(["sbatch",
                                 "maketable_submit",
                                 "-f", fname,
                                 "-x", refname+'.fa',
                                 "-c", chr_name,
                                 "-o", outfname,
                                 "-w", str(Nlen),
                                 "-b", ':'.join(bs_input),
                                 "-h", ':'.join(chip_input)])

### make profile
