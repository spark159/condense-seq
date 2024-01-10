import subprocess
import glob
import sys, re

### parameters
#path = "/home/spark159/scratch/2022_12_17_Park97lib/"
path = "/home/spark159/scratch/2022_12_17_PlusonelibSD/"

# library name
#lib_choice = 'Park97mmlib'
#lib_choice = 'Park97IDlib'
lib_choice = 'Plusone'

# reference prefix name
if lib_choice == 'Park97IDlib':
    ref_fname = path + 'Park97_Insertion'
elif lib_choice == 'Park97mmlib':
    ref_fname = path + 'Park97_Mismatch'
elif lib_choice == 'Plusone':
    ref_fname = path + 'plusonelib'

# sort mode name
if lib_choice.startswith('Park97'):
    sort_mode = '--direct'
elif lib_choice.startswith('Plusone'):
    sort_mode = '--dw'


# NEB adapter sequence for trimming
adapt_seq1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
adapt_seq2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'


### gather all trimmed fastq files
if True:
    lib_time_rep_pair_fname = {}
    for fname in glob.glob(path + '*.fastq.gz'):
        if fname.startswith(path+"Undetermined"):
            continue
        if 'trim' in fname:
            continue
        if 'combined' in fname:
            continue
        cols = fname.split('/')[-1].split('_')

        if len(cols) == 6:
            lib, time = cols[:2]
            rep = '1rep'

        elif len(cols) == 7:
            lib, time, rep  = cols[:3]
        elif len(cols) == 8:
            libname, libtype, time, timeunit  = cols[:4]
            lib = libname + libtype
            time = time + timeunit
            rep = '1rep'
            
        pair = cols[-2]
        #lane = cols[-3]

        if lib != lib_choice:
            continue

        if lib not in lib_time_rep_pair_fname:
            lib_time_rep_pair_fname[lib] = {}
        if time not in lib_time_rep_pair_fname[lib]:
            lib_time_rep_pair_fname[lib][time] = {}
        if rep not in lib_time_rep_pair_fname[lib][time]:
            lib_time_rep_pair_fname[lib][time][rep] = {}
        lib_time_rep_pair_fname[lib][time][rep][pair] = fname


### trim adapter
if False:    
    lib_time_rep_pair_trimfname = {}
    for lib in lib_time_rep_pair_fname:
        for time in lib_time_rep_pair_fname[lib]:
            for rep in lib_time_rep_pair_fname[lib][time]:
                pair_fname = lib_time_rep_pair_fname[lib][time][rep]
                read1_fname = pair_fname['R1']
                read2_fname = pair_fname['R2']
                read1_outfname = read1_fname.replace('.fastq.gz', '_trim.fastq.gz')
                read2_outfname = read2_fname.replace('.fastq.gz', '_trim.fastq.gz')

                subprocess.call(["sbatch",
                                 "cutadapt-submit_edit",
                                 "-f", read1_fname,
                                 "-g", read2_fname,
                                 "-a", adapt_seq1,
                                 "-b", adapt_seq2,
                                 "-o", read1_outfname,
                                 "-p", read2_outfname])

                if lib not in lib_time_rep_pair_trimfname:
                    lib_time_rep_pair_trimfname[lib] = {}
                if time not in lib_time_rep_pair_trimfname[lib]:
                    lib_time_rep_pair_trimfname[lib][time] = {}
                if rep not in lib_time_rep_pair_trimfname[lib][time]:
                    lib_time_rep_pair_trimfname[lib][time][rep] = {}
                lib_time_rep_pair_trimfname[lib][time][rep]['R1'] = read1_outfname
                lib_time_rep_pair_trimfname[lib][time][rep]['R2'] = read2_outfname


### sort fastq file
if True:
    for lib in lib_time_rep_pair_fname:
        for time in lib_time_rep_pair_fname[lib]:
            for rep in lib_time_rep_pair_fname[lib][time]:
                pair_fname = lib_time_rep_pair_fname[lib][time][rep]
                read1_fname = pair_fname['R1'].replace('.fastq.gz', '_trim.fastq.gz')
                read2_fname = pair_fname['R2'].replace('.fastq.gz', '_trim.fastq.gz')
                out_fname = path + '_'.join([lib, time, rep])

                subprocess.call(["sbatch",
                                 "slide-submit_edit",
                                 "-f", read1_fname,
                                 "-g", read2_fname,
                                 "-x", ref_fname,
                                 "-m", sort_mode,
                                 "-o", out_fname])
