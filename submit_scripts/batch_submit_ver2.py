import subprocess
import glob
import sys, re

def read_ENCODE_table (fnames):
    cell_name_info = {}
    for fname in fnames:
        First = True
        for line in open(fname):
            cols = line.strip().split('\t')
            if First:
                First = False
                continue
            name, cell, source, data_list, control_list, peak_list = cols
            data_list = data_list.split(',')
            control_list = control_list.split(',')
            peak_list = peak_list.split(',')

            if cell not in cell_name_info:
                cell_name_info[cell] = {}
            assert name not in cell_name_info[cell]
            cell_name_info[cell][name] = {}
            cell_name_info[cell][name]['source'] = source
            cell_name_info[cell][name]['data'] = [data.strip() for data in data_list]
            cell_name_info[cell][name]['control'] = [control.strip() for control in control_list]
            cell_name_info[cell][name]['peak'] = [peak.strip() for peak in peak_list]
    return cell_name_info

def to_cn (fname, new_name, GC=False):
    f = open(new_name, 'w')
    include_GC = False
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            if cols[-1].startswith('GC') and GC:
                s = 'SNP\tChromosome\tPhysicalPosition\t' + "\t".join(cols[4:])
                include_GC = True
            else:
                s = 'SNP\tChromosome\tPhysicalPosition\t' + "\t".join(cols[4:-1])
            print >> f, s
            First = False
            continue
        s = "\t".join(cols[:2])
        s += "\t" + str(int(round(0.5*(float(cols[2]) + float(cols[3]))))) + "\t"
        if include_GC:
            s += "\t".join(cols[4:])
        else:
            s += "\t".join(cols[4:-1])
        print >> f, s
    f.close()

  
### basic parameters
#path = "/home/spark159/scratch/2022_11_24_H1_qc_reanlaysis/"
#path = "/home/spark159/scratch/2022_04_27_H1_protein_newqc/"
#path = "/home/spark159/scratch/2022_09_08_H1_HP1a_deep/"
#path = "/home/spark159/scratch/2022_12_13_H1_LKH_deep/"
#path = "/home/spark159/scratch/2022_10_28_H1_Ki67_deep/"
#path = "/home/spark159/scratch/2023_01_25_2nd_replicatesQC/"
#path = "/home/spark159/scratch/2023_02_07_sp_replicates_deep/"
path = "/home/spark159/scratch/2023_02_09_rep_fillin_qc/"


org_refname = {'human':path+"hg38", 'mouse':path+"mm10"}
cell_org = {'H1':'human', 'GM':'human', 'mCD8T':'mouse'}
cell_chrnames = {'H1':['chr%s' % (i) for i in range(1, 23)] + ['chrX', 'chrY'],
                 'GM':['chr%s' % (i) for i in range(1, 23)] + ['chrX'],
                 'mCD8T':['chr%s' % (i) for i in range(1, 20)] + ['chrX']}

#cell_chrnames = {'H1':['chr%s' % (i) for i in range(1, 6)]}

ENCODE_path = '/home/spark159/scratch/ENCODE_data/'
bs_fnames = [ENCODE_path+"ENCODE_BS_information.csv"]
chip_fnames = [ENCODE_path+"ENCODE_histone_chip_information.csv",
               ENCODE_path+"ENCODE_TF_chip_information.csv"]
cell_gtfname = {'GM':'ENCFF159KBI.gtf', 'H1':"ENCFF159KBI.gtf"}


### gather all fastq files
if False:
    exp_lane_pair_fname = {}
    for fname in glob.glob(path + '*.fastq.gz'):
        if fname.startswith(path+"Undetermined"):
            continue
        cols = fname.split('/')[-1].split('_')

        cell, sample, agent, tnum, rep = cols[:5]
        tnum = int(tnum)
        rep = int(rep[:-3])

        exp = (cell, sample, agent, tnum, rep)
        lane, pair = cols[-3], cols[-2]

        if exp not in exp_lane_pair_fname:
            exp_lane_pair_fname[exp] = {}
        if lane not in exp_lane_pair_fname[exp]:
            exp_lane_pair_fname[exp][lane] = {}

        exp_lane_pair_fname[exp][lane][pair] = fname

### Bowtie2 alignment
if False:
    for exp in exp_lane_pair_fname.keys():
        lane_pair_fname = exp_lane_pair_fname[exp]
        read_fnames = [[],[]]
        for lane in sorted(lane_pair_fname):
            pair_fname = lane_pair_fname[lane]
            pairs = sorted(pair_fname)
            for i in range(len(pairs)):
                pair = pairs[i]
                read_fnames[i].append(pair_fname[pair])
        read1_fnames, read2_fnames = read_fnames

        cell, sample, agent, tnum, rep = exp
        org = cell_org[cell]
        refname = org_refname[org]
        outfname = path + '_'.join([cell, sample, agent, str(tnum), str(rep)])
        
        subprocess.call(["sbatch",
                         "bowtie-submit",
                         "-f", ','.join(read1_fnames),
                         "-g", ','.join(read2_fnames),
                         "-x", refname,
                         "-o", outfname])


### gather all bam files
if True:
    act_tnums = {}
    for fname in glob.glob(path + '*.bam'):
        cols = fname.split('/')[-1].rsplit('.', 1)[0].split('_')
        cell, sample, agent, tnum, rep = cols

        tnum = int(tnum)
        rep = int(rep)

        act = (cell, sample, agent, rep)
        if act not in act_tnums:
            act_tnums[act] = []
        act_tnums[act].append(tnum)

    for act in act_tnums:
        act_tnums[act] = sorted(act_tnums[act])


### some basic QCs
if True:
    #bin_size = 1000
    #bin_size = 5000
    bin_size = 10000
    #bin_size = 25000

    for act in act_tnums:
        cell, sample, agent, rep = act
        tnums = act_tnums[act]

        # put input at last
        if tnums[0] == 0:
            tnums = tnums[1:] + [0]

        fnames = []
        for tnum in tnums:
            fname = path + '_'.join([cell, sample, agent, str(tnum), str(rep)]) +'.bam'
            fnames.append(fname)
            
        fnames = ','.join(fnames)
        org = cell_org[cell]
        refname = org_refname[org]
        
        # get read length distribution
        #outfname = path + '_'.join([cell, sample, agent])
        outfname = path
        subprocess.call(["sbatch",
                         "lendist_submit",
                         "-f", fnames,
                         "-o", outfname])

        # get read counts for each genomic bins
        outfname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', str(rep)])
        subprocess.call(["sbatch",
                         "bincount_submit",
                         "-f", fnames,
                         "-o", outfname,
                         "-x", refname+'.fa',
                         "-w", str(bin_size)])
                

### get condensability scores
if False:
    for cell in cell_sample_agent_tnum_fname:
        for sample in cell_sample_agent_tnum_fname[cell]:
            agent_tnum_fname = cell_sample_agent_tnum_fname[cell][sample]
            for agent in agent_tnum_fname:
                fname_list = []
                tnum_list = sorted(agent_tnum_fname[agent])[1:] + [0] 
                for tnum in tnum_list:
                    fname = agent_tnum_fname[agent][tnum]
                    fname_list.append(fname)
                fname_list = ','.join(fname_list)
                refname = org_refname[cell_org[cell]]
                
                chr_names = cell_chrnames[cell]
                for chr_name in chr_names:    
                    outfname = path + '_'.join([cell, sample, agent, chr_name])
                    # get condensabiltiy score
                    subprocess.call(["sbatch", "condense_seq-submit_edit", "-f", fname_list, "-o", outfname, "-x", refname+'.fa', "-c", chr_name])



### check nucleosome positioning signal (QC)
if False:
    for cell in cell_sample_agent_tnum_fname:
        for sample in cell_sample_agent_tnum_fname[cell]:
            agent_tnum_fname = cell_sample_agent_tnum_fname[cell][sample]
            for agent in agent_tnum_fname:
                for chr_name in chr_names:
                    fname = path + '_'.join([cell, sample, agent, chr_name, 'peak.cn'])
                    outfname = path + '_'.join([cell, sample, agent, chr_name])
 
                    # get nucleosome motif sequence
                    subprocess.call(["sbatch", "motif-submit", "-f", fname, "-x", refname+'.fa', "-c", chr_name, "-o", outfname])




### write annotation file
if False:
    # read BS data
    cellname_bs_info = read_ENCODE_table(bs_fnames)
    # read chip seq data
    cellname_chip_info = read_ENCODE_table(chip_fnames)
    
    for cell in cell_sample_agent_tnum_fname:
        if cell.startswith('GM'):
            cellname = 'GM12878'
        elif cell.startswith('H1'):
            cellname = 'H1-hESC'
        elif cell.startswith('mCD8T'):
            cellname = 'Mouse CD8 T cell (invitro activated)'

        bs_info = cellname_bs_info[cellname]
            
        bs_fname = {}
        for bs in bs_info:
            bed_files = bs_info[bs]['peak']
            assert len(bed_files) == 1
            bs_fname[bs] = bed_files[0] + '.bed'

        bs_files = []
        for bs, fname in bs_fname.items():
            bs_files.append(bs)
            bs_files.append(ENCODE_path+fname)
        bs_files = ','.join(bs_files)

        chip_info = cellname_chip_info[cellname]

        chip_fname = {}
        for chip in chip_info:
            bed_files = chip_info[chip]['peak']
            assert len(bed_files) == 1
            chip_fname[chip] = bed_files[0] + '.bed'

        chip_files = []
        for chip, fname in chip_fname.items():
            chip_files.append(chip)
            chip_files.append(ENCODE_path+fname)
        chip_files = ','.join(chip_files)
        
        for sample in cell_sample_agent_tnum_fname[cell]:
            agent_tnum_fname = cell_sample_agent_tnum_fname[cell][sample]
            for agent in agent_tnum_fname:
                chr_names = cell_chrnames[cell]
                for chr_name in chr_names:
                    refname = org_refname[cell_org[cell]]
                    fname = path + '_'.join([cell, sample, agent, chr_name])
                    outfname = path + '_'.join([cell, sample, agent, chr_name])
                    subprocess.call(["sbatch", "annot_submit_edit", "-f", fname, "-x", refname+'.fa', "-n", chr_name, "-o", outfname, "-b", bs_files, "-c", chip_files])
                    #subprocess.call(["sbatch", "annot_submit_edit", "-f", fname, "-x", refname+'.fa', "-n", chr_name, "-o", outfname, "-c", chip_files])
            

# make profile
if False:
    for cell in cell_sample_agent_tnum_fname:
        gtfname = path + cell_gtfname[cell]
        for sample in cell_sample_agent_tnum_fname[cell]:
            agent_tnum_fname = cell_sample_agent_tnum_fname[cell][sample]
            for agent in agent_tnum_fname:
                for chr_name in chr_names:
                    # make profile for all annotations
                    fname = path + '_'.join([cell, sample, agent, chr_name, '167win25step', 'anot']) + '.cn'
                    outfname = path + '_'.join([cell, sample, agent, chr_name])
                    subprocess.call(["sbatch", "profile_submit", "-f", fname, "-g", gtfname, "-x", refname+'.fa', "-c", chr_name, "-o", outfname])

                    # make profile for occupancy
                    fname = path + '_'.join([cell, sample, agent, chr_name, 'occ']) + '.cn'
                    outfname = path + '_'.join([cell, sample, agent, chr_name, 'occ'])
                    subprocess.call(["sbatch", "profile_submit", "-f", fname, "-g", gtfname, "-x", refname+'.fa', "-c", chr_name, "-o", outfname])
                    
