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
## Human cell lines
#path = "/home/spark159/scratch/2022_11_23_some_replicates_QC/"
#refname = path + "hg38"
#chr_names = ['chr1']
#chr_names += ['chr%s' % (i) for i in range(2, 23)]
#chr_names += ['chrX', 'chrY']
#ENCODE_path = '/home/spark159/scratch4-tha4/sangwoo/ENCODE_data/'
#bs_fname = ENCODE_path+"ENCODE_BS_information.csv"
#chip_fname = ENCODE_path+"ENCODE_histone_chip_information.csv"
#cell_gtfname = {'GM':'ENCFF159KBI.gtf', 'H1':"ENCFF159KBI.gtf"}

## Mouse cell lines
#path = "/home/spark159/scratch/2022_11_08_mCD8T_KO_detail/"
#refname = path + "mm10"
#chr_names = ['chr1']
#chr_names = ['chr2']
#chr_names += ['chr%s' % (i) for i in range(2, 20)]
#chr_names = ['chr%s' % (i) for i in range(3, 20)]
#chr_names = ['chr%s' % (i) for i in range(2, 10)]
#chr_names += ['chrX', 'chrY']
#ENCODE_path = '/home/spark159/scratch/MouseCD8TcellData/'
#chip_fname = ENCODE_path + 'MouseCD8Tcell_histone_chip_seq.csv'
#cell_gtfname = {'mCD8T':'gencodeM21pri-UCSC-tRNAs-ERCC-phiX.gtf'}

### parameters
#path = "/home/spark159/scratch/2022_11_24_H1_qc_reanlaysis/"
#path = "/home/spark159/scratch/2022_04_27_H1_protein_newqc/"
#path = "/home/spark159/scratch/2022_09_08_H1_HP1a_deep/"
#path = "/home/spark159/scratch/2022_12_13_H1_LKH_deep/"
#path = "/home/spark159/scratch/2022_10_28_H1_Ki67_deep/"
#path = "/home/spark159/scratch/2023_01_25_2nd_replicatesQC/"
path = "/home/spark159/scratch/2023_02_07_sp_replicates_deep/"
#path = "/home/spark159/scratch/2023_02_09_rep_fillin_qc/"
#path = "/home/spark159/scratch/2023_02_11_mCD8T_KO_tnum2_fillin/"

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
if True:
    cell_sample_agent_tnum_pair_lane_fname = {}
    for fname in glob.glob(path + '*.fastq.gz'):
        if fname.startswith(path+"Undetermined"):
            continue
        cols = fname.split('/')[-1].split('_')

        if 're' in cols:
            cell, sample = cols[:2]
            sample = 're:' + sample
            tnum = int(cols[-6])
            agent = cols[-7]
            pair = cols[-2]
            lane = cols[-3]
        else:
            cell, sample = cols[:2]
            tnum = int(cols[-5])
            agent = cols[-6]
            pair = cols[-2]
            lane = cols[-3]

        #if fname.endswith("trimmed.fastq.gz"):
        #    continue
        #cols = re.split('_|-', fname.split('/')[-1])
        #cell, sample = cols[:2]
        #tnum = int(cols[-3])
        #agent = cols[-4]
        #pair = cols[-2]
        #lane = 'L001'

        if cell not in cell_sample_agent_tnum_pair_lane_fname:
            cell_sample_agent_tnum_pair_lane_fname[cell] = {}
        if sample not in cell_sample_agent_tnum_pair_lane_fname[cell]:
            cell_sample_agent_tnum_pair_lane_fname[cell][sample] = {}
        if agent not in cell_sample_agent_tnum_pair_lane_fname[cell][sample]:
            cell_sample_agent_tnum_pair_lane_fname[cell][sample][agent] = {}
        if tnum not in cell_sample_agent_tnum_pair_lane_fname[cell][sample][agent]:
            cell_sample_agent_tnum_pair_lane_fname[cell][sample][agent][tnum] = {}
        if pair not in cell_sample_agent_tnum_pair_lane_fname[cell][sample][agent][tnum]:
            cell_sample_agent_tnum_pair_lane_fname[cell][sample][agent][tnum][pair] = {}

        assert lane not in cell_sample_agent_tnum_pair_lane_fname[cell][sample][agent][tnum][pair]
        cell_sample_agent_tnum_pair_lane_fname[cell][sample][agent][tnum][pair][lane] = fname


### Bowtie2 alignment
if True:
    for cell in cell_sample_agent_tnum_pair_lane_fname:
        for sample in cell_sample_agent_tnum_pair_lane_fname[cell]:
            for agent in cell_sample_agent_tnum_pair_lane_fname[cell][sample]:
                tnum_pair_lane_fname = cell_sample_agent_tnum_pair_lane_fname[cell][sample][agent]
                for tnum in sorted(tnum_pair_lane_fname.keys()):
                    lanes1 = sorted(tnum_pair_lane_fname[tnum]['R1'].keys())
                    lanes2 = sorted(tnum_pair_lane_fname[tnum]['R2'].keys())
                    read1_fnames = [tnum_pair_lane_fname[tnum]['R1'][lane] for lane in lanes1]
                    read2_fnames = [tnum_pair_lane_fname[tnum]['R2'][lane] for lane in lanes2]

                    refname = org_refname[cell_org[cell]]
                    outfname = path + '-'.join([cell, sample, agent, str(tnum)])
                                        
                    # Bowtie2 alignment
                    if len(lanes1) == len(lanes2) == 1:
                        assert lanes1[0] == lanes2[0] == 'L000'
                        #subprocess.call(["sbatch", "bowtie-submit-trim", "-f", ','.join(read1_fnames), "-g", ','.join(read2_fnames), "-x", refname, "-o", outfname])
                    else:
                        subprocess.call(["sbatch", "bowtie-submit", "-f", ','.join(read1_fnames), "-g", ','.join(read2_fnames), "-x", refname, "-o", outfname])


### gather all bam files
if False:
    cell_sample_agent_tnum_fname = {}
    for fname in glob.glob(path + '*.bam'):
        cols = fname.split('/')[-1].rsplit('.', 1)[0].split('-')        
        try:
            cell, sample, agent, tnum = cols
        except:
            try:
                cell, option, sample, agent, tnum = cols
                sample = '-'.join([option, sample])
            except:
                continue

        # need to fix
        if sample.startswith('re'):
            continue

        try:
            tnum = int(tnum)
        except:
            continue

        if cell not in cell_sample_agent_tnum_fname:
            cell_sample_agent_tnum_fname[cell] = {}
        if sample not in cell_sample_agent_tnum_fname[cell]:
            cell_sample_agent_tnum_fname[cell][sample] = {}
        if agent not in cell_sample_agent_tnum_fname[cell][sample]:
            cell_sample_agent_tnum_fname[cell][sample][agent] = {}

        assert tnum not in cell_sample_agent_tnum_fname[cell][sample][agent]
        cell_sample_agent_tnum_fname[cell][sample][agent][tnum] = fname


### some basic QCs
if False:
    #bin_size = 5000
    #bin_size = 1000
    bin_size = 10000
    #bin_size = 25000
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
                
                # get read length distribution
                outfname = path + '_'.join([cell, sample, agent])
                subprocess.call(["sbatch", "lendist_submit", "-f", fname_list, "-o", outfname])

                # get read counts for each genomic bins
                outfname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb'])
                subprocess.call(["sbatch", "bincount_submit", "-f", fname_list, "-o", outfname, "-x", refname+'.fa', "-w", str(bin_size)])

                # get profile file for QC
                if False:
                    #oldcn = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb']) + '_bin.cn'
                    #newcn = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb']) + '.cn'
                    #outfname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb'])
                    #to_cn(oldcn, newcn)

                    #subprocess.call(["sbatch", "binannot_submit", "-f", newcn, "-x", refname+'.fa', "-w", str(bin_size), "-s", str(bin_size), "-o", outfname])

                    fname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb', 'anot']) + '.cn'
                    gtfname = path + cell_gtfname[cell]
                    outfname = path + '_'.join([cell, sample, agent, str(int(bin_size/1000.0)) + 'kb'])
                    #subprocess.call(["sbatch", "profile_submit", "-f", fname, "-g", gtfname, "-x", refname+'.fa', "-c", "chr1", "-o", outfname])

                    chr_names = cell_chrnames[cell]
                    for chr_name in chr_names:
                        outfname = path + '_'.join([cell, sample, agent, chr_name])

                        # get nucleosome occupancy along the genome
                        subprocess.call(["sbatch", "occupancy-submit", "-f", fname_list, "-o", outfname, "-x", refname+'.fa', "-c", chr_name])
                

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
                    
