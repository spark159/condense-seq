import subprocess
import glob
import sys, re

def read_ENCODE_table (fname):
    cell_name_info = {}
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


"""
index_num = {"TTACCGAC~CGAATACG":1,"TCGTCTGA~GTCCTTGA":2,"TTCCAGGT~CAGTGCTT":3,"TACGGTCT~TCCATTGC":4,"AAGACCGT~GTCGATTG":5,"CAGGTTCA~ATAACGCC":6,"TAGGAGCT~GCCTTAAC":7,"TACTCCAG~GGTATAGG":8,"AGTGACCT~TCTAGGAG":9,"AGCCTATC~TGCGTAAC":10,"TCATCTCC~CTTGCTAG":11,"CCAGTATC~AGCGAGAT":12,"TTGCGAGA~TATGGCAC":13,"GAACGAAG~GAATCACC":14,"CGAATTGC~GTAAGGTG":15,"GGAAGAGA~CGAGAGAA":16,"TCGGATTC~CGCAACTA":17,"CTGTACCA~CACAGACT":18,"GAGAGTAC~TGGAAGCA":19,"TCTACGCA~CAATAGCC":20,"GCAATTCC~CTCGAACA":21,"CTCAGAAG~GGCAAGTT":22,"GTCCTAAG~AGCTACCA":23,"GCGTTAGA~CAGCATAC":24,"CGTATCTC~CAAGGTAC":25,"TTACGTGC~AGACCTTG":26,"AGCTAAGC~GTCGTTAC":27,"AAGACACC~GTAACCGA":28,"CAACTCCA~GAATCCGT":29,"GATCTTGC~CATGAGCA":30,"CTTCACTG~CTTAGGAC":31,"CTCGACTT~ATCTGACC":32}


for index in index_num:
    num = index_num[index]
    files = [[], []]
    for pair in range(2):
        for lane_num in range(2):
            name = "data/sp_spd_tests/" + "H72M2BCX3_" + str(lane_num+1) + '_' + index + '_' + str(pair+1) + '.fastq.gz'
            #subprocess.call(["fastqc", name])
            #subprocess.call(["mv", "data/" + "H72M2BCX3_" + str(lane_num+1) + '_' + index + '_' + str(pair+1) + '_fastqc.html', "data/" + "sp_spd_test" + str(num) +  '_' + str(lane_num+1) + '_' + str(pair+1) + '_fastqc.html'])
            files[pair].append(name)
    # Bowtie2 alignment
    files1 = ','.join(files[0])
    files2 = ','.join(files[1])
    file_list = files1 + ',' + files2
    outfname = "work/sp_spd_tests/" + "sp_spd_test" + str(num)
    #subprocess.call(["sbatch", "bowtie-submit", "-f", files1, "-g", files2, "-o", outfname])
    subprocess.call(["sbatch", "qc_submit", "-f", file_list, "-o", outfname])
    #break


for i in range(4):
    file_list = []
    outfname = ""
    if i % 2 == 0:
        outfname += "NCP_"
    else:
        outfname += "DNA_"
    if i < 2:
        outfname += "Spermine(4+)"
    else:
        outfname += "Spermidine(3+)"
    for j in range(8*i, 8*(i+1)):
        file_list.append("work/sp_spd_tests/" + "sp_spd_test" + str(j+1) +'.bam')
    file_list = ','.join(file_list)
    outfname = "data/sp_spd_tests/" + outfname
    subprocess.call(["sbatch", "bincount_submit", "-f", file_list, "-o", outfname])
"""



"""
def write_temp_para (fname, energy_list, valency_list):
    f = open(fname + '_para.cn', 'w')
    print >> f, "%s\t%s\t%s\t%s\t%s" % ("SNP", "Chromosome", "PhysicalPosition", "Score", "Valency")
    for i in range(len(energy_list)):
        print >> f, "%s\t%s\t%s\t%s\t%s" % (str(i), "chr1", str(10*i), energy_list[i], valency_list[i])
    f.close()

# almost pure case
N = 100
f = 23
E = -10
for v in [9]:     
    f_single = 3 + 10*v
    for i in range(10)[::-1]:
        E_single = - 2*i
        energy_list = [E_single] + [E]*(N-1)
        valency_list = [f_single] + [f]*(N-1)
        fname = 'single_' + str(N) + '_' + str(f) + '_' + str(E) + ':' + str(f_single) + '_' + str(E_single)
        write_temp_para (fname, energy_list, valency_list)
        subprocess.call(["sbatch", "rgs-submit", "-f", fname + '_para.cn', "-o", fname])
"""

#path = "work/2021_02_10_H1_sp_spd_test/"
#adapt_seq = 'AGATCGGAAGAGCACACGTC'
#refname = path + 'hg19'
#cells = ['H1']
#stypes = ['DNA', 'NCP-new']
#agents = ['sp', 'spd']

#for cell in cells:
#    for stype in stypes:
#        for agent in agents:
#            for i in range(6):
#                raw_name_pair = []
#                #trimmed_name_pair = []
#                for j in [1,2]:
#                    raw_name = path+cell+'-'+stype+'-'+agent+'-'+str(i)+'_R'+str(j)+'_001.fastq.gz'
#                    raw_name_pair.append(raw_name)

#                    #trimmed_name = path+cell+'-'+stype+'-'+agent+'-'+str(i)+'_R'+str(j)+'_001_trimmed.fastq.gz'
#                    #trimmed_name_pair.append(trimmed_name)

#                # read QC
#                #subprocess.call(["sbatch", "qc_submit", "-f", ','.join(file_pair), "-o", outfname])

#                ## adaptor trimming
#                ##subprocess.call(["sbatch", "cutadapt-submit", "-f", raw_name_pair[0], "-g", raw_name_pair[1], "-a", adapt_seq, "-o", trimmed_name_pair[0], "-p", trimmed_name_pair[1]])

#                # Bowtie2 alignment
#                #outfname = path+cell+'-'+stype+'-'+agent+'-'+str(i)
#                #outfname = path+cell+'-'+stype+'-'+agent+'-'+str(i)
#                #subprocess.call(["sbatch", "bowtie-submit", "-f", raw_name_pair[0], "-g", raw_name_pair[1], "-x", refname, "-o", outfname])


## ATAC-seq binning go through all chromosomes
#chr_names = ['chr' + str(i+1) for i in range(22)]
#chr_names += ['chrX', 'chrY']
#chr_names = ['chr' + str(i+1) for i in range(3)]
#chr_names = ['chr12']
#chr_names = ['chr' + str(i) for i in range(1, 12) + range(13, 23)]
#chr_names = ['chr' + str(i) for i in range(1, 23)]
#chr_names += ['chrX', 'chrY']
#chr_names = ['chr13']

#for chr_name in chr_names:
#    subprocess.call(["sbatch", "condense_seq-submit", "-c", chr_name])

#sys.exit(1)

# H1 spermine4+ go through all chromosomes
#chr_names = ['chr' + str(i+1) for i in range(22)]
#chr_names += ['chrX', 'chrY']
#chr_names = ['chr' + str(i+1) for i in range(3)]
#chr_names = ['chr12']
#chr_names = ['chr' + str(i) for i in range(1, 12) + range(13, 23)]
#chr_names = ['chr' + str(i) for i in range(1, 23)]
#chr_names += ['chrX', 'chrY']
#chr_names = ['chr13']

#for chr_name in chr_names:
#    subprocess.call(["sbatch", "condense_seq-submit", "-c", chr_name])
    
#sys.exit(1)

# Miseq protein qc samples
#path = "work/2022_04_04_H1_protein_qc_again/"
#refname = path + 'hg38'
#sample_num = 12

#for i in range(1, sample_num+1):
#    read1_fname = path + '%d_S%d_L001_R1_001.fastq.gz' % (i, i)
#    read2_fname = path + '%d_S%d_L001_R2_001.fastq.gz' % (i, i)
#    outfname = path + '%d' % (i)
#
#    # Bowtie2 alignment
#    subprocess.call(["sbatch", "bowtie-submit", "-f", read1_fname, "-g", read2_fname, "-x", refname, "-o", outfname])

#sys.exit(1)

## check bin counting
#win_size = 100000
#for i in range(1, sample_num+1):
#    fname = path + '%d.bam' % (i)
#    outfname = path + '%d_%dkb' % (i, win_size/1000)
#    subprocess.call(["sbatch", "bincount_submit", "-f", fname + ',' + fname, "-o", outfname, "-x", refname+'.fa', "-w", str(win_size)])

#sys.exit(1)


## check length
#for i in range(1, sample_num+1):
#    fname = path + '%d.bam' % (i)
#    outfname = path + '%d' % (i)
#    subprocess.call(["sbatch", "lendist_submit", "-f", fname, "-o", outfname])

#sys.exit(1) 



"""
#inpath = "work/2021_06_07_H1_sp_detail/"
#outpath = "work/2021_06_07_H1_sp_detail/"
#inpath = "work/2022_04_27_H1_protein_filling/"
#outpath = "work/2022_04_27_H1_protein_filling/"
inpath = "/home/spark159/scratch4-tha4/sangwoo/2022_07_04_mouseCD8/"
outpath = "/home/spark159/scratch4-tha4/sangwoo/2022_07_04_mouseCD8/"
refname = outpath + 'mm10'

cell_sample_agent_tnum_pair_lane_fname = {}
for fname in glob.glob(inpath + '*.fastq.gz'):
    if fname.startswith(inpath+"Undetermined"):
        continue
    cols = fname.split('/')[-1].split('_')
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

for cell in cell_sample_agent_tnum_pair_lane_fname:
    for sample in cell_sample_agent_tnum_pair_lane_fname[cell]:
        for agent in cell_sample_agent_tnum_pair_lane_fname[cell][sample]:
            tnum_pair_lane_fname = cell_sample_agent_tnum_pair_lane_fname[cell][sample][agent]
            for tnum in sorted(tnum_pair_lane_fname.keys()):
                lanes1 = sorted(tnum_pair_lane_fname[tnum]['R1'].keys())
                lanes2 = sorted(tnum_pair_lane_fname[tnum]['R2'].keys())
                read1_fnames = [tnum_pair_lane_fname[tnum]['R1'][lane] for lane in lanes1]
                read2_fnames = [tnum_pair_lane_fname[tnum]['R2'][lane] for lane in lanes2]

                # Bowtie2 alignment
                outfname = outpath + '-'.join([cell, sample, agent, str(tnum)])
                subprocess.call(["sbatch", "bowtie-submit", "-f", ','.join(read1_fnames), "-g", ','.join(read2_fnames), "-x", refname, "-o", outfname])

sys.exit(1) 
"""

# gather all bam files
#path = "work/2021_05_13_H1_test/"
#refname = path + 'hg19'
#path = "work/2022_04_27_H1_protein_filling/"
#path = "work/2022_03_15_H1_protein_test/"
path = "/home/spark159/scratch4-tha4/sangwoo/2022_07_04_GM_NCP/"
#path = "/home/spark159/scratch4-tha4/sangwoo/2022_07_04_mouseCD8/"
refname = path + 'hg38'
#refname = path + "mm10"

cell_sample_agent_tnum_fname = {}
for fname in glob.glob(path + '*.bam'):
    cols = fname.split('/')[-1].rsplit('.', 1)[0].split('-')
    try:
        cell, sample, agent, tnum = cols
    except:
        cell, option, sample, agent, tnum = cols
        sample = '-'.join([option, sample])
        
    tnum = int(tnum)

    if cell not in cell_sample_agent_tnum_fname:
        cell_sample_agent_tnum_fname[cell] = {}
    if sample not in cell_sample_agent_tnum_fname[cell]:
        cell_sample_agent_tnum_fname[cell][sample] = {}
    if agent not in cell_sample_agent_tnum_fname[cell][sample]:
        cell_sample_agent_tnum_fname[cell][sample][agent] = {}
        
    assert tnum not in cell_sample_agent_tnum_fname[cell][sample][agent]
    cell_sample_agent_tnum_fname[cell][sample][agent][tnum] = fname


# check length dist
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
            outfname = path + '_'.join([cell, sample, agent])
            #subprocess.call(["sbatch", "lendist_submit", "-f", fname_list, "-o", outfname])


# bincount
#win_size = 10000
win_size = 1000
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
            outfname = path + '_'.join([cell, sample, agent, str(int(win_size/1000.0)) + 'kb'])
            #subprocess.call(["sbatch", "bincount_submit", "-f", fname_list, "-o", outfname, "-x", refname+'.fa', "-w", str(win_size)])
            #subprocess.call([". ", "bincount_submit", "-f", fname_list, "-o", outfname, "-x", refname+'.fa', "-w", str(win_size)])

#sys.exit(1)


# make bin annotation file
cellname = "GM12878"
win_size = 1000
bin_width = win_size
bin_step = win_size

# read BS data
bs_path = '/home/spark159/scratch4-tha4/sangwoo/ENCODE_data/'
bs_info = read_ENCODE_table(bs_path+"ENCODE_BS_information.csv")[cellname]

bs_fname = {}
for bs in bs_info:
    bed_files = bs_info[bs]['peak']
    assert len(bed_files) == 1
    bs_fname[bs] = bed_files[0] + '.bed'

bs_files = []
for bs, fname in bs_fname.items():
    bs_files.append(bs)
    bs_files.append(bs_path+fname)
bs_files = ','.join(bs_files)    


# read chip seq data
chip_path = '/home/spark159/scratch4-tha4/sangwoo/ENCODE_data/'
chip_info = read_ENCODE_table(chip_path+"ENCODE_histone_chip_information.csv")[cellname]

chip_fname = {}
for chip in chip_info:
    bed_files = chip_info[chip]['peak']
    assert len(bed_files) == 1
    chip_fname[chip] = bed_files[0] + '.bed'

chip_files = []
for chip, fname in chip_fname.items():
    chip_files.append(chip)
    chip_files.append(chip_path+fname)
chip_files = ','.join(chip_files)    


#chip_fname = {"H2AZ":"H2AZ.bedGraph",
#              "H3k27ac":"H3k27ac.bedGraph",
#              "H3k27me3":"H3k27me3.bedGraph",
#              "H3k36me3":"H3k36me3.bedGraph",
#              "H3k4me1":"H3k4me1.bedGraph",
#              "H3k4me3":"H3k4me3.bedGraph",
#              "H3k9ac":"H3k9ac.bedGraph",
#              "H3k9me3":"H3k9me3.bedGraph"}

#chip_files = []
#for chip, fname in chip_fname.items():
#    chip_files.append(chip)
#    chip_files.append(path+fname)

#chip_files = ','.join(chip_files)

for cell in cell_sample_agent_tnum_fname:
    for sample in cell_sample_agent_tnum_fname[cell]:
        agent_tnum_fname = cell_sample_agent_tnum_fname[cell][sample]
        for agent in agent_tnum_fname:
            # convert bin.cn to .cn
            oldcn = path + '_'.join([cell, sample, agent, str(int(win_size/1000.0)) + 'kb']) + '_bin.cn'
            newcn = path + '_'.join([cell, sample, agent, str(int(win_size/1000.0)) + 'kb']) + '.cn'
            #to_cn(oldcn, newcn)
            # write annotation file
            outfname = path + '_'.join([cell, sample, agent, str(int(win_size/1000.0)) + 'kb'])
            #subprocess.call(["sbatch", "binannot_submit", "-f", newcn, "-x", refname+'.fa', "-w", str(bin_width), "-s", str(bin_step), "-o", outfname, "-b", bs_files, "-c", chip_files])
            #subprocess.call(["sbatch", "binannot_submit", "-f", newcn, "-x", refname+'.fa', "-w", str(bin_width), "-s", str(bin_step), "-o", outfname, "-c", chip_files])
            # exclude other ENCODE annotation
            #subprocess.call(["sbatch", "binannot_submit", "-f", newcn, "-x", refname+'.fa', "-w", str(bin_width), "-s", str(bin_step), "-o", outfname])

#sys.exit(1)


# make profile
gtfname =  path + "ENCFF159KBI.gtf"
#gtfname = "work/condense_seq/Homo_sapiens.GRCh37.87.gtf"
win_size = 1000
for cell in cell_sample_agent_tnum_fname:
    for sample in cell_sample_agent_tnum_fname[cell]:
        agent_tnum_fname = cell_sample_agent_tnum_fname[cell][sample]
        for agent in agent_tnum_fname:
            fname = path + '_'.join([cell, sample, agent, str(int(win_size/1000.0)) + 'kb', 'anot']) + '.cn'
            outfname = path + '_'.join([cell, sample, agent, str(int(win_size/1000.0)) + 'kb', 'TSS-TTS'])
            subprocess.call(["sbatch", "profile_submit", "-f", fname, "-g", gtfname, "-x", refname+'.fa', "-o", outfname])
            
                
sys.exit(1)

## check length
#win_size = 1000
#for cell in cells:
#    for stype in stypes:
#        for agent in agents:
#            fname_list = []
#            for i in range(6):
#                fname = path+cell+'-'+stype+'-'+agent+'-'+str(i)+'.bam'
#                fname_list.append(fname)
#            fname_list = ','.join(fname_list)
#            outfname = path + '_'.join([cell, stype, agent, str(int(win_size/1000.0)) + 'kb'])
#            subprocess.call(["sbatch", "bincount_submit", "-f", fname_list, "-o", outfname, "-x", refname+'.fa', "-w", str(win_size)])


## check length dist
#for cell in cells:
#    for stype in stypes:
#        for agent in agents:
#            fname_list = []
#            for i in range(6):
#                fname = path+cell+'-'+stype+'-'+agent+'-'+str(i)+'.bam'
#                fname_list.append(fname)
#            fname_list = ','.join(fname_list)
#            outfname = path + '_'.join([cell, stype, agent])
#            subprocess.call(["sbatch", "lendist_submit", "-f", fname_list, "-o", outfname])


## bin count
#win_size = 1000
#for cell in cells:
#    for stype in stypes:
#        for agent in agents:
#            fname_list = []
#            for i in range(6):
#                fname = path+cell+'-'+stype+'-'+agent+'-'+str(i)+'.bam'
#                fname_list.append(fname)
#            fname_list = ','.join(fname_list)
#            outfname = path + '_'.join([cell, stype, agent, str(int(win_size/1000.0)) + 'kb'])
#            subprocess.call(["sbatch", "bincount_submit", "-f", fname_list, "-o", outfname, "-x", refname+'.fa', "-w", str(win_size)])

## make bin annotation file
#win_size = 1000
#bin_width = win_size
#bin_step = win_size
#
#chip_fname = {"H2AZ":"H2AZ.bedGraph",
#              "H3k27ac":"H3k27ac.bedGraph",
#              "H3k27me3":"H3k27me3.bedGraph",
#              "H3k36me3":"H3k36me3.bedGraph",
#              "H3k4me1":"H3k4me1.bedGraph",
#              "H3k4me3":"H3k4me3.bedGraph",
#              "H3k9ac":"H3k9ac.bedGraph",
#              "H3k9me3":"H3k9me3.bedGraph"}


#chip_names, chip_files = [], []
#for chip, fname in chip_fname.items():
#    chip_names.append(chip)
#    chip_files.append(path+fname)

#chip_names = ','.join(chip_names)
#chip_files = ','.join(chip_files)

#cells = ['H1']
#stypes = ['DNA']
#agents = ['sp']


#for cell in cells:
#    for stype in stypes:
#        for agent in agents:
#            fname = path + '_'.join([cell, stype, agent, str(int(win_size/1000.0)) + 'kb']) + '.cn'
#            subprocess.call(["./binannot_submit", "-f", fname, "-x", refname+'.fa', "-w", str(bin_width), "-s", str(bin_step), "-n", chip_names, "-c", chip_files], shell=True)

