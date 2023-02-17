import sys
import subprocess
import glob

def read_data_info (fname):
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

def make_download_file (fname, data_list, extension='bam'):
    f = open(fname, 'w')
    for data in data_list:
        print >> f, "https://www.encodeproject.org/files/%s/@@download/%s.%s" % (data, data, extension)
    f.close()
    return

# load data information
chipname_info = read_data_info("MouseCD8Tcell_histone_chip_seq.csv")["Mouse CD8 T cell (invitro activated)"]

# download all fastq files
fastq_fnames = set([])
for chipname in chipname_info:
    fastq_fnames |= set(chipname_info[chipname]['data'])
    fastq_fnames |= set(chipname_info[chipname]['control'])

fastq_fnames = list(fastq_fnames)
for fastq_fname in fastq_fnames:
    if glob.glob(fastq_fname + '*.fastq'):
        continue
    subprocess.call("fastq-dump --split-3" + fastq_fname, shell=True)
    


# Bowtie2 alignment to make bam files
#path = './'
refname = 'mm10'
for fastq_fname in fastq_fnames:
    if glob.glob(fastq_fname + '.bam'):
        continue
    else:
        print fastq_fname
        if len(glob.glob(fastq_fname+'*.fastq')) == 1:
            #continue
            #print fastq_fname
            subprocess.call(["sbatch", "/home/spark159/submit_scripts/bowtie-submit_single", "-f", fastq_fname+'.fastq', "-x", refname, "-o", fastq_fname])
        elif len(glob.glob(fastq_fname+'*.fastq')) == 2:
            #print fastq_fname
            subprocess.call(["sbatch", "/home/spark159/submit_scripts/bowtie-submit", "-f", fastq_fname+'_1.fastq', "-g", fastq_fname+'_2.fastq', "-x", refname, "-o", fastq_fname])
        else:
            assert False

#sys.exit(1)

# make ChromHMM input table
cell_name = "Mouse CD8 T cell (invitro activated)"
chipnames = ["H3K4me1",
             "H3K4me3",
             "H3K27ac",
             "H3K27me3",
             "H3K36me3",
             "H3K9me3"]

path = ""
fname = cell_name + "_input_table.txt"
f = open(fname, "w")
doneset = set([])
for chipname in chipnames:
    # merge replicates
    data_list = sorted(chipname_info[chipname]["data"])
    combined_data  = "-".join(data_list)
    data_input = [path + data + '.bam' for data in data_list]
    data_output = path + combined_data + '.bam'
    if len(glob.glob(path + data_output)) <= 0 and data_output not in doneset:
    #if data_output not in doneset:
        merge_cmd = ["samtools", "merge", data_output] + data_input
        #subprocess.call(merge_cmd, shell=True)
        #subprocess.call(["sbatch", "/home/spark159/submit_scripts/samtools-submit", '-o', data_output, '-f', ','.join(data_input)])
        doneset.add(data_output)
    
    control_list = sorted(chipname_info[chipname]["control"])
    combined_control = "-".join(control_list)
    control_input = [path + control + '.bam' for control in control_list]
    control_output = path + combined_control + '.bam'
    if len(glob.glob(path + control_output)) <= 0 and control_output not in doneset:
    #if control_output not in doneset:
        merge_cmd = ["samtools", "merge", control_output] + control_input
        #subprocess.call(merge_cmd, shell=True)
        #subprocess.call(["sbatch", "/home/spark159/submit_scripts/samtools-submit", '-o', control_output, '-f', ','.join(control_input)])
        doneset.add(control_output)

    # make input table
    print >> f, "%s\t%s\t%s\t%s" % (cell_name, chipname, combined_data + '.bam', combined_control + '.bam')
    #print >> f, "%s\t%s\t%s\t%s" % (cell_name, chipname, data_list[0] + '.bam', control_list[0] + '.bam')

#f.close()




## make batch submit file (load bam files)
#cell_name = "GM12878"
#submit_fname = cell_name + "_chip_bam_files.txt"
#extension = 'bam'
#data_list = set([])
#for chipname in cell_chipname_info[cell_name]:
#    data_list |= set(cell_chipname_info[cell_name][chipname]['data'])
#    data_list |= set(cell_chipname_info[cell_name][chipname]['control'])
#data_list = list(data_list)
#make_download_file (submit_fname, data_list, extension=extension)


"""
# make batch submit file (load bed files)
cell_name = "GM12878"
submit_fname = cell_name + "_chip_bed_files.txt"
extension = 'bed.gz'
data_list = set([])
for chipname in cell_chipname_info[cell_name]:
    data_list |= set(cell_chipname_info[cell_name][chipname]['peak'])
data_list = list(data_list)
make_download_file (submit_fname, data_list, extension=extension)


# download ENCODE data
subprocess.call("xargs -n 1 curl -O -L < " + submit_fname, shell=True)
"""

"""
# make ChromHMM input table
cell_name = "GM12878"
chipnames = ["CTCF",
             "H3K27ac",
             "H3K27me3",
             "H3K36me3",
             "H3K4me1",
             "H3K4me2",
             "H3K4me3",
             "H3K9ac",
             "H3K9me3",
             "H4K20me1"]

path = "GM_chip_seq_files/"
fname = cell_name + "_input_table.txt"
f = open(fname, "w")

for chipname in chipnames:
    # merge replicates
    data_list = sorted(cell_chipname_info[cell_name][chipname]["data"])
    combined_data  = "-".join(data_list)
    data_input = [path + data + '.bam' for data in data_list]
    data_output = path + combined_data + '.bam'
    merge_cmd = ["samtools", "merge", data_output] + data_input
    #subprocess.call(merge_cmd)
    
    control_list = sorted(cell_chipname_info[cell_name][chipname]["control"])
    combined_control = "-".join(control_list)
    control_input = [path + control + '.bam' for control in control_list]
    control_output = path + combined_control + '.bam'
    merge_cmd = ["samtools", "merge", control_output] + control_input
    #subprocess.call(merge_cmd)

    # make input table
    #print >> f, "%s\t%s\t%s\t%s" % (cell_name, chipname, combined_data + '.bam', combined_control + '.bam')
    print >> f, "%s\t%s\t%s\t%s" % (cell_name, chipname, data_list[0] + '.bam', control_list[0] + '.bam')

f.close()
"""
    

"""
def read_data_info (fname):

    name_info = {}
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

        assert name not in name_info
        name_info[name] = {}
        name_info[name]['cell'] = cell
        name_info[name]['source'] = source
        name_info[name]['data'] = [data.strip() for data in data_list]
        name_info[name]['control'] = [control.strip() for control in control_list]
        name_info[name]['peak'] = [peak.strip() for peak in peak_list]
    return name_info
chipname_info = read_data_info("ENCODE_data_information.csv")

data_list = set([])
for name in chipname_info:
    data_list |= set(chipname_info[name]['peak'])
data_list = list(data_list)

make_download_file ("GM12878_histone_chip_seq_bedfiles.txt", data_list, extension='bed.gz')


data_list = set([])
for name in chipname_info:
    data_list |= set(chipname_info[name]['data'])
    data_list |= set(chipname_info[name]['control'])
data_list = list(data_list)

make_download_file ("H1_histone_chip_seq_files.txt", data_list)

def make_ChromHMM_input (fname, name_info):
    f = open(fname, 'w')
    for name in name_info:
        data_list = [str(data.strip()) for data in sorted(name_info[name]['data'])]
        control_list = [str(control.strip()) for control in sorted(name_info[name]['control'])]
        new_data_fname = ':'.join(data_list)
        new_control_fname = ':'.join(control_list)
        data_cmd = ['samtools merge', new_data_fname+'.bam'] + [data + '.bam' for data in data_list]
        subprocess.call(data_cmd)
        subprocess.call(['rm'] + [data + '.bam' for data in data_list])
        control_cmd = ['samtools merge', new_control_fname+'.bam'] + [control + '.bam' for control in control_list]
        subprocess.call(control_cmd)
        subprocess.call(['rm'] + [control + '.bam' for control in control_list])
        print >> f, "H1\t%s\t%s" % (new_data_fname, new_control_fname)
    f.close()     
"""
