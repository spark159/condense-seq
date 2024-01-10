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
        name, cell, source, data_list, control_list, pipeline, peak_list = cols
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

# parameters
path = '/home/spark159/scratch/MouseEpigeneticData/TFChipSeq/'
table_fname = "MouseCD8Tcell_tf_chip_seq.csv"
#cell_name = "Mouse CD8 T cell (invitro activated)"

#cell_name = "Mouse embryonic stem cell (Smarca4V5_WT)"
#cell_name = "Mouse embryonic stem cell (Arid1aflfl_EtOH)"
cell_name = "Mouse CD4 reg T cell"

# load data information
chipname_info = read_data_info(path + table_fname)[cell_name]

# download all fastq files
fastq_fnames = set([])
for chipname in chipname_info:
    fastq_fnames |= set(chipname_info[chipname]['data'])
    fastq_fnames |= set(chipname_info[chipname]['control'])

fastq_fnames = list(fastq_fnames)
for fastq_fname in fastq_fnames:
    if glob.glob(fastq_fname + '*.fastq.*'):
        continue
    subprocess.call(["sbatch", "sra-submit", '-f', fastq_fname, '-o', path])
    
sys.exit(1)
