import subprocess
import glob

def read_data_info (fname):
    name_info = {}
    First = True
    for line in open(fname):
        cols = line.strip().split('\t')
        if First:
            First = False
            continue
        name, cell, source, data_list, control_list, bed_files = cols
        data_list = data_list.split(',')
        control_list = control_list.split(',')
        bed_files = bed_files.split(',')
        assert name not in name_info
        name_info[name] = {}
        name_info[name]['cell'] = cell
        name_info[name]['source'] = source
        name_info[name]['data'] = [data.strip() for data in data_list]
        name_info[name]['control'] = [control.strip() for control in control_list]
        name_info[name]['bed'] = [bed.strip() for bed in bed_files]
    return name_info

# read BS data
bs_path = '/home-4/spark159@jhu.edu/work/ENCODE_data/'
bs_info = read_data_info(bs_path+"ENCODE_BS_data_info.csv")

bs_fname = {}
for bs in bs_info:
    bed_files = bs_info[bs]['bed']
    assert len(bed_files) == 1
    bs_fname[bs] = bed_files[0] + '.bed'

bs_files = []
for bs, fname in bs_fname.items():
    bs_files.append(bs)
    bs_files.append(bs_path+fname)
bs_files = ','.join(bs_files)    


# read chip seq data
chip_path = '/home-4/spark159@jhu.edu/work/ENCODE_data/'
chip_info = read_data_info(chip_path+"ENCODE_Histone_Chipseq_data_info.csv")

chip_fname = {}
for chip in chip_info:
    bed_files = chip_info[chip]['bed']
    assert len(bed_files) == 1
    chip_fname[chip] = bed_files[0] + '.bed'

chip_files = []
for chip, fname in chip_fname.items():
    chip_files.append(chip)
    chip_files.append(chip_path+fname)
chip_files = ','.join(chip_files)    


# make annotation file
path = '/home-4/spark159@jhu.edu/work/2021_06_07_H1_sp_detail/'
refname = path + 'hg38.fa'
cell = 'H1'
sample = 'NCP'
agent = 'sp'
chr = 'chr1'

infname = path + '_'.join([cell, sample, agent, chr, 'Ncov']) + '.cn'
outfname = path + '_'.join([cell, sample, agent, chr, 'extended'])


subprocess.call(["sbatch", "annot_submit", "-f",infname, "-x", refname, "-c", chr, "-o", outfname, "-b", bs_files, "-h", chip_files])
