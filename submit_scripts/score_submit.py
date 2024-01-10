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


### parameters
path = "/home/spark159/scratch/2022_12_30_H1_scores/"
#path = "/home/spark159/scratch/2022_11_27_mouse_scores/"

# sample information
exp_list = [('H1', 'NCP', 'sp'),
            ('H1', 'NCP', 'HP1a'),
            ('H1', 'NCP', 'LKH'),
            ('H1', 'NCP', 'Ki67')]


for cell, sample, agent in exp_list:

    # set species and gender
    if cell in ['H1', 'GM']:
        species = 'human'
    elif cell in ['mCD8T']:
        species = 'mouse'

    if cell in ['H1']:
        gender = 'male'
    elif cell in ['GM', 'mCD8T']:
        gender = 'female'

    # set chromosome list
    if species == 'human':
        chr_list = ['chr' + str(i) for i in range(1, 23)]
    elif species == 'mouse':
        chr_list = ['chr' + str(i) for i in range(1, 20)]
    chr_list += ['chrX']

    if gender == 'male':
        chr_list += ['chrY']

    # gtf file
    if cell in ['H1', 'GM']:
        gtfname = path + "ENCFF159KBI.gtf"
    elif cell in ['mCD8T']:
        gtfname = path + "gencodeM21pri-UCSC-tRNAs-ERCC-phiX.gtf"

    # reference file
    if cell in ['H1', 'GM']:
        refname = path + "hg38"
    elif cell in ['mCD8T']:
        refname = path + "mm10"

    # ENCODE file
    if species == 'human':
        ENCODE_path = '/home/spark159/scratch/ENCODE_data/'
        bs_fnames = [ENCODE_path+"ENCODE_BS_information.csv"]
        chip_fnames = [ENCODE_path+"ENCODE_histone_chip_information.csv",
                       ENCODE_path+"ENCODE_TF_chip_information.csv"]

    elif species == 'mouse':
        ENCODE_path = '/home/spark159/scratch/MouseCD8TcellData/'
        chip_fnames = [ENCODE_path + 'MouseCD8Tcell_histone_chip_seq.csv']


    # make annotation file with score/zscore files
    if False:

        if cell.startswith('GM'):
            cellname = 'GM12878'
        elif cell.startswith('H1'):
            cellname = 'H1-hESC'
        elif cell.startswith('mCD8T'):
            cellname = 'Mouse CD8 T cell (invitro activated)'


        # read BS data
        cellname_bs_info = read_ENCODE_table(bs_fnames)

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


        # read chip seq data
        cellname_chip_info = read_ENCODE_table(chip_fnames)

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

        # call annotation script
        for chr_name in chr_list:

            # score
            fname = path + '_'.join([cell, sample, agent, chr_name]) + '_score.cn'
            outfname = path + '_'.join([cell, sample, agent, chr_name])

            subprocess.call(["sbatch", "annot_submit_edit_edit",
                             "-f", fname,
                             "-x", refname+'.fa',
                             "-n", chr_name,
                             "-o", outfname,
                             "-b", bs_files,
                             "-c", chip_files])

            # zscore
            fname = path + '_'.join([cell, sample, agent, chr_name]) + '_zscore.cn'
            outfname = path + '_'.join([cell, sample, agent, chr_name]) + '_zscore'

            subprocess.call(["sbatch", "annot_submit_edit_edit",
                             "-f", fname,
                             "-x", refname+'.fa',
                             "-n", chr_name,
                             "-o", outfname,
                             "-b", bs_files,
                             "-c", chip_files])



    # make profile using annotation files
    if False:
        for chr_name in chr_list:

            # score annotation file
            fname = path + '_'.join([cell, sample, agent, chr_name]) + '_anot.cn'
            outfname = path + '_'.join([cell, sample, agent, chr_name])

            subprocess.call(["sbatch", "profile_submit_edit",
                             "-f", fname,
                             "-g", gtfname,
                             "-x", refname+'.fa',
                             "-c", chr_name,
                             "-o", outfname])

            # zscore annotation file
            fname = path + '_'.join([cell, sample, agent, chr_name]) + '_zscore_anot.cn'
            outfname = path + '_'.join([cell, sample, agent, chr_name]) + '_zscore'

            subprocess.call(["sbatch", "profile_submit_edit",
                             "-f", fname,
                             "-g", gtfname,
                             "-x", refname+'.fa',
                             "-c", chr_name,
                             "-o", outfname])




    # make profile using score files
    if False:
        for chr_name in chr_list:
            fname = path + '_'.join([cell, sample, agent, chr_name]) + '_score.cn'
            outfname = path + '_'.join([cell, sample, agent, chr_name])

            subprocess.call(["sbatch", "profile_submit_edit",
                             "-f", fname,
                             "-g", gtfname,
                             "-x", refname+'.fa',
                             "-c", chr_name,
                             "-o", outfname])

    # make profile using zscore files
    if False:
        for chr_name in chr_list:
            fname = path + '_'.join([cell, sample, agent, chr_name]) + '_zscore.cn'
            outfname = path + '_'.join([cell, sample, agent, chr_name]) + '_zscore'

            subprocess.call(["sbatch", "profile_submit_edit",
                             "-f", fname,
                             "-g", gtfname,
                             "-x", refname+'.fa',
                             "-c", chr_name,
                             "-o", outfname])


