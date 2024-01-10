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
#exp_list = [('H1', 'NCP', 'sp'),
#            ('H1', 'NCP', 'spd'),
#            ('H1', 'NCP', 'CoH'),
#            ('H1', 'NCP', 'PEG'),
#            ('H1', 'NCP', 'Mg'),
#            ('H1', 'NCP', 'Ca'),
#            ('H1', 'NCP', 'HP1a'),
#            ('H1', 'NCP', 'HP1bSUV'),
#            ('H1', 'NCP', 'LKH'),
#            ('H1', 'NCP', 'Ki67'),
#            ('H1', 'NCP', 'FUS')]

exp_list = [('H1', 'NCP', 'sp'),
            ('H1', 'NCP', 'HP1a'),
            ('H1', 'NCP', 'LKH'),
            ('H1', 'NCP', 'Ki67')]



#bin_size
#bin_size = 10000
bin_size = 5000
#bin_size = 25000
#bin_size = 1000

note = ""


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
        #bs_fnames = []
        chip_fnames = [ENCODE_path+"ENCODE_histone_chip_information.csv",
                       ENCODE_path+"ENCODE_TF_chip_information.csv"]

    elif species == 'mouse':
        ENCODE_path = '/home/spark159/scratch/MouseCD8TcellData/'
        bs_fnames = []
        chip_fnames = [ENCODE_path + 'MouseCD8Tcell_histone_chip_seq.csv']

    # other data
    if species == 'human':
        bg_fname = {'SON':'4DNFI625PP2A.bedgraph',
                    'LaminB1':'4DNFIXNBG8L1.bedgraph',
                    'Nucleolar':'K562_NucleolarDamID.bedgraph',
                    'eigen':'eigen_H1_100kb.bedgraph'}


    # make annotation file with score/zscore files
    if True:

        if cell.startswith('GM'):
            cellname = 'GM12878'
        elif cell.startswith('H1'):
            cellname = 'H1-hESC'
        elif cell.startswith('mCD8T'):
            cellname = 'Mouse CD8 T cell (invitro activated)'

        # read BS data
        bs_files = []
        if bs_fnames:
            cellname_bs_info = read_ENCODE_table(bs_fnames)
            bs_info = cellname_bs_info[cellname]

            for bs in bs_info:
                bed_files = bs_info[bs]['peak']
                assert len(bed_files) == 1
                bs_files.append(bs)
                bs_files.append(ENCODE_path + bed_files[0]+'.bed')
            bs_files = ','.join(bs_files)


        # read chip seq data
        chip_files = []
        if chip_fnames:
            cellname_chip_info = read_ENCODE_table(chip_fnames)
            chip_info = cellname_chip_info[cellname]

            for chip in chip_info:
                bed_files = chip_info[chip]['peak']
                assert len(bed_files) == 1
                chip_files.append(chip)
                chip_files.append(ENCODE_path + bed_files[0]+'.bed')
            chip_files = ','.join(chip_files)

        
        # read additional bedgraph files
        bg_files = []
        for bg, fname in bg_fname.items():
            bg_files.append(bg)
            bg_files.append(ENCODE_path+fname)
        bg_files = ','.join(bg_files)    


        # call annotation script
        cmd = []

        if bs_files:
            cmd.append("-b")
            cmd.append(bs_files)

        if chip_files:
            cmd.append("-c")
            cmd.append(chip_files)
            
        if bg_files:
            cmd.append("-g")
            cmd.append(bg_files)

        cmd.append('-w')
        cmd.append(str(bin_size))
        cmd.append('-s')
        cmd.append(str(bin_size))
        

        # score
        fname = path + '_'.join([cell, sample, agent, str(bin_size/1000) + 'kb']) + '_score.cn'
        outfname = [cell, sample, agent, str(bin_size/1000) + 'kb']
        if note:
            outfname.append(note)
        outfname = path + '_'.join(outfname)

        subprocess.call(["sbatch",
                         "annot_submit_edit_edit_bin",
                         "-f", fname,
                         "-x", refname+'.fa',
                         "-o", outfname] + cmd)
        
        # zscore
        fname = path + '_'.join([cell, sample, agent, str(bin_size/1000) + 'kb']) + '_zscore.cn'
        outfname = [cell, sample, agent, str(bin_size/1000) + 'kb']
        if note:
            outfname.append(note)
        outfname = path + '_'.join(outfname) + '_zscore'


        subprocess.call(["sbatch",
                         "annot_submit_edit_edit_bin",
                         "-f", fname,
                         "-x", refname+'.fa',
                         "-o", outfname] + cmd)



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


