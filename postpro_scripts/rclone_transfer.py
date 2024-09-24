import glob
import subprocess
from datetime import date

# rclone_copy
def rclone_copy (from_path,
                 to_path,
                 dir_name):

    today = date.today().strftime("%b-%d-%Y")
    if len(glob.glob('./' + today)) == 0:
        subprocess.call(['mkdir', today])

    cmd = ['rclone',
           'copy',
           '-vP',
           from_path + dir_name,
           to_path + dir_name]

    log_name = '-'.join(from_path.strip('/').split('/') + [dir_name]) 

    subprocess.call(cmd,
                    stdout=open("./" + today + "/" + log_name + "_out.txt", 'w'),
                    stderr=open("./" + today + "/" + log_name + "_err.txt", 'w'))
    
    return

# get all subdirectories in the path
def get_subdirs (path):

    cmd = ['rclone',
           '-lsf',
           '--dirs-only',
           path]

    rclone_ls_dir = subprocess.Popen(cmd,
                                     stdout=subprocess.PIPE,
                                     stderr=open("/dev/null", 'w'))

    subdirs = []
    for line in rclone_ls_dir.stdout:
        line = line.strip()
        subdirs.append(line[:-1])

    return subdirs

# copy the whole directory by rclone_copy recursively
def rclone_dir (from_path,
                to_path,
                dir_name):

    rclone_copy (from_path,
                 to_path,
                 dir_name)
    
    subdir_names = get_subdirs (from_path + dir_name)

    if len(subdir_names) == 0:
        return

    for subdir_name in subdir_names:
        rclone_dir (from_path + dir_name + '/',
                    to_path + dir_name + '/',
                    subdir_name)
        

### parameters
#from_path = 'rockfish_jhu:/home/spark159/data/'
#to_path = 'onedrive_jhu:/Ha-SPark/Condense-seq_project/data/fastq_files/'

#dir_names = ['spark205_172495',
#             'spark205_172936',
#             'spark205_180247',
#             'spark205_181643',
#             'spark205_183500',
#             'spark205_185389',
#             'spark205_185767',
#             'spark205_186596',
#             'spark205_186600',
#             'spark205_186759',
#             'spark205_187243',
#             'spark205_188157',
#             'spark205_188441',
#             'spark205_188633',
#             'spark205_188810',
#             'spark205_H1spqc']

#dir_names = ['pancancer_data']
#dir_names = ['spark205_172495']
#dir_names = ['rclone_test']
#dir_names = ['old_analysis']

#from_path = 'rockfish_jhu:/home/spark159/'
#to_path = 'onedrive_jhu:/Ha-SPark/Condense-seq_project/scripts/SLRUM_submit/'
#dir_names = ['submit_scripts']

#from_path = 'onedrive_jhu:/'
#to_path = 'sangwoo_dropbox:/'
#dir_names = ['Ha-SPark']

#from_path = 'azenta_halab:'
#to_path = 'my_dropbox:My\ PhD\ works/Research\ projects/Condense-seq_project/data/'
#to_path = 'my_dropbox:My PhD works/Research projects/Condense-seq_project/data/'
#to_path = 'jhu_rockfish:data/ODC_mouse_RNAseq/'
#dir_names = ['30-916504364']

#for dir_name in dir_names:
#    rclone_copy (from_path, to_path, dir_name)

#### data transfer for jonghan
from_path = 'jhu_rockfish:data/2024_01_05_GEO/processed_files/'
to_path = 'my_dropbox:Condense-seq data share/'

cells = ['H1', 'GM']
reps = [1, 2]

for rep in reps:
    for cell in cells:
        fname = '%s_NCP_sp_%drep_deep_1kb_score_table.gtab.gz' % (cell, rep)

        cmd = ['rclone',
               'copy',
               '-v',
               from_path + fname,
               to_path]

        subprocess.call(cmd,
                        stdout=open("./logs/%s-%d_out.txt"
                                    % (cell, rep), 'w'),
                        stderr=open("./logs/%s-%d_err.txt"
                                    % (cell, rep), 'w'))

"""
cell_chrnames = {'H1':['chr%s' % (i) for i in range(1, 23)] + ['chrX', 'chrY'],
                 'GM':['chr%s' % (i) for i in range(1, 23)] + ['chrX'],
                 'mCD8T:WT':['chr%s' % (i) for i in range(1, 20)] + ['chrX'],
                 'mCD8T:DFMO':['chr%s' % (i) for i in range(1, 20)] + ['chrX'],
                 'mCD8T:ODCKO':['chr%s' % (i) for i in range(1, 20)] + ['chrX']}

cells = ['H1', 'GM']
reps = [1, 2]

for ext in ['cov', 'peak', 'Ncov']:
    for rep in reps:
        for cell in cells:
            chrs = cell_chrnames[cell]
            for chr in chrs:
                for tnum in [0, 4, 8]:
                    #fname = '%s_NCP_sp_%drep_deep_%s_score_table.gtab.gz' % (cell, rep, chr)
                    fname = '%s_NCP_sp_%d_%drep_deep_%s_%s.gtab.gz' % (cell, tnum, rep, chr, ext)
                    cmd = ['rclone',
                           'copy',
                           '-v',
                           from_path + fname,
                           to_path]

                    subprocess.call(cmd,
                                    stdout=open("./logs/%d-%s-%s-%d-%s_out.txt"
                                                % (rep, cell, chr, tnum, ext), 'w'),
                                    stderr=open("./logs/%d-%s-%s-%d-%s_err.txt"
                                                % (rep, cell, chr, tnum, ext), 'w'))
"""
