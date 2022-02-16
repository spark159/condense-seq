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
chipname_info = read_data_info("ENCODE_Histone_Chipseq_data_info.csv")

def make_download_file (fname, data_list, extension='bam'):
    f = open(fname, 'w')
    for data in data_list:
        print >> f, "https://www.encodeproject.org/files/%s/@@download/%s.%s" % (data, data, extension)
    f.close()
    return

data_list = set([])
for name in chipname_info:
    data_list |= set(chipname_info[name]['peak'])
data_list = list(data_list)

make_download_file ("H1_histone_chip_seq_bedfiles.txt", data_list, extension='bed.gz')

"""
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
