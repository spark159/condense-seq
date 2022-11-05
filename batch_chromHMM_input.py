import subprocess
import glob

def make_download_file (fname, data_list, extension='bam'):
    f = open(fname, 'w')
    for data in data_list:
        print >> f, "https://www.encodeproject.org/files/%s/@@download/%s.%s" % (data, data, extension)
    f.close()
    return

def read_data_info (fname):
    name_info = {}
    First = True
    for line in open(fname):
        cols = line.strip().split('\t')
        if First:
            First = False
            continue
        name, cell, source, data_list, control_list = cols
        data_list = data_list.split(',')
        control_list = control_list.split(',')
        assert name not in name_info
        name_info[name] = {}
        name_info[name]['cell'] = cell
        name_info[name]['source'] = source
        name_info[name]['data'] = [data.strip() for data in data_list]
        name_info[name]['control'] = [control.strip() for control in control_list]
    return name_info

# read chip seq data list
name_info = read_data_info("All_Histone_Chipseq_data_info.csv")

# combine controls
for name in name_info:
    print name
    control_list = sorted(name_info[name]['control'])
    combined_control_name = '_'.join(control_list)
    if not glob.glob(combined_control_name + '.bam'):
        cmd = ['samtools', 'merge', combined_control_name + '.bam'] + [control + '.bam' for control in control_list]
        subprocess.call(cmd)
    name_info[name]['control'] = combined_control_name

# write chromHMM input file
f = open("H1_ENCODE_Histone_input.txt", 'w')
for name in name_info:
    control = name_info[name]['control']
    for data in name_info[name]['data']:
        print >> f, 'H1\t%s\t%s.bam\t%s.bam' % (name, data, control)
f.close()
