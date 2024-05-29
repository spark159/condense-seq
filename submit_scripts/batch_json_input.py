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

# parameters
path = '/home/spark159/scratch/MouseEpigeneticData/TFChipSeq/'
table_fname = "MouseCD8Tcell_tf_chip_seq.csv"

cell_names = ["Mouse CD8 T cell (invitro activated)",
              "Mouse embryonic stem cell (Smarca4V5_WT)",
              "Mouse embryonic stem cell (Arid1aflfl_EtOH)",
              "Mouse CD4 reg T cell"]

#cell_name = "Mouse CD8 T cell (invitro activated)"
#cell_name = "Mouse embryonic stem cell (Smarca4V5_WT)"
#cell_name = "Mouse embryonic stem cell (Arid1aflfl_EtOH)"

#cell_names = ["Mouse CD4 reg T cell"]

chip_type = "tf"
ref_name = "mm10"
extension = '.fastq.gz'

optional = True

# load data information
cell_chipname_info = read_data_info(path + table_fname)

for cell_name in cell_names:
    chipname_info = cell_chipname_info[cell_name]

    # write input json file for chip-seq analysis
    for chipname in chipname_info:
        data_list = chipname_info[chipname]['data']
        control_list = chipname_info[chipname]['control']

        fname = "%s_input_json.json" % (chipname)
        f = open(path+fname, 'w')
        print >> f, '{'

        if optional:
            print >> f,  '\t"chip.call_peak_cpu": 10.0,'
            print >> f,  '\t"chip.call_peak_spp_mem_factor": 50.0,'
            print >> f,  '\t"chip.call_peak_macs2_mem_factor": 50.0,'
            print >> f,  '\t"chip.bam2ta_time_hr": 10,'
            print >> f,  '\t"chip.jsd_time_hr": 10,'
            print >> f,  '\t"chip.xcor_time_hr": 10,'

        print >> f, '\t"chip.pipeline_type" : "%s",' % (chip_type)
        print >> f, '\t"chip.genome_tsv" : "https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v4/%s.tsv",' % (ref_name)

        for k in range(len(data_list)):
            print >> f, '\t"chip.fastqs_rep%d_R1" : ["%s"],' % (k+1, path+data_list[k]+extension)

        for k in range(len(control_list)):
            print >> f, '\t"chip.ctl_fastqs_rep%d_R1" : ["%s"],' % (k+1, path+control_list[k]+extension)

        print >> f, '\t"chip.paired_end" : false,'
        print >> f, '\t"chip.title" : "%s %s",' % (cell_name, chipname)
        print >> f, '\t"chip.description" : "%s %s"' % (cell_name, chipname)

        print >> f, '}'
        f.close()
