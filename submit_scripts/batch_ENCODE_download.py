import sys
import csv, codecs
import subprocess
import glob

def read_csv_file (fname,
                   mode='row',
                   delim=',',
                   header=True,
                   rowID=True,
                   jump=None):
    if rowID:
        col_st = 1
    else:
        col_st = 0
        
    ID_field_value = {}
    First = True
    counter = -1

    #for cols in csv.reader(open(fname), delimiter=delim):
    for cols in csv.reader(codecs.EncodedFile(open(fname), 'utf-8', 'utf-8-sig'),
                           delimiter=delim):
        if First and header:
            field_names = cols[col_st:]
            First = False
            continue
        elif First and not header:
            field_names = range(len(cols[col_st:]))
            First = False
            pass

        counter += 1
        if jump and counter % jump != 0:
            continue

        if rowID:
            ID = cols[0]
        else:
            ID = counter

        if ID not in ID_field_value:
            ID_field_value[ID] = {}

        cols = cols[col_st:]
        #print cols
        for i in range(len(cols)):
            field = field_names[i]
            try:
                value = float(cols[i])
            except:
                value = cols[i]
            if field not in ID_field_value[ID]:
                ID_field_value[ID][field] = value
            else:
                if type(ID_field_value[ID][field]) != list:
                    ID_field_value[ID][field] = [ID_field_value[ID][field]]
                ID_field_value[ID][field].append(value)

    if mode == 'row':
        return ID_field_value

    if mode == 'col' or mode == 'both':
        field_ID_value = {}
        for ID in ID_field_value:
            field_value = ID_field_value[ID]
            for field in field_value:
                value = field_value[field]
                if field not in field_ID_value:
                    field_ID_value[field] = {}
                field_ID_value[field][ID] = value

    if mode == 'col':
        return field_ID_value

    if mode == 'both':
        return ID_field_value, field_ID_value

def make_ENCODE_download_file (fname,
                               data_list,
                               extension='bam'):
    f = open(fname, 'w')
    for data in data_list:
        print >> f, "https://www.encodeproject.org/files/%s/@@download/%s.%s" % (data, data, extension)
    f.close()
    return

# data information files
path = '/home/spark159/data/HumanEpigeneticData/HistoneChipseq/'
#info_fnames = ['human_HistoneChip_info.csv',
#               'human_BSseq_info.csv']
info_fnames = ['human_HistoneChipseq_info.csv']



extension = 'bed.gz'

# read data information and download files from ENCODE
for fname in info_fnames:
    ID_field_value = read_csv_file(path + fname, rowID=False)
    data_list = []
    for ID in ID_field_value:
        bed_files = ID_field_value[ID]['Bed files'].split(',')
        data_list += bed_files

    # download ENCODE data
    #ENCODE_down_fname = path + fname.rsplit('.', 1)[0] + '_ENCODE_down.txt'
    #make_ENCODE_download_file (ENCODE_down_fname,
    #                           data_list,
    #                           extension=extension)
    #subprocess.call("xargs -n 1 curl -O -L < " + ENCODE_down_fname, shell=True)
    #subprocess.call("(cd %s && xargs -n 1 curl -O -L < %s)" % (path, ENCODE_down_fname), shell=True)
