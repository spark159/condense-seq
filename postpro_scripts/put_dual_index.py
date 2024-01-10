import load_file

def rev_comp (seq):
    rev_seq=seq[::-1]; new_seq=''
    for nt in rev_seq:
        nt = nt.upper()
        if nt == 'A':
            new_seq += 'T'
        if nt == 'T':
            new_seq += 'A'
        if nt == 'C':
            new_seq += 'G'
        if nt == 'G':
            new_seq += 'C'
    return new_seq


# read well and index information
well_field_value = load_file.read_tabular_file("SampleSheet_dual_index_info.csv")


# read well list
id_well = load_file.read_tabular_file("input_well.csv",
                                       mode='col',
                                       header=False,
                                       rowID=False)[0]
wells = [id_well[id] for id in sorted(id_well)]


# print out i7 sequence and i5 sequence
f = open("dual_out.txt", 'w')
for well in wells:
    I7_ID = well_field_value[well]['I7_Index_ID']
    I7_seq = well_field_value[well]['index']
    I5_ID = well_field_value[well]['I5_Index_ID']
    I5_seq = well_field_value[well]['index2']
    I5_seq = rev_comp(I5_seq) # rev comp of I5 (for JHU sequencing core)
    print >> f, '\t'.join([well, I7_ID, I7_seq, I5_ID, I5_seq])
f.close()
    
