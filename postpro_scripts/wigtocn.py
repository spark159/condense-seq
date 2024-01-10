
def convert_file (fname, chr_target='chr1'):
    f = open(fname.split('.')[0] + '.cn', 'w')
    s ="SNP"+"\t"+"Chromosome"+"\t"+"PhysicalPosition"+"\t"+"Weight"
    print >> f, s
    ID = 0
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        if cols[0] == "fixedStep":
            chr = cols[1][6:]
            pos = 0
            continue
        if chr != chr_target:
            continue
        count = int(cols[0])
        if count <= 0:
            pos += 1
            continue
        s = str(ID) + "\t" + chr + "\t" + str(pos) + "\t" + str(count)
        print >> f, s
        pos +=1
        ID +=1
    f.close()
    return None

convert_file ("GSM907784_mnase_mids_NA18508_126_to_184.wig")


        
            
            
        
        
