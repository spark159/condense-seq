
def to_cn (fname, new_name, GC=False):
    f = open(new_name, 'w')
    include_GC = False
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            if cols[-1].startswith('GC') and GC:
                s = 'SNP\tChromosome\tPhysicalPosition\t' + "\t".join(cols[4:])
                include_GC = True
            else:
                s = 'SNP\tChromosome\tPhysicalPosition\t' + "\t".join(cols[4:-1])
            print >> f, s
            First = False
            continue
        s = "\t".join(cols[:2])
        s += "\t" + str(int(round(0.5*(float(cols[2]) + float(cols[3]))))) + "\t"
        if include_GC:
            s += "\t".join(cols[4:])
        else:
            s += "\t".join(cols[4:-1])
        print >> f, s
    f.close()

#to_cn("/home/spark159/Downloads/DNA_Spermine(4+)_1kb_bin.cn", "/home/spark159/Downloads/DNA_Spermine(4+)_1kb.cn")
#to_cn("/home/spark159/Downloads/NCP_Spermine(4+)_1kb_bin.cn", "/home/spark159/Downloads/NCP_Spermine(4+)_1kb.cn")
#to_cn("/home/spark159/Downloads/DNA_Spermidine(3+)_1kb_bin.cn", "/home/spark159/Downloads/DNA_Spermidine(3+)_1kb.cn")
#to_cn("/home/spark159/Downloads/NCP_Spermidine(3+)_1kb_bin.cn", "/home/spark159/Downloads/NCP_Spermidine(3+)_1kb.cn")

to_cn("H1_DNA_sp_10kb_bin.cn", "H1_DNA_sp_10kb.cn")
to_cn("H1_NCP-new_sp_10kb_bin.cn", "H1_NCP-new_sp_10kb.cn")
to_cn("H1_DNA_spd_10kb_bin.cn", "H1_DNA_spd_10kb.cn")
to_cn("H1_NCP-new_spd_10kb_bin.cn", "H1_NCP-new_spd_10kb.cn")


        
