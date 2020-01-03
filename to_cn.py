
def to_cn (fname, new_name):
    f = open(new_name, 'w')
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            s = 'SNP\tChromosome\tPhysicalPosition\t' + "\t".join(cols[4:])
            print >> f, s
            First = False
            continue
        s = "\t".join(cols[:2])
        s += "\t" + str(int(round(0.5*(float(cols[2]) + float(cols[3]))))) + "\t"
        s += "\t".join(cols[4:])
        print >> f, s
    f.close()

to_cn("/home/spark159/Downloads/DNA_Spermine(4+)_1kb_bin.cn", "/home/spark159/Downloads/DNA_Spermine(4+)_1kb.cn")
to_cn("/home/spark159/Downloads/NCP_Spermine(4+)_1kb_bin.cn", "/home/spark159/Downloads/NCP_Spermine(4+)_1kb.cn")
to_cn("/home/spark159/Downloads/DNA_Spermidine(3+)_1kb_bin.cn", "/home/spark159/Downloads/DNA_Spermidine(3+)_1kb.cn")
to_cn("/home/spark159/Downloads/NCP_Spermidine(3+)_1kb_bin.cn", "/home/spark159/Downloads/NCP_Spermidine(3+)_1kb.cn")

        
