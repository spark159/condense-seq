from gseapy import Biomart
bm = Biomart()    

# note the dataset and attribute names are different
h2m = bm.query(dataset='hsapiens_gene_ensembl',
               attributes=['ensembl_gene_id',
                           'external_gene_name',
                           'mmusculus_homolog_ensembl_gene',
                           'mmusculus_homolog_associated_gene_name'])

# get a dict symbol mappings
h2m_dict = {}
for i, row in h2m.loc[:,["external_gene_name", "mmusculus_homolog_associated_gene_name"]].iterrows():
    if row.isna().any():
        continue
    h2m_dict[row['external_gene_name']] = row["mmusculus_homolog_associated_gene_name"]

print ("Human to Mouse gene symbol mapping is done")

# parameters
fname_list = ["c2.all.v2022.1.Hs.symbols.gmt",
              "c2.cp.kegg.v2022.1.Hs.symbols.gmt",
              "c7.all.v2022.1.Hs.symbols.gmt",
              "c7.immunesigdb.v2022.1.Hs.symbols.gmt"]

#fname_list = ["c2.all.v2022.1.Hs.symbols.gmt"]


for fname in fname_list:
    print ("Converting %s" % (fname))
    
    outfname = fname.rsplit('.', 1)[0] + '.' + 'mouse_homolog' + '.gmt'

    f = open(outfname, 'w')

    lost_gs_num = 0
    lost_gene_num = 0

    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        gname, source, gene_list = cols[0], cols[1], cols[2:]

        homolog_list = []
        for gene in gene_list:
            try:
                homolog = h2m_dict[gene]
                homolog_list.append(homolog)
            except:
                continue

        lost_gene_num += len(gene_list) - len(homolog_list) 

        if len(homolog_list) <= 0:
            lost_gs_num +=1
            continue

        s = '\t'.join([gname, source] + homolog_list)
        print (s, file=f)

        #s = '\t'.join([gname, source] + homolog_list)
        #print >> f, s

    print ("lost gene-set number: %d" % (lost_gs_num))
    print ("lost gene number: %d" % (lost_gene_num))


    f.close()

