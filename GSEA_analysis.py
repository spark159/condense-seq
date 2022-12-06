import pandas as pd
import gseapy as gp
from gseapy import gseaplot
import matplotlib.pyplot as plt


# read gene set libraries
#mlib_names = gp.get_library_name(organism='mouse')
#gmt_dict = gp.parser.gsea_gmt_parser('GO_Biological_Process_2018', organism='mouse')

# parameters
rnk_fname = 'KO-WT.rnk'
gmt_fnames = ['m5.go.bp.v2022.1.Mm.symbols.gmt']

# prerank run
pre_res = gp.prerank(rnk=rnk_fname, # or rnk = rnk,
                     gene_sets=gmt_fnames,
                     threads=4,
                     min_size=5,
                     max_size=1000,
                     permutation_num=1000, # reduce number to speed up testing
                     outdir='.', # don't write to disk
                     seed=6,
                     verbose=True, # see what's going on behind the scenes
                    )


