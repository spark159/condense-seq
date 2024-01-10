import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math

genome_size = load_file.read_genome_size("data/hg19.fa")
ID_field_values1, field_ID_values1 = load_file.read_hgtable ("data/hgTables", "chr1", mode="both")
ID_field_values2, field_ID_values2 = load_file.read_GTF ("data/Homo_sapiens.GRCh37.87.gtf", "chr1", mode="both")
