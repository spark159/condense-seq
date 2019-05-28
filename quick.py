import load_file
import graphics
import statis
import sys
import copy
import Interval_dict
import matplotlib.pyplot as plt
import numpy as np
import math
from sklearn import linear_model

gID_field_values, field_gID_values = load_file.read_GTF ("data/Homo_sapiens.GRCh37.87.gtf", "chr1", mode="both")

gID_length = {}
for gID in gID_field_values:
    field_values = gID_field_values[gID]
    length = abs(field_values['TTS'] - field_values['TSS'])
    gID_length[gID] = length


    
