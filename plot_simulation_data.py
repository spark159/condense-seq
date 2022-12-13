import sys
import copy
import numpy as np
import matplotlib.pyplot as plt
import load_file
import statis

### parameters
path = "/media/spark159/sw/advait_simulation_data/2022-11-25_results/H1_hg38_WT_spermine_tp4_chr1_14M-64M_2000beads_25000bp/"

score_fname = "condensability_scores.npy"
exp_fname = "experiment_contact_map.npz"
sim_fname = "simulation_contact_map.npz"
analysis_fname = "analysis_results.npz"

scores = np.load(path+score_fname)
exp_data = np.load(path+exp_fname)
sim_data = np.load(path+sim_fname)

#analysis = np.load(path+analysis_fname, allow_pickle=True)
