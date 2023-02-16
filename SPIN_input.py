import numpy as np
import matplotlib.pyplot as plt
import load_file

def read_hic_matrix (fname, binsize):
    pair_value = {}

    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        start1, start2, value = cols
        start1, start2 = int(start1), int(start2)
        value = float(value)

        assert start1 % binsize == 0
        assert start2 % binsize == 0

        idx1 = start1/binsize
        idx2 = start2/binsize

        pair = (idx1, idx2)
        pair_value[pair] = value

    return pair_value

### parameters
binsize = 25000
chr_choice = 'chr1'
hicm_fname = 'H1_hicm_%s_%dkb.txt' % (chr_choice, binsize/1000)
anot_fname = "H1_NCP_sp_%dkb_anot.txt" % (binsize/1000)
domain_fields = {'NSpeckle':['SON'],
                 'Trnx':['POLR2A', 'POLR2AphosphoS5', 'H3K9ac', 'H3K4me3', 'H3K27ac'],
                 'Polycomb':['CBX8', 'EZH2', 'RNF2', 'SUZ12'],
                 'Hetero':['H3K9me3', 'CBX5'],
                 'Nucleolus':['Nucleolar'],
                 'Lamin':['LaminB1'],
                 'Compartment':['eigen'],
                 'other':['ATcontent']}

domains = ['NSpeckle', 'Trnx', 'Polycomb', 'Hetero', 'Nucleolus', 'Lamin']
#domains = ['NSpeckle', 'Nucleolus', 'Lamin']

fields = []
for domain in domains:
    fields += domain_fields[domain]

    
### load file and reading data
print "reading data"

pair_value = read_hic_matrix (hicm_fname, binsize)
field_ID_value = load_file.read_tabular_file(anot_fname, mode='col')
ID_chr = field_ID_value['Chromosome']
ID_start = field_ID_value['Start']

idx_sigs = {}
for ID in sorted(ID_chr.keys()):
    chr, start = ID_chr[ID], ID_start[ID]
    if chr != chr_choice:
        continue
    sigs = []
    for field in fields:
        sigs.append(field_ID_value[field][ID])
    idx = start/binsize
    idx_sigs[idx] = sigs
del ID_chr
del ID_start
del field_ID_value

hic_idxes = []
for idx1, idx2 in pair_value:
    hic_idxes.append(idx1)
    hic_idxes.append(idx2)

min_idx = min(hic_idxes + idx_sigs.keys())
max_idx = max(hic_idxes + idx_sigs.keys())


### standardization of sigs
sigs_list = np.asarray(idx_sigs.values())
mean_sigs = np.mean(sigs_list, axis=0)
std_sigs = np.std(sigs_list, axis=0)

for idx in idx_sigs:
    for i in range(len(idx_sigs[idx])):
        idx_sigs[idx][i] = float(idx_sigs[idx][i] - mean_sigs[i])/std_sigs[i]


### write input files
print "writing data"
headname = "SPIN_%s_%dkb" % (chr_choice, binsize/1000)

f = open(headname + "_sig_input.txt", 'w')
for idx in range(min_idx, max_idx+1):
    try:
        sigs = idx_sigs[idx]
    except:
        sigs = [0.0]*len(fields)
    sigs = [str(sig) for sig in sigs]
    print >> f, '\t'.join(sigs)
f.close()

f = open(headname + "_hic_input.txt", 'w')
for pair in sorted(pair_value):
    idx1, idx2 = pair
    newidx1, newidx2 = idx1-min_idx, idx2-min_idx
    value = pair_value[pair]
    print >> f, '%d\t%d\t%f' % (newidx1, newidx2, value)
f.close()

f = open(headname + "_bin_input.txt", 'w')
for idx in range(min_idx, max_idx+1):
    start = idx*binsize
    end = start + binsize
    newidx = idx - min_idx
    print >> f, '%s\t%d\t%d\t%d' % (chr_choice, start, end, newidx)
f.close()








    




"""
idx_list = []
for idx1, idx2 in pair_value:
    idx_list.append(idx1)
    idx_list.append(idx2)

total_bins = max(idx_list) + 1
img = np.zeros((total_bins, total_bins))

for pair in pair_value:
    idx1, idx2 = pair
    value = pair_value[pair]

    img[idx1][idx2] = value
    img[idx2][idx1] = value

fig = plt.figure()
plt.imshow(np.log(img), cmap='jet')
plt.imshow(img)
plt.show()
plt.close()
"""
