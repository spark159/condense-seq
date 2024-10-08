import numpy as np
import matplotlib.pyplot as plt

# read titration file
def read_titration (fname):
    conc_list = []
    mean_list = []
    std_list = []
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        if not line:
            continue
        cols = line.strip().split()
        conc, mean, std = float(cols[0]), float(cols[-3]), float(cols[-2])
        conc_list.append(conc)
        mean_list.append(mean)
        std_list.append(std)
    return conc_list, mean_list, std_list

### parameters
agent_unit = {'sp':'mM',
              'spd':'mM',
              'CoH':'mM',
              'PEG':'%',
              'Mg':'mM',
              'Ca':'mM',
              'HP1a':'$\\mu$M',
              'HP1bSUV':'$\\mu$M',
              'LKH':'$\\mu$M',
              'Ki67':'$\\mu$M'}

agent_logbase = {'sp':10,
                 'spd':10,
                 'CoH':10,
                 'PEG':None,
                 'Mg':None,
                 'Ca':None,
                 'HP1a':2,
                 'HP1bSUV':2,
                 'LKH':2,
                 'Ki67':2}

agent_fullname = {'sp':'Spermine(4+)',
                  'spd':'Spermidine(3+)',
                  'CoH':'Cobalt Hexammine(3+)',
                  'PEG':'PEG 8000',
                  'Mg':'Magnesium',
                  'Ca':'Calcium',
                  'HP1a':'HP1$\\alpha$',
                  'HP1bSUV':'HP1$\\beta$/SUV39H1',
                  'LKH':'Linker histone1',
                  'Ki67':'Ki67'}

#path = ""
#path = "/home/spark159/../../media/spark159/sw/"
path = "/data/"

### experiment information
# mouse data
agent = "sp"
exp_list = [('mCD8T', 'WT-NCP'),
            ('mCD8T', 'inht-NCP'),
            ('mCD8T', 'KO-NCP')]
labels = ['WT', '+inht', 'ODC KO']
colors = ['tab:blue', 'tab:orange', 'tab:green']

title = "Mouse CD8 T cell Nucleosome condensation"

# proteins
agent = 'HP1bSUV'
exp_list = [('H1', 'NCP'),
            ('H1', 'DNA')]
labels = ['NCP', 'DNA']
colors = ['red', 'blue']
title = "H1-hESC %s" % (agent_fullname[agent])





### plot titration
#fig = plt.figure()
#fig = plt.figure(figsize=(4, 3))

fig = plt.figure(figsize=(2.4,1.8))
#fig = plt.figure(figsize=(2, 1.4))


for i in range(len(exp_list)):
    cell, sample  = exp_list[i]
    fname = path + '_'.join([cell, sample, agent, 'titration']) + '.csv'
    conc_list, mean_list, std_list = read_titration (fname)

    if agent_logbase[agent] != None:
        conc_list = conc_list[1:]
        mean_list = mean_list[1:]
        std_list = std_list[1:]

    label, color = labels[i], colors[i]

    if color.startswith('tab'):
        ecolor = color
    else:
        ecolor = 'tab:' + color

    #plt.plot(conc_list, mean_list, 'o-', lw=2, color=color, label=label)
    #plt.errorbar(conc_list, mean_list, yerr=std_list, fmt='.',
    #             mfc=color, mec=color, color=ecolor, alpha=0.8)

    plt.plot(conc_list, mean_list, 'o-', color=color, lw=1.5, markersize=3.5, label=label)
    plt.errorbar(conc_list, mean_list, yerr=std_list, fmt='.',
                 markersize=3.5, mfc=color, mec=color, lw=1.5, color='tab:'+color, alpha=0.8)

if agent_logbase[agent] != None:
    plt.xscale("log", basex=agent_logbase[agent])

plt.gca().tick_params(axis='both', which='major', labelsize=5)
plt.gca().tick_params(axis='both', which='minor', labelsize=5)
plt.xlabel("%s concentration (%s)" % (agent_fullname[agent], agent_unit[agent]), fontsize=8)
plt.ylabel("Soluble fraction", fontsize=8)
plt.title(title, fontsize=10)
plt.legend(fontsize=8)
plt.savefig('titration.svg', format='svg', bbox_inches='tight')
plt.close()
