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
path = ""
#path = "/home/spark159/../../media/spark159/sw/"

agent_unit = {'sp':'mM',
              'spd':'mM',
              'CoH':'mM',
              'PEG':'%',
              'HP1a':'$\\mu$M',
              'Mg':'mM',
              'Ca':'mM'}

agent_logbase = {'sp':10,
                 'spd':10,
                 'CoH':10,
                 'PEG':None,
                 'HP1a':2,
                 'Mg':None,
                 'Ca':None}

agent_fullname = {'sp':'Spermine(4+)',
                  'spd':'Spermidine(3+)',
                  'CoH':'Cobalt Hexammine(3+)',
                  'PEG':'PEG 8000',
                  'HP1a':'HP1 $\\alpha$',
                  'Mg':'Magnesium',
                  'Ca':'Calcium'}


### experiment information
agent = "sp"
cell1, cell2, cell3 = "mCD8T", "mCD8T", "mCD8T"
sample1, sample2, sample3 = "WT-NCP", "inht-NCP", "KO-NCP"

exp_list = [(cell1, sample1),
            (cell2, sample2),
            (cell3, sample3)]

labels = ['WT', '+inht', 'ODC KO']
colors = ['tab:blue', 'tab:orange', 'tab:green']

title = "Mouse CD8 T cell Nucleosome condensation"

### plot titration
fig = plt.figure()
#fig = plt.figure(figsize=(3, 2))
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

    plt.plot(conc_list, mean_list, 'o-', color=color, label=label)
    plt.errorbar(conc_list, mean_list, yerr=std_list, fmt='.',
                 mfc=color, mec=color, color=color, alpha=0.8)

    #plt.plot(conc_list, mean_list, 'o-', color=color, lw=1.5, markersize=3.5, label=label)
    #plt.errorbar(conc_list, mean_list, yerr=std_list, fmt='.',
    #             markersize=3, mfc=color, mec=color, lw=1, color=color, alpha=0.8)

if agent_logbase[agent] != None:
    plt.xscale("log", basex=agent_logbase[agent])

plt.xlabel("%s concentration (%s)" % (agent_fullname[agent], agent_unit[agent]))
plt.ylabel("Soluble fraction")
plt.title(title)
plt.legend(fontsize=10)
plt.savefig('titration.svg', format='svg', bbox_inches='tight')
plt.close()
