import glob
import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import copy
import re
import numpy as np
import random
from scipy.optimize import curve_fit

# sorting agents
def cmp_agent (agent1, agent2):
    agent_idx = {'sp':0,
                 'spd':1,
                 'CoH':2,
                 'PEG':3,
                 'HP1a':4}

    idx1 = agent_idx[agent1]
    idx2 = agent_idx[agent2]

    if idx1 < idx2:
        return -1
    elif idx2 > idx1:
        return 1
    else:
        return 0

# get GC content of sequence
def GC_content(seq):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    return (num/float(len(seq)))*100

# parsing the histone modifications
def mhistone_parser (mhistone):
    hname, mutations = re.split('(H2A(?:[.]\w)?|H2B|H3(?:[.]\d)?|H4)', mhistone)[1:]
    #all_mutations = set([])
    pos_mutation = {}
    pattern = '([A-Z])(\d+(?:,\d+)*)(ac|me2[as]|me[1-3]|me|ub|ph|cr|GlcNAc|[A-Z])'  
    for find in re.findall(pattern, mutations):
        aa, pos_list, mutation = find
        assert aa in aa_info.keys()
        if mutation in aa_info.keys():
            mtype = 'mut'
        elif mutation.startswith('me'):
            mtype = 'me'
        else:
            mtype = mutation
        for pos in pos_list.split(','):
            pos = int(pos)
            assert pos not in pos_mutation
            pos_mutation[pos] = (mtype, aa, mutation)
        #all_mutations.add((hname, pos, mtype, aa, mutation))
    #return hname, all_mutations, pos_mutation
    return hname, pos_mutation

# read NGS library index information
def read_index (fname):
    index_name, name_index = {}, {}
    name_titr = {}
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        cols = line.strip().split('\t')
        index = cols[-1]
        name = cols[0]
        titr = int(cols[3])
        assert index not in index_name
        assert name not in name_index
        index_name[index] = name
        name_index[name] = index
        name_titr[name] = titr
    return index_name, name_index, name_titr

# read nucleosome library table
def read_table (fname):
    subunit_list = ['H2A', 'H2B', 'H3', 'H4']
    ID_BC, BC_ID = {}, {}
    ID_minfo = {}
    #ID_mutations = {}
    ID_shortname = {}
    First = True
    for line in open(fname):
        line = line.strip()
        if First:
            First = False
            continue
        if not line:
            continue
        cols = line.split('\t')
        ID, H2A, H2B, H3, H4, DNA, BC = cols
        ID = int(ID)
        assert ID not in ID_BC
        ID_BC[ID] = BC
        assert BC not in BC_ID
        BC_ID[BC] = ID

        assert ID not in ID_minfo
        ID_minfo[ID] = {}
        mutations = []
        shortname = ""
        for subunit, mhistone in zip(subunit_list, [H2A, H2B, H3, H4]):
            mhistone = mhistone.strip()
            ID_minfo[ID][subunit] = {}
            if mhistone == 'NA':
                ID_minfo[ID][subunit]['name'] = None
                ID_minfo[ID][subunit]['mutations'] = {}
            else:
                hname, pos_mutation = mhistone_parser(mhistone)
                ID_minfo[ID][subunit]['name'] = hname
                ID_minfo[ID][subunit]['mutations'] = pos_mutation

            if 'KpolyAc' in mhistone:
                mutation = subunit + 'KpolyAc'
            elif 'Acidic Patch Mutant' in mhistone:
                mutation = subunit + ' AP mutant'
            else:
                mutation = mhistone
            if mutation != subunit:
                mutations.append(mutation)

        if set(mutations) == set(['NA']):
            shortname += 'freeDNA'
        elif len(mutations) == 0:
            shortname += 'WT'
        else:
            shortname += '/'.join(mutations)

        if DNA == 'NA':
            ID_minfo[ID]['DNA'] = None
        else:
            ID_minfo[ID]['DNA'] = DNA
            shortname += ' (' + DNA + ')'

        assert ID not in ID_shortname
        ID_shortname[ID] = shortname

    return ID_BC, BC_ID, ID_minfo, ID_shortname
ID_BC, BC_ID, ID_minfo, ID_shortname = read_table('PTMlibTable.csv')

# read titration data file
def read_titration_data (fname):
    conc_list = []
    rep_fracts = [[], [], []]
    mean_list = []
    std_list = []

    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        conc, mean, std = float(cols[0]), float(cols[-2]), float(cols[-1])
        #fract_cols = cols[1:4]
        fract_cols = cols[4:4+3]
        conc_list.append(conc)
        mean_list.append(mean)
        std_list.append(std)
        for i in range(len(fract_cols)):
            rep_fracts[i].append(float(fract_cols[i]))

    return conc_list, rep_fracts, mean_list, std_list

# read sort file
def read_sort (fname,
               idex_name,
               BC_ID):
    vtype_count = {}
    name_ID_count = {}
    for line in open(fname):
        line = line.strip()
        if line.startswith('@'):
            vtype, detail = line.split('::')[1].split(':')

            if vtype == 'invalid':
                vtype += detail

            if vtype not in vtype_count:
                vtype_count[vtype] = 0
            vtype_count[vtype] +=1

            if vtype != 'valid':
                continue
            
            BC, index = detail[1:-1].split('][')
            
            try:
                name = index_name[index]
                ID = BC_ID[BC]
            except:
                continue

            if name not in name_ID_count:
                name_ID_count[name] = {}
            if ID not in name_ID_count[name]:
                name_ID_count[name][ID] = 0
            name_ID_count[name][ID] +=1
    return vtype_count, name_ID_count

# 4-parameter logistic function (sigmoid type)
def sigmoid_func (x, top, rate, chalf, bottom):
    y = bottom + float(top-bottom)/(1+np.exp(rate*(x-chalf)))
    return y

# 4-parameter logistic function (Hill type)
def hill_func (x, top, rate, chalf, bottom):
    y = bottom + float(top-bottom)/(1.0 + (x/float(chalf))**rate)
    return y

# objective function to optimize parameters
def obj_func (parms, func, X, Y):
    warnings.filterwarnings("ignore")
    Y_pred = func(X, *parms)
    return np.sum((Y_pred - Y)**2)

# compute CP value from 4PL data
def get_CP (top, hill, chalf, bottom, percent):
    surv_frac = 1 - percent/100.0
    CP = chalf*(((float(top-bottom)/(surv_frac-bottom)) - 1)**(1.0/float(hill)))
    return CP

### parameters
agent_fullname = {'sp':'Spermine(4+)',
                  'spd':'Spermidine(3+)',
                  'CoH':'Cobalt Hexammine(3+)',
                  'PEG':'PEG 8000',
                  'HP1a':'HP1 $\\alpha$'}

agent_unit = {'sp':'mM',
              'spd':'mM',
              'CoH':'mM',
              'PEG':'%',
              'HP1a':'$\\mu$M'}

repnum = 3
agent_list = ['sp', 'spd', 'CoH', 'PEG', 'HP1a']
# remove IDs already susceptible to PstI digest
# BCs with PstI site (ID 71, 78, 113) and free DNA with PstI site (ID 117)
all_IDs = sorted(ID_BC.keys())
all_good_IDs = sorted(list(set(all_IDs) - set([71, 78, 113, 117])))


### data fnames
# titration files
agent_tfname = {'sp':"PTMlib_sp_titration_corrected.csv", #nano-drop bg corrected
                'spd':"PTMlib_spd_titration_corrected.csv", #nano-drop bg corrected
                'CoH':"PTMlib_CoH_titration_corrected.csv", #nano-drop bg corrected
                'PEG':"PTMlib_PEG_titration.csv",
                'HP1a':"PTMlib_HP1a_titration.csv"}

# NGS index library file
index_fname = "PTMlib_NGS_information.csv"

# Nucleosome library table
table_fname = 'PTMlibTable.csv'

# sort files
sort_fnames = ['PTMlib_1rep.sort',
               'PTMlib_2rep.sort',
               'PTMlib_3rep.sort']



def Nuclibrary_analysis (sort_fnames,
                         table_fname,
                         index_fname,
                         agent_tfname,
                         QC_check,
                         model,
                         method,
                         min_rsq,
                         min_top,
                         max_top,
                         min_bottom,
                         max_bottom,
                         min_rate,
                         max_rate,
                         min_chalf,
                         max_chalf,
                         out_fname):

    ### read table
    ID_BC, BC_ID, ID_minfo, ID_shortname = read_table(table_fname)

    ### read index
    index_name, name_index, name_titr = read_index(index_fname)    

    ### read titration file
    agent_tinfo = {}
    for agent in agent_tfname:
        tfname = agent_tfname[agent]
        conc_list, rep_fracs, mean_list, std_list = read_titration_data_new(tfname)
        if agent not in agent_tinfo:
            agent_tinfo[agent] = {}
        agent_tinfo[agent]['conc'] = conc_list
        agent_tinfo[agent]['survival'] = mean_list

    ### read sort file and get scores
    rep_agent_ID_titration = []
    rep_agent_ID_score = []
    rep_agent_ID_logistic = []
    repnum = len(sort_fnames)
    for rep in range(repnum):
        sort_fname = sort_fnames[rep]

        print >> sys.stderr, "processing %s" % (sort_fname)

        ## read sort file
        vtype_count, name_ID_count = read_sort (sort_fname,
                                                index_name,
                                                BC_ID)

        ## get read counts for each conditions
        agent_tnum_ID_count = {}
        for name in name_ID_count:
            agent = name[:-1]
            tnum = int(name_titr[name]) - 1
            ID_count = name_ID_count[name]
            if agent not in agent_tnum_ID_count:
                agent_tnum_ID_count[agent] = {}
            agent_tnum_ID_count[agent][tnum] = ID_count

        agent_list = sorted(agent_tnum_ID_count.keys(), cmp=cmp_agent)

        ## QC sequencing
        if QC_check:
            # check data type
            vtypes = sorted(vtype_count.keys())
            counts = []
            for vtype in vtypes:
                counts.append(count)

            fig = plt.figure()
            plt.pie(counts,
                    labels=vtypes,
                    shadow=True,
                    startangle=90,
                    autopct='%1.1f%%')
            plt.savefig('sort_type_rep-%d.svg' % (rep+1),
                        format='svg',
                        bbox_inches='tight')
            plt.close()

            # check coverage by sample (agent, tnum)
            names, counts = [], []
            for agent in agent_list:
                tnum_ID_count = agent_tnum_ID_count[agent]
                tnums = sorted(tnum_ID_count.keys())
                for tnum in tnums:
                    name = '%s #%d' % (agent, tnum)
                    names.append(name)
                    count = sum(tnum_ID_count[tnum].values())
                    counts.append(count)

            fig = plt.figure()
            plt.bar(names, counts)
            plt.xticks(range(len(names)),
                       names,
                       ha="right",
                       rotation_mode="anchor",
                       rotation=70)
            plt.ylabel('Read counts')
            plt.savefig('coverage_by_sample_rep-%d.svg' % (rep+1),
                        format='svg',
                        bbox_inches='tight')
            plt.close()

            # check coverage by BC (input sample only)
            IDs = set([])
            for agent in agent_list:
                ID_count = agent_tnum_ID_count[agent][0]
                IDs |= set(ID_count.keys())
            IDs = sorted(list(IDs))
            ID_idx = {IDs[i]:i for i in range(len(IDs))}
            
            fig = plt.figure()
            for agent in agent_list:
                ID_count = agent_tnum_ID_count[agent][0]            
                X, Y = [], []
                for ID, count in ID_count.items():
                    X.append(ID_idx[ID])
                    Y.append(count)
                plt.plot(X, Y, '.-', alpha=0.8, label=agent)

            xticks, xticklabels = [], []
            for i in range(len(IDs)/5):
                xticks.append(5*i)
                xticklabels.append(IDs[5*i])

            plt.xticks(xticks,
                       xticklabels,
                       ha="right",
                       rotation_mode="anchor",
                       rotation=70)
            
            plt.xlabel('Library ID')
            plt.ylabel('Read counts')
            plt.title('Input read counts')
            plt.legend()
            plt.savefig('coverage_by_ID_rep-%d.svg' % (rep+1),
                        format='svg',
                        bbox_inches='tight')
            plt.close()
            
        ## compute the survival probability over titrations
        agent_ID_titration = {}
        for agent in agent_tnum_ID_count:
            conc_list = agent_tinfo[agent]['concentration']
            msurv_list = agent_tinfo[agent]['survival']
            
            tnum_ID_count = agent_tnum_ID_count[agent]
            tnums = sorted(tnum_ID_count.keys())
            assert tnums[0] == 0 # should include input control

            ID_control = tnum_ID_count[0]
            total_control = sum(ID_control.values())
            for tnum in tnums:
                ID_count = tnum_ID_count[tnum]
                total = sum(ID_count.values())

                conc = conc_list[tnum] # agent concentration
                msurv = msurv_list[tnum] # mean survival probability

                for ID in ID_control:
                    control = ID_control[ID]
                    if control <= 0: # skip when control is empty
                        continue

                    frac_control = float(control) / total_control
                    
                    count = ID_count[ID]
                    frac = float(count) / total

                    fold = frac / frac_control # get fold-change                    
                    surv = msurv * fold # get survival probability

                    if agent not in agent_ID_titration:
                        agent_ID_titration[agent] = {}
                    if ID not in agent_ID_titration[agent]:
                        agent_ID_titration[agent][ID] = {'conc':[],
                                                         'foldchange':[],
                                                         'survival':[]}
                        
                    agent_ID_titration[agent][ID]['conc'].append(conc)
                    agent_ID_titration[agent][ID]['foldchange'].append(fold)
                    agent_ID_titration[agent][ID]['survival'].append(surv)

        ## compute the condensability metrics
        agent_ID_score = {}
        agent_ID_logistic = {}
        for agent in agent_list:
            for ID in all_good_IDs:
                X = agent_ID_titration[agent][ID]['conc']
                Y = agent_ID_titration[agent][ID]['survival']

                ## compute the condensability score
                score = np.mean(-np.log2(np.asarray(Y)))

                if agent not in agent_ID_score:
                    agent_ID_score[agent] = {}
                agent_ID_score[agent][ID] = score

                ## fitting survival curve with logistic function
                # set guess and boundary of parameters
                top_guess = max(Y)
                bottom_guess = min(Y)
                chalf_guess = np.mean(X)
                rate_guess = 1.0

                if min_top == None:
                    top_min = top_guess * (1-0.01)
                else:
                    top_min = min_top

                if max_top == None:
                    top_max = top_guess * (1+0.01)
                else:
                    top_max = max_top

                if min_bottom == None:
                    bottom_min = bottom_guess * (1-0.01)
                else:
                    bottom_min = min_bottom

                if max_bottom == None:
                    bottom_max = bottom_guess * (1+0.01)
                else:
                    bottom_max = max_bottom

                if min_chalf == None:
                    chalf_min = min(X)
                else:
                    chalf_min = min_chalf

                if max_chalf == None:
                    chalf_max = max(X)
                else:
                    chalf_max = max_chalf

                if min_rate == None:
                    rate_min = 0.0
                else:
                    rate_min = min_rate

                if max_rate == None:
                    rate_max = 100.0
                else:
                    rate_max = max_rate

                # fitting the data with a logistic function
                if method == 'curve_fit':
                    # set initial guess of parameters
                    p0 = [top_guess, rate_guess, chalf_guess, bottom_guess]
                    # set boundary values of parameters
                    bounds = ([top_min, rate_min, chalf_min, bottom_min],
                              [top_max, rate_max, chalf_max, bottom_max])

                    try:
                        if model == 'sigmoid':                        
                            p_opt, p_cov = curve_fit(sigmoid_func,
                                                     X,
                                                     Y,
                                                     p0,
                                                     bounds = bounds,
                                                     method='dogbox')
                            
                        elif model == 'hill':
                            p_opt, p_cov = curve_fit(hill_func,
                                                     X,
                                                     Y,
                                                     p0,
                                                     bounds = bounds,
                                                     method='dogbox')

                        success = True
                    except:
                        success = False

                elif method == 'evolution':
                    # set boundary values of parameters
                    bounds = [(top_min, top_max),
                              (rate_min, rate_max),
                              (chalf_min, chalf_max),
                              (bottom_min, bottom_max)]

                    if model == 'sigmoid':
                        result = differential_evolution(obj_func,
                                                        args=(sigmoid_func, X, Y),
                                                        bounds=bounds,
                                                        seed=3)
                    elif model == 'hill':
                        result = differential_evolution(obj_func,
                                                        args=(hill_func, X, Y),
                                                        bounds=bounds,
                                                        seed=3)

                    p_opt = result.x
                    success = result.success

                # check fitting quality
                if success:
                    if model == 'sigmoid':
                        residuals = np.asarray(Y)- sigmoid_func(X, *p_opt)
                    elif model == 'hill':
                        residuals = np.asarray(Y)- hill_func(X, *p_opt)

                    ss_res = np.sum(residuals**2)
                    ss_tot = np.sum((np.asarray(Y)-np.mean(Y))**2)
                    r_squared = 1 - (ss_res / ss_tot)

                    # too poor fitting
                    if r_squared < min_rsq:
                        success = False
                        pass

                if success:
                    top, rate, chalf, bottom = p_opt
                    rsq_list.append(r_squared)

                else:
                    top, rate, chalf, bottom = np.NaN, np.NaN, np.NaN, np.NaN
                    r_squared = np.NaN
                    fail_count +=1

                if agent not in agent_ID_logistic:
                    agent_ID_logistic[agent] = {}
                agent_ID_logistic[agent][ID] = {}
                agent_ID_logistic[agent][ID]['top'] = top
                agent_ID_logistic[agent][ID]['rate'] = rate
                agent_ID_logistic[agent][ID]['chalf'] = chalf
                agent_ID_logistic[agent][ID]['bottom'] = bottom

        rep_agent_ID_titration.append(agent_ID_titration)
        rep_agent_ID_score.append(agent_ID_score)
        rep_agent_ID_logistic.append(agent_ID_logistic)


    ### writing the output files
    print >> sys.stderr, "writing output"

    ## save survival probability over titration
    repnum = len(rep_agent_ID_titration)
    for rep in range(repnum):
        agent_ID_titration = rep_agent_ID_titration[rep]

        for agent in agent_list:
            f = open("PTMlib_%s_rep-%d_survival.txt" % (agent, rep+1), 'w')

            # write down the header
            s = ""
            s += "ID" + '\t'
            s += 'Name' + '\t'
            s += '\t'.join(['[%s] %f %s' %
                            (agent, value, agent_unit[agent])
                            for value in agent_ID_titration[agent][1]['conc']])
            print >> f, s

            # record library data
            for ID in ID_list:
                s = ""
                s += str(ID) + '\t'
                s += ID_shortname[ID] + '\t'
                Y = agent_ID_titration[agent][ID]['survival']
                s += '\t'.join([str(y) for y in Y])
                print >> f, s
            f.close()    

    ## save condensability scores
    repnum = len(rep_agent_ID_score)
    for agent in agent_list:
        f = open("PTMlib_%s_score.txt" % (agent), 'w')

        ## write down the header
        s = ""
        s += "ID" + '\t'
        s += 'Name' + '\t'
        s += '\t'.join(['Score (rep-%d)' % (rep+1) for rep in range(repnum)])
        s += '\t' + 'Score (Mean)'
        s += '\t' + 'Score (Std)'
        print >> f, s

        ## record library data
        for ID in ID_list:
            s = ""
            s += str(ID) + '\t'
            s += ID_shortname[ID] + '\t'
            scores = []
            for rep in range(repnum):
                score = rep_agent_ID_score[rep][agent][ID]
                scores.append(score)
            s += '\t'.join([str(score) for score in scores])
            s += '\t' + str(np.mean(scores))
            s += '\t' + str(np.std(scores))
            print >> f, s
        f.close()

     ## save logistic parameters
     repnum = len(rep_agent_ID_logistic)
     for agent in agent_list:
        f = open("PTMlib_%s_logistic.txt" % (agent), 'w')

        ## write down the header
        s = ""
        s += "ID" + '\t'
        s += 'Name' + '\t'
        s += '\t'.join(['Top (rep-%d)' % (rep) for rep in range(1, repnum+1)])
        s += '\t'.join(['Rate (rep-%d)' % (rep) for rep in range(1, repnum+1)])
        s += '\t'.join(['C-half (rep-%d)' % (rep) for rep in range(1, repnum+1)])
        s += '\t'.join(['Bottom (rep-%d)' % (rep) for rep in range(1, repnum+1)])
        s += '\t'.join(['Top (Mean)', 'Rate (Mean)', 'C-half (Mean)', 'Bottom (Mean)'])
        s += '\t'.join(['Top (Std)', 'Rate (Std)', 'C-half (Std)', 'Bottom (Std)'])
        print >> f, s

        ## record library data
        for ID in ID_list:
            s = ""
            s += str(ID) + '\t'
            s += ID_shortname[ID] + '\t'

            tops, rates, chalfs, bottoms = [], [], [], []
            for rep in range(repnum):
                top = rep_agent_ID_logistic[rep][agent][ID]['top']
                tops.append(top)
                rate = rep_agent_ID_logistic[rep][agent][ID]['rate']
                rates.append(rate)
                chalf = rep_agent_ID_logistic[rep][agent][ID]['chalf']
                chalfs.append(chalf)
                bottom = rep_agent_ID_logistic[rep][agent][ID]['bottom']
                bottoms.append(bottom)

            s += '\t'.join([str(top) for top in tops])
            s += '\t'.join([str(rate) for rate in rates])
            s += '\t'.join([str(chalf) for chalf in chalfs])
            s += '\t'.join([str(bottom) for bottom in bottoms])

            s += '\t' + str(np.mean(tops))
            s += '\t' + str(np.mean(rates))
            s += '\t' + str(np.mean(chalfs))
            s += '\t' + str(np.mean(bottoms))

            s += '\t' + str(np.std(tops))
            s += '\t' + str(np.std(rates))
            s += '\t' + str(np.std(chalfs))
            s += '\t' + str(np.std(bottoms))
            print >> f, s
            
        f.close()

    print >> sys.stderr, "writing output"
    
if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')        

    parser = ArgumentParser(description='Analyze Nucleosome library sort files')
    parser.add_argument(metavar='-f',
                        dest="sort_fnames",
                        type=str,
                        nargs='+',
                        help='sort files')
    parser.add_argument('--table',
                        dest='table_fname',
                        type=str,
                        help='Nucleosome library table file')
    parser.add_argument('--index',
                        dest='index_fname',
                        type=str,
                        help='NGS index file')
    parser.add_argument('-t',
                        dest='tfnames',
                        type=str,
                        nargs='+',
                        help='titration files (agent:filename)')
    parser.add_argument('--exclude',
                        dest='exclude_IDs',
                        type=str,
                        nargs='+',
                        help='exclude IDs for analysis')
    parser.add_argument('--QC',
                        dest="QC_check",
                        type=str2bool,
                        default=True,
                        help='check sequencing QC')
    parser.add_argument('--model',
                        dest="model",
                        type=str,
                        default='sigmoid',
                        help='logistic model for fitting data (sigmoid or hill)')
    parser.add_argument('--method',
                        dest="method",
                        type=str,
                        default='evolution',
                        help='logistic regression method (curve_fit or evolution)')    
    parser.add_argument('--min_rsq',
                        dest="min_rsq",
                        type=float,
                        default=0.5,
                        help='minimum R-squared value for fitting quality')
    parser.add_argument('--min_top',
                        dest="min_top",
                        type=float,
                        help='lower bound of Top parameter in 4PL model')
    parser.add_argument('--max_top',
                        dest="max_top",
                        type=float,
                        help='upper bound of Top parameter in 4PL model')
    parser.add_argument('--min_bottom',
                        dest="min_bottom",
                        type=float,
                        help='lower bound of Bottom parameter in 4PL model')
    parser.add_argument('--max_bottom',
                        dest="max_bottom",
                        type=float,
                        help='upper bound of Bottom parameter in 4PL model')
    parser.add_argument('--min_rate',
                        dest="min_rate",
                        type=float,
                        help='lower bound of Rate parameter in 4PL model')
    parser.add_argument('--max_rate',
                        dest="max_rate",
                        type=float,
                        help='upper bound of Rate parameter in 4PL model')
    parser.add_argument('--min_chalf',
                        dest="min_chalf",
                        type=float,
                        help='lower bound of C-half parameter in 4PL model')
    parser.add_argument('--max_chalf',
                        dest="max_chalf",
                        type=float,
                        help='upper bound of C-half parameter in 4PL model')
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()

    # set titration files
    agent_tfname = {}
    for pair in args.tfnames:
        agent, tfname = pair
        agent_tfname[agent] = tfname

    # exclude IDs
    exclude_IDs = []
    for ID in args.exclude_IDs:
        try:
            ID = int(ID)
        except:
            pass
        exclude_IDs.append(ID)

    # set logistic model
    model = args.model.lower()

    # set regression method
    method = args.method.lower()

    Nuclibrary_analysis (sort_fnames,
                         table_fname,
                         index_fname,
                         agent_tfname,
                         exclude_IDs,
                         QC_check,
                         model,
                         method,
                         args.min_rsq,
                         args.min_top,
                         args.max_top,
                         args.min_bottom,
                         args.max_bottom,
                         args.min_rate,
                         args.max_rate,
                         args.min_chalf,
                         args.max_chalf,
                         args.out_fname
                         )
