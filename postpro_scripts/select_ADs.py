import Interval_dict

def chr_cmp (chr1, chr2):
    try:
        chrnum1 = int(chr1[3:])
    except:
        chrnum1 = chr1[3:]
    try:
        chrnum2 = int(chr2[3:])
    except:
        chrnum2 = chr2[3:]
    if type(chrnum1) == type(chrnum2):
        if chrnum1 < chrnum2:
            return -1
        elif chrnum1 > chrnum2:
            return 1
        else:
            return 0
    else:
        if type(chrnum1) != int:
            return 1
        return -1



def read_bedgraph(fname):
    chr_win_value = {}
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        if not line.startswith('chr'):
            continue
        cols = line.split()
        chr, start, end, value = cols
        start = int(start)
        end = int(end)
        try:
            value = float(value)
        except:
            value = value
        if chr not in chr_win_value:
            chr_win_value[chr] = {}
        win = (start, end)
        chr_win_value[chr][win] = value
    return chr_win_value

def read_chromHMM(fname, chr_choice=None, change=False):
    chr_state_intervals = {}
    for line in open(fname):
        cols = line.strip().split()
        chr, st, ed, state = cols[:4]
        if chr_choice!= None and chr not in chr_choice:
            continue
        st, ed = int(st), int(ed)
        if change:
            state = change[state]
        if chr not in chr_state_intervals:
            chr_state_intervals[chr] = {}
        if state not in chr_state_intervals[chr]:
            chr_state_intervals[chr][state] = []
        chr_state_intervals[chr][state].append((st,ed))
    return chr_state_intervals



# parameters
#path = "/home/spark159/../../storage/"
path = "./data/"

# bedgraph file information
name_bedgraph = {"LAD":"Lamin_score.bedgraph",
                 "SPAD":"SON_score.bedgraph",
                 "NAD":"K562_NucleolarDamID.bedgraph"}

# HMM bed file information
name_bed = {"LAD":"LAD_Gaussian_2.bed",
            "SPAD":"SPAD_Gaussian_2.bed",
            "NAD":"NAD_Gaussian_2.bed"}

# change HMM state name
name_state = {"SPAD":'E2',
              "NAD":'E1',
              "LAD":'E1'}

# cut off
name_cutoff = {"SPAD":0.03,
               "NAD":0.3,
               "LAD":0.3}


# set names to be analyzed
names = ['SPAD', 'NAD', 'LAD']


for name in names:
    chr_win_value = read_bedgraph(path+name_bedgraph[name])
    chr_state_intervals = read_chromHMM(name_bed[name])

    chr_good_sIDs = {}
    chr_sID_interval = {}
    for chr in chr_win_value:
        state_intervals = chr_state_intervals[chr]

        sID_interval = {}
        for state in state_intervals:
            intervals = state_intervals[state]
            for k in range(len(intervals)):
                sID = (state, k)
                sID_interval[sID] = intervals[k]

        state_dict = Interval_dict.double_hash(sID_interval,
                                               domain_size=10000,
                                               max_pos=5*10**8)

        win_value = chr_win_value[chr]
        for win, value in win_value.items():
            st, ed = win
            state_dict.insert_range(st, ed, value)
        sID_total = state_dict.get()

        sID_score = {}
        for sID in sID_interval:
            try:
                interval = sID_interval[sID]
                score = float(sID_total[sID])/(interval[1] - interval[0])
            except:
                score = 0.0
            sID_score[sID] = score

        score_sID = sorted([(score, sID) for sID, score in sID_score.items()], reverse=True)
        score_sID = score_sID[:int(len(score_sID)*name_cutoff[name])]
        good_sIDs = [sID for _, sID in score_sID]

        chr_good_sIDs[chr] = good_sIDs
        chr_sID_interval[chr] = sID_interval

    f = open(name + '_selected.bed', 'w')
    for chr in sorted(chr_good_sIDs, cmp=chr_cmp):
        good_sIDs = chr_good_sIDs[chr]
        sID_interval = chr_sID_interval[chr]
        
        
        for sID in sorted(good_sIDs):
            state, _ = sID
            start, end = sID_interval[sID]
            print >> f, "%s\t%d\t%d\t%s" % (chr, start, end, state)
    f.close()
        
        
