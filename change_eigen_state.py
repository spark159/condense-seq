from functools import cmp_to_key

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
        if not line.startswith("chr"):
            continue
        cols = line.split()
        chr, start, end, value = cols
        start = int(start)
        end = int(end)
        value = float(value)
        if chr not in chr_win_value:
            chr_win_value[chr] = {}
        win = (start, end)
        chr_win_value[chr][win] = value
    return chr_win_value

## change eigen to state
fname = "eigen_H1_100kb.bedgraph"
chr_win_value = read_bedgraph(fname)

chr_win_state = {}
for chr in chr_win_value:
    for win in chr_win_value[chr]:
        value = chr_win_value[chr][win]
        if value >=0:
            state = 0
        else:
            state = 1
        if chr not in chr_win_state:
            chr_win_state[chr] = {}
        chr_win_state[chr][win] = state


## write HMM annotation
def write_line (sinterval, f):
    state, chr, start, end = sinterval
    s = '%s\t%d\t%d\t%s' % (chr, start, end, 'E'+str(state+1))
    print (s, end='\n', file=f)

f = open(fname.rsplit('.', 1)[0] + '_state.bed', 'w')
prev_state = None
sinterval = []
chr_list = sorted(chr_win_state.keys(), key=cmp_to_key(chr_cmp))
for chr in chr_list:
    win_state = chr_win_state[chr]
    for win in sorted(win_state):
        start, end = win
        state = win_state[win]
        if prev_state == None:
            sinterval = [state, chr, start]
            prev_state = state
            prev_chr = chr
            prev_end = end
        if prev_state != state or prev_chr != chr:
            sinterval.append(prev_end)
            write_line (sinterval, f)
            sinterval = [state, chr, start]
        prev_state = state
        prev_chr = chr
        prev_end = end
sinterval.append(prev_end)
write_line (sinterval, f)
f.close()
