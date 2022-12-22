# read whole small genome reference fasta file
# fname: fasta file name
# output: dictionary of {chromosome name (key): chromosome sequence (value)}
def read_genome(fname):
    chr_dic = {}
    chr_name, sequence = "", ""
    for line in open(fname):
        if line.startswith(">"):
            if chr_name and sequence:
                chr_dic[chr_name] = sequence
            chr_name = line.strip().split()[0][1:]
            sequence = ""
        else:
            sequence += line.strip()
    if chr_name and sequence:
        chr_dic[chr_name] = sequence
    return chr_dic


# read reference file and extract sequence for each NCP
# ref_fname: name of fast file
# chr: chromosome (eg. chr1, chr2)
# st_ID: dictionary of {start location of nucleosome (key): nucleosome ID (value)}
# win: nucleosome size (eg. 145 bp)
# output: dictionary of {nucleosome ID (key): corresponding nucleosome sequence}
def get_seq(ref_fname, chr, st_ID, win):
    seq = ""
    pt = -1
    k = 0
    left = []
    Find = False
    stID = [[st, st_ID[st]] for st in sorted(st_ID.keys())]
    pos, ID = stID[k]
    ID_seq = {}
    for line in open(ref_fname):
        line = line.strip()
        if line.startswith(">"):
            if Find:
                break
            if line[1:] == chr:
                Find = True
            continue
        if Find:
            if len(left) == 0 and pt + len(line) < pos:
                pt += len(line)
                continue
            for i in range(len(left)):
                leftID, seq = left.pop(0)
                ed = min(len(line), win-len(seq))
                seq += line[:ed]
                if len(seq) == win:
                    #AT = AT_content(seq)
                    #ID_AT[leftID] = AT
                    ID_seq[leftID] = seq
                else:
                    left.append([leftID, seq])
            while pt + len(line) >= pos and k < len(stID):
                loc = pos - pt - 1
                seq = line[loc:min(loc+win,len(line))]
                if len(seq) == win:
                    #AT = AT_content(seq)
                    #ID_AT[ID] = AT
                    ID_seq[ID] = seq
                else:
                    left.append([ID, seq])
                k += 1
                try:
                    pos, ID = stID[k]
                except:
                    None
            if len(left) == 0 and len(ID_seq) == len(stID):
                break
            pt += len(line)
    while len(left) > 0:
        leftID, seq = left.pop(0)
        #AT = AT_content(seq)
        #ID_AT[leftID] = AT
        ID_seq[leftID] = seq
    assert len(ID_seq) == len(stID)
    return ID_seq
