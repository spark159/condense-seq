def AT_content(seq):
    count = 0.0
    for nt in seq:
        if nt in "AT":
            count += 1
    return count/len(seq)

def pick (num, data):
    data = sorted(data)
    return

def read_lib (fname):
    id_seq = {}
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            id = line[1:]
            continue
        assert id not in id_seq
        id_seq[id] = line.strip()
    return id_seq

id_seq = read_lib("Blib.ref")

AT_seq = {}
for id, seq in id_seq.items():
    AT = AT_content(seq)
    if AT not in AT_seq:
        AT_seq[AT] = []
    AT_seq[AT].append(id)

