letter_code = {'B':['C','G','T'], 'D':['A','G','T'], 'H':['A','C','T'], 'K':['G','T'], 'M':['A','C'],'N':['A','C','G','T'],
               'R':['A','G'],'S':['C','G'],'V':['A','C','G'], 'W':['A','T'], 'Y':['C','T']}

def key_cmp(a, b):
    if a[0] <= b[0]:
        return -1
    else:
        return 1

def value_cmp(a, b):
    if a[1] <= b[1]:
        return -1
    else:
        return 1

def read_enzyme(filename):
    dic = {}
    for line in open(filename):
        if line.strip():
            enzyme, site = line.strip().split()
            site_list = [""]
            for i in range(len(site)):
                nt = site[i]
                if nt in "ATCG":
                    for i in range(len(site_list)):
                        site_list[i] += nt
                else:
                    new_site_list = []
                    alts = letter_code[nt]
                    for i in range(len(site_list)):
                        for next in alts:
                            new_site_list.append(site_list[i] + next)
                    site_list = new_site_list
            assert enzyme not in dic
            dic[enzyme] = site_list
    return dic

def read_ref(ref_fname):
    id_seq = {}
    for line in open(ref_fname):
        line = line.strip()
        if line.startswith('>'):
            id = int(line[4:])
            continue
        if line:
            assert id not in id_seq
            id_seq[id] = line
    return id_seq

enzyme_sites = read_enzyme("MeReEnzyme.txt")
id_seq = read_ref("Blib.ref")


enzyme_count = {}
for id, seq in id_seq.items():
    for enzyme, sites in enzyme_sites.items():
        for site in sites:
            if site in seq:
                if enzyme not in enzyme_count:
                    enzyme_count[enzyme] = 0
                enzyme_count[enzyme] +=1
                break
            
enzyme_count = [[enzyme, count] for enzyme, count in enzyme_count.items()]
enzyme_count = sorted(enzyme_count, cmp=value_cmp, reverse=True)

for enzyme, count in enzyme_count:
    print enzyme, count
