import random
import analysis

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

def GC_content(seq):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    return (num/float(len(seq)))*100

def Amer_len(seq):
    num = []
    i = 0
    while i < len(seq):
        if seq[i] == 'A':
            count = 1
            j = i + 1
            while j < len(seq):
                if seq[j] != 'A':
                    break
                count +=1
                j +=1
            num.append(count)
            i = j + 1
        else:
            i +=1
    if len(num) == 0:
        return 0
    return max(num)


# select the subset of inserts
def sampling (cutlocs_list, dyadlocs_list, num_choice, sample_mode):
    assert len(cutlocs_list) == len(dyadlocs_list)
    sample_list = []
    
    for i in range(len(cutlocs_list)):
        cutlocs = cutlocs_list[i]
        dyadlocs = dyadlocs_list[i]
        assert len(cutlocs) == len(dyadlocs)
        allint_list = cutlocs.keys()
        int_list = []

        
        if num_choice > len(allint_list):
            num_choice = len(allint_list)

        if sample_mode == 'r':
            int_list = random.sample(allint_list, num_choice)

        elif sample_mode.startswith('s'):
            input_list = sample_mode.split(':')[1].strip().split(',')
            for insert in input_list:
                insert = insert.upper()
                if insert in allint_list:
                    int_list.append(insert)
                else:
                    print "insert " + insert +" is missing"
            num_choice = len(int_list)

        elif sample_mode.startswith('cluster'):
            index, num = sample_mode.split(":")[1].split('/')
            index, num = int(index), int(num)
            assert index > 0 and index <= num
            cutmaps_list, dyadmaps_list = analysis.get_maps([cutlocs], [dyadlocs], norm_choice=True)
            insert_cdx_list, cdx_insert_list = analysis.Kmeans(dyadmaps_list, num)
            cdx_insert = cdx_insert_list[0]
            if num_choice < len(cdx_insert[index-1]):
                int_list = random.sample(cdx_insert[index-1], num_choice)
            else:
                int_list = cdx_insert[index-1]
            
        else:
            mode, rank = sample_mode.split(':')
            rank = int(rank) -1
            st, en = rank*num_choice, (rank+1)*num_choice
            if st >= len(allint_list):
                st, en = len(allint_list)-num_choice, len(allint_list)
            elif en > len(allint_list):
                en = len(allint_list)

            if mode in ['t','b']:
                insert_count = {}
                for insert in dyadlocs:
                    insert_count[insert] = len(dyadlocs[insert])
                count_insert = [[count, insert] for insert, count in insert_count.items()]
                if mode == 't':
                    count_insert = sorted(count_insert, cmp=key_cmp, reverse=True)[st:en]
                else:
                    count_insert = sorted(count_insert, cmp=key_cmp)[st:en]
                for count, insert in count_insert:
                    int_list.append(insert)

            elif mode in ['GCt', 'GCb']:
                insert_GC = {}
                for insert in dyadlocs:
                    insert_GC[insert] = GC_content(insert)
                GC_insert = [[GC, insert] for insert, GC in insert_GC.items()]
                if mode == 'GCt':
                    GC_insert = sorted(GC_insert, cmp=key_cmp, reverse=True)[st:en]
                else:
                    GC_insert = sorted(GC_insert, cmp=key_cmp)[st:en]
                for GC, insert in GC_insert:
                    int_list.append(insert)

            elif mode in ['at','ab']:
                insert_Amer = {}
                for insert in dyadlocs:
                    insert_Amer[insert] = Amer_len(insert)
                Amer_insert = [[Amer, insert] for insert, Amer in insert_Amer.items()]
                if mode == 'at':
                    Amer_insert = sorted(Amer_insert, cmp=key_cmp, reverse=True)[st:en]
                else:
                    Amer_insert = sorted(Amer_insert, cmp=key_cmp)[st:en]
                for Amer, insert in Amer_insert:
                    int_list.append(insert)

            elif mode == 'sort':
                int_list = allint_list
                int_list = sorted(int_list)[st:en]

            else:
                assert mode == 'rsort'
                int_list = allint_list
                int_list = sorted(int_list, reverse=True)[st:en]

        sample_list.append(int_list)
        for int_list in sample_list:
            print "cond" + str(i+1) + " selected: " + str(len(int_list)) 

    return sample_list
