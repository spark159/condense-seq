import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import copy
import pickle

def mean(data_list):
    assert len(data_list) > 0
    return float(sum(data_list))/len(data_list)

def std(data_list):
    assert len(data_list) > 1
    total = 0.0
    mean = mean(data_list)
    for data in data_list:
        total += (data - mean)**2
    return math.sqrt(float(total)/(len(data_list)-1))

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

def GC_content(seq):
    seq = seq.upper()
    output = 0.0
    for nt in seq:
        if nt in "GC":
            output += 1.0
        elif nt in 'N':
            output += 0.5
    return output/len(seq)

#Turn Phred+33 ASCII-encoded quality into Phred-scaled integer
def phred33_to_q(qual):
    return ord(qual)-33

#Turn Phred-scaled integer into Phred+33 ASCII-encoded quality
def q_to_phred33(Q):
    return chr(Q + 33)

#Turn Phred-scaled integer into error probability
def q_to_p(Q):
    return 10.0 ** (-0.1 * Q)

#Turn error probability into Phred-scaled integer
def p_to_q(p):
    return int(round(-10.0 * math.log10(p)))

def mean_prob_error (quality):
    total = 0.0
    for qual in quality:
        total += q_to_p(phred33_to_q(qual))
    return total/len(quality)

def read_QC (fnames,
             out_fname):

    # read the fastq file
    labels = []
    pt_key = {1:'seq', 2:'opt', 3:'quality'}
    for i in range(len(fnames)):
        filename = fnames[i]
        label = filename.rsplit('/')[-1]
        labels.append(label)
        print >> sys.stderr, "reading %s" % (label)

        exten = filename.rsplit('.',2)
        if exten[-1] == 'fastq':
            cmd = ["cat", filename]
        elif exten[-1] == 'gz' and exten[-2] == 'fastq':
            cmd = ["gunzip",  "-c", filename]
        else:
            print >> sys.stderr, "Error: %s is not fastq file." % (fname)
            sys.exit(1)
            
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=open("/dev/null", 'w'))

        # make qc data list
        # per position
        pos_nt_count = {} # nts per position
        pos_probs= {} # error probabilities per position
        # per read
        length_list = [] # read length list
        GC_list = [] # GC content list
        prob_list = [] # mean error probablity list
        
        line_pt = None
        for line in proc.stdout:
            if line.startswith('@'):
                read_id=':'.join(re.split(':|/', line)[3:5])
                if not read_id:
                    continue
                line_pt = 0
                continue

            if line_pt == None:
                continue

            line_pt += 1
            key = pt_key[line_pt]

            if key == "seq":
                seq = line.strip()
                GC_list.append(GC_content(seq))
                length_list.append(len(seq))

                for k in range(len(seq)):
                    if k not in pos_nt_count:
                        pos_nt_count[k] = {}
                    nt = seq[k].upper()
                    if nt not in 'ATCG':
                        continue
                    if nt not in pos_nt_count[k]:
                        pos_nt_count[k][nt] = 0
                    pos_nt_count[k][nt] += 1

            elif key == "quality":
                quality = line.strip()
                prob_list.append(mean_prob_error(quality))

                for k in range(len(quality)):
                    if k not in pos_probs:
                        pos_probs[k] = []
                    prob = q_to_p(phred33_to_q(quality[k]))
                    pos_probs[k].append(prob)

            else:
                 continue

        #temporal pickling
        f1 = open(filename + "_length_list_" + str(i) + ".pickle", "wb")
        pickle.dump(length_list, f1)
        f1.close()
        f2 = open(filename + "_GC_list_" + str(i) + ".pickle", "wb")
        pickle.dump(GC_list, f2)
        f2.close()
        f3 = open(filename + "_prob_list_" + str(i) + ".pickle", "wb")
        pickle.dump(prob_list, f3)
        f3.close()
        f4 = open(filename + "_pos_nt_count_" + str(i) + ".pickle", "wb")
        pickle.dump(pos_nt_count, f4)
        f4.close()
        f5 = open(filename + "_pos_probs_" + str(i) + ".pickle", "wb")
        pickle.dump(pos_probs, f5)
        f5.close()


    # writing output file
    print >> sys.stderr, "writing qc file"

    f = open(out_fname + '_qc.txt', 'w')

    s = "@Input:"
    for name in labels:
        s += name + ','
    s = s[:-1]
    print >> f, s

    print >> f, "@length"
    for i in range(len(fnames)):
        f1 = open(fnames[i] + "_length_list_" + str(i) + ".pickle", "rb")
        length_list = pickle.load(f1)
        print >> f, ",".join([str(value) for value in length_list])
        f1.close()
        os.remove(fnames[i] + "_length_list_" + str(i) + ".pickle")
        
    print >> f, "@GCcontent"
    for i in range(len(fnames)):
        f2 = open(fnames[i] + "_GC_list_" + str(i) + ".pickle", "rb")
        GC_list = pickle.load(f2)
        print >> f, ",".join([str(value) for value in GC_list])
        f2.close()
        os.remove(fnames[i] + "_GC_list_" + str(i) + ".pickle")
        
    print >> f, "@MeanErrorProb"
    for i in range(len(fnames)):
        f3 = open(fnames[i] + "_prob_list_" + str(i) + ".pickle", "rb")
        prob_list = pickle.load(f3)
        print >> f, ",".join([str(round(value, 5)) for value in prob_list])
        f3.close()
        os.remove(fnames[i] + "_prob_list_" + str(i) + ".pickle")

    total_pos_nt_count = {}
    for i in range(len(fnames)):
        f4 = open(fnames[i] + "_pos_nt_count_" + str(i) + ".pickle", "rb")
        pos_nt_count = pickle.load(f4)
        for pos in pos_nt_count:
            for nt in pos_nt_count[pos]:
                count = pos_nt_count[pos][nt]
                if pos not in total_pos_nt_count:
                    total_pos_nt_count[pos] = {}
                if nt not in total_pos_nt_count[pos]:
                    total_pos_nt_count[pos][nt] = 0
                total_pos_nt_count[pos][nt] += 1
        f4.close()
        os.remove(fnames[i] + "_pos_nt_count_" + str(i) + ".pickle")

    max_len = max(total_pos_nt_count.keys())
    nt_fracs = {nt:[0.0]*max_len for nt in 'ATCG'}
    for pos in xrange(max_len):
        try:
            total = sum(total_pos_nt_count[pos].values())
        except:
            continue
        for nt in 'ATCG':
            try:
                count = total_pos_nt_count[pos][nt]
            except:
                count = 0.0
            frac = float(count)/total
            nt_fracs[nt][pos] = round(frac, 5)            

    for nt in "ATCG":
        print >> f, "@%sFractionPerPosition" % (nt)
        print >> f, ",".join([str(value) for value in nt_fracs[nt]])

    pos_total, pos_dcount = [], []
    for i in range(len(fnames)):
        f5 = open(fnames[i] + "_pos_probs_" + str(i) + ".pickle", "rb")
        for pos in xrange(max(pos_probs.keys())):
            try:
                total = sum(pos_probs[pos])
                dcount = len(pos_probs[pos])
            except:
                total = 0.0
                dcount = 0
            if pos >= len(pos_total):
                pos_total.append(0.0)
            if pos >= len(pos_dcount):
                pos_dcount.append(0)
            pos_total[pos] += total
            pos_dcount[pos] += dcount
        f5.close()
        os.remove(fnames[i] + "_pos_probs_" + str(i) + ".pickle")

    mean_probs = []
    for i in range(len(pos_total)):
        total = pos_total[i]
        dcount = pos_dcount[i]
        if total <= 0:
            mean_probs.append('NA')
            continue
        mean = float(total) / dcount
        mean_probs.append(mean)

    print >> f, "@MeanErrorProbPerPosition"
    print >> f, ",".join([str(round(value, 5)) for value in mean_probs])

    f.close()

    print >> sys.stderr, "Done"
    

if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = ArgumentParser(description='Quality control on read files')
    parser.add_argument(metavar='-f',
                        dest="fnames",
                        type=str,
                        nargs='+',
                        help='fastq filenames')
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()

    if not args.out_fname:
        def common (strings):
            end = min([len(string) for string in strings])
            if end <= 0:
                return None
            while end > 0:
                check = strings[0][:end]
                find = True
                for i in range(1, len(strings)):
                    part = strings[i][:end]
                    if check != part:
                        find = False
                        break
                if find:
                    break
                end -= 1
            if find:
                return check
            else:
                return None            
        common_name = common([fname.rsplit('.',1)[0] for fname in args.sort_filenames])

        if common_name != None:
            out_fname = common_name
        else:
            out_fname = 'out'
    else:
        out_fname = args.out_fname

    read_QC (args.fnames,
             out_fname)
