import os, sys, subprocess, re
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy
import numpy as np

# IUPAC code
code_nts = {'A':['A'], 'T':['T'], 'C':['C'], 'G':['G'], 'R':['A','G'],
            'Y':['C','T'], 'S':['G','C'], 'W':['A','T'], 'K':['G','T'], 'M':['A','C'],
            'B':['C', 'G', 'T'], 'D':['A', 'G', 'T'], 'H':['A', 'C', 'T'], 'V':['A', 'C', 'G'],
            'N':['A', 'T', 'C', 'G']}

class GenomeCutter:
    def __init__(self,
                 genome_fname,
                 enzyme_fname,
                 chr_choice='all'):

        self.genome_fname = genome_fname
        self.enzyme_fname = enzyme_fname
        self.chr_choice = chr_choice
        self.enzyme_chr_cutsites = {}

        # read geneome sequence
        self.chr_seq = {}
        for line in open(self.genome_fname):
            line = line.strip()
            if line.startswith('>'):
                chr = line[1:]
                if (self.chr_choice != 'all') and chr not in self.chr_choice:
                    continue
                if chr not in self.chr_seq:
                    self.chr_seq[chr] = ""
                continue
            try:
                self.chr_seq[chr] += line
            except:
                continue
            
        # read enzyme information
        self.enzyme_dict = {}
        First = True
        for line in open(self.enzyme_fname):
            if First:
                First = False
                continue
            cols = line.strip().split()
            recogsite, enzymes = cols[0], cols[1:]

            cutseq = ""
            for substring in re.split('(\(-?\d+/-?\d+\))',recogsite):
                if substring == '':
                    continue
                if substring.startswith('('):
                    top, bottom = substring[1:-1].split('/')
                    top, bottom = int(top), int(bottom)
                    continue
                else:
                    cutseq += ''.join(substring.split('/'))

            for enzyme in enzymes:
                enzyme = enzyme.split('\xc2\xae')[0]
                self.enzyme_dict[enzyme] = {}
                self.enzyme_dict[enzyme]['cutseq'] = cutseq

    def cut (self,
             enzymes):
        
        # build bowtie2 index files
        if self.chr_choice != 'all':
            ref_fname = self.genome_fname.rsplit('.', 1)[0] + "_sub"
            f = open(subgenome_fname + '.fa', 'w')
            for chr in chr_choice:
                print >> f, chr
                print >> f, chr_seq[chr]
            f.close()
        else:
            ref_fname = self.genome_fname.rsplit('.', 1)[0]

        subprocess.call(["bowtie2-build", ref_fname+'.fa', ref_fname], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

        
        # make a fastq file with restriction enzyme information
        def all_path (states_list):
            N = len(states_list)
            if N==1:
                return states_list[0]
            output = []
            for path in all_path(states_list[:-1]):
                for state in states_list[-1]:
                    output.append(path+state)
            return output

        read_fname = "enzyme_list"
        f = open(read_fname + '.fastq', 'w')

        for enzyme in enzymes:
            states_list = [code_nts[code] for code in self.enzyme_dict[enzyme]['cutseq']]
            all_cutseqs = all_path(states_list)
            for i in range(len(all_cutseqs)):
                read_id = enzyme + ':' + str(i)
                read_seq = all_cutseqs[i]
                print >> f, '@' + read_id 
                print >> f, read_seq
                print >> f, '+'
                print >> f, 'G'*len(read_seq)

        f.close()

        
        # find out restriction sites by Bowtie2 alignment
        aligner_cmd = ["bowtie2", '-x', ref_fname, '-U', read_fname+'.fastq' ]
        aligner_cmd += ['--score-min', 'L,' + str(0) + ',' +str(0)] # consider only exact matches
        aligner_cmd += ['-N', str(1), '-L', str(1), '-i',  'S,1,0'] # turn off heuristic seeding
        aligner_cmd += ['--np', '1'] # ignore alignments with 'N' 
        aligner_cmd += ['-a'] # report all alignments
        align_proc = subprocess.Popen(aligner_cmd, stdout=subprocess.PIPE, stderr=open("/dev/null", 'w'))        
        enzyme_chr_cutsites = {}
        chr_cutsites = {}
        for line in align_proc.stdout:
            if line.startswith('@'):
                continue
            #print line
            cols = line.strip().split()
            read_id, flag, ref_id, pos, mapQ, cigar_str = cols[:6]
            flag, pos = int(flag), int(pos)
            pos-=1

            # invalid: mapping failure
            if pos < 0 or flag & 0x4 != 0:
                print >> sys.stderr, "Error: mapping failure."
                continue

            AS,NM,MD = None, None, None
            for i in range(11, len(cols)):
                col = cols[i]
                if col.startswith('AS'):
                    AS = int(col[5:])
                elif col.startswith('NM'):
                    NM = int(col[5:])
                elif col.startswith('MD'):
                    MD = col[5:]

            # Not exact match
            if NM > 0:
                print >> sys.stderr, "Error: not exact match."
                continue

            enzyme = read_id.split(':')[0]
            chr = ref_id.strip()    
 
            if enzyme not in enzyme_chr_cutsites:
                enzyme_chr_cutsites[enzyme] = {}
            if chr not in enzyme_chr_cutsites[enzyme]:
                enzyme_chr_cutsites[enzyme][chr] = []
            enzyme_chr_cutsites[enzyme][chr].append(pos)

            if chr not in chr_cutsites:
                chr_cutsites[chr] = []
            chr_cutsites[chr].append(pos)
            

        # remove intermediate files
        if self.chr_choice != 'all':
            subprocess.call(['rm', ref_fname + '.fa'])
        subprocess.call('rm ' + ref_fname + '*bt2', shell=True)
        subprocess.call(['rm', read_fname+'.fastq'])

        # update cutsite information
        self.enzyme_chr_cutsites.update(enzyme_chr_cutsites)

        # make fragment list
        all_fragments = []
        for chr, cutsites in chr_cutsites.items():
            cutsites = sorted(cutsites)
            fragments = []
            for i in range(len(cutsites)-1):
                fragments.append(cutsites[i+1] - cutsites[i])
            fragments.append(cutsites[0])
            fragments.append(len(self.chr_seq[chr]) - cutsites[-1])
            all_fragments += fragments
                
        return all_fragments


test = GenomeCutter("sacCer3.fa", "NEB_enzymes.csv")
#test = GenomeCutter("hg38.fa", "NEB_enzymes.csv")
fragments = test.cut(['BbvCI'])
