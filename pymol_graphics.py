import os, sys, subprocess, re
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import random  
import math
import copy
import pickle
import pymol
import pytms

AA_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
           'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
           'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
           'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def read_FASTA_file (fname):
    id_seq = {}
    for line in open(fname):
        if not line.strip():
            continue
        if line.startswith('>'):
            id = line[1:].strip()
        else:
            seq = line.strip()
            id_seq[id] = seq
    return id_seq

def read_FASTA_string (s):
    id_seq = {}
    for line in open(fname):
        if not line.strip():
            continue
        if line.startswith('>'):
            id = line[1:].strip()
        else:
            seq = line.strip()
            id_seq[id] = seq
    return id_seq



class Molecules:
    def __init__(self, code):
        self.code = code
        self.chain_resi_resn = {}
        self.chain_resi_index_atom = {}

        # read data from pymol loaded struture
        #pymol.pymol_argv = ['pymol','-qc']
        #pymol.finish_launching()
        #pymol.cmd.log_open()
        pymol.cmd.fetch(code)
        
        for atom in pymol.cmd.get_model(code).atom:
            chain = atom.chain
            name = atom.name
            resn = atom.resn.upper()
            resi = int(atom.resi)
            x, y, z = atom.coord
            x, y, z = float(x), float(y), float(z)
            index = int(atom.index)
            b = float(atom.b)
            
            if chain not in self.chain_resi_resn:
                self.chain_resi_resn[chain] = {}
            if resi not in self.chain_resi_resn[chain]:
                self.chain_resi_resn[chain][resi] = resn

            if chain not in self.chain_resi_index_atom:
                self.chain_resi_index_atom[chain] = {}
            if resi not in self.chain_resi_index_atom[chain]:
                self.chain_resi_index_atom[chain][resi] = {}
            assert index not in self.chain_resi_index_atom[chain][resi]
            self.chain_resi_index_atom[chain][resi][index] = {"name":name,
                                                             "resn":resn,
                                                             "coord":(x, y, z),
                                                             "b-factor":b}
            
        self.chain_seq = {}
        for chain in self.chain_resi_resn:
            seq = ""
            for resi in sorted(self.chain_resi_resn[chain].keys()):
                resn = self.chain_resi_resn[chain][resi]
                if len(resn) == 3:
                    try:
                        resn = AA_dict[resn]
                    except:
                        continue
                        #resn = '[' + resn + ']'
                elif len(resn) == 2 and resn.startswith('D'):
                    resn = resn[1:]
                else:
                    continue
                    #print (resn)
                    #pass
                    #assert len(resn) == 1
                seq += resn
            self.chain_seq[chain] = seq

        print("%s is loaded" % (code), file=sys.stderr)


    def print_seq(self):
        for chain in sorted(self.chain_seq.keys()):
            print("chain %s" % (chain), file=sys.stderr)
            print(self.chain_seq[chain], file=sys.stderr)

    def selection_string (self,
                          chain_resi):
        selections = []
        for chain in chain_resi:
            select = []
            for resi in chain_resi[chain].keys():
                if resi < 0:
                    select.append('\\' + str(resi))
                else:
                    select.append(str(resi))
            select = ','.join(select)
            selections.append("(" + "chain " + chain + " and " + "resi " + select + ")")
        selections = " or ".join(selections)
        return selections
        

    def remove_ions(self):
        # remove waters and ions
        pymol.cmd.remove('resn hoh')
        pymol.cmd.remove('resn mn')
        pymol.cmd.remove('resn cl')
        return 
        

    def stylish(self):
        # stylish DNA
        #pymol.cmd.cartoon('oval')
        #pymol.cmd.set('cartoon_oval_length', '1')
        #pymol.cmd.set('cartoon_oval_width', '0.2')
        pymol.cmd.set('cartoon_ring_mode', '3')
        pymol.cmd.set('cartoon_ring_finder', '2')
        pymol.cmd.set('cartoon_ring_width', '0.2')
        
        # turn off all reflections
        #pymol.cmd.set('reflect', '0')
        #pymol.cmd.set('light_count', '1')
        #pymol.cmd.set('ambient', '1')

        # cartoonish ray setting
        pymol.cmd.hide('lines')
        pymol.cmd.set('ray_trace_mode', '3')
        pymol.cmd.bg_color('white')
        pymol.cmd.set('antialias', '5')
        ##pymol.cmd.set('ray_trace_fog', '0')
        pymol.cmd.set('ray_shadows', '0')

        
        #pymol.cmd.space('pymol')
        #pymol.cmd.cartoon('oval', string selection)        
        return


    def make_sphere (self, chain_resi):

        selections = []
        for chain in chain_resi:
            select = []
            for resi in chain_resi[chain].keys():
                if resi < 0:
                    select.append('\\' + str(resi))
                else:
                    select.append(str(resi))
            select = ','.join(select)
            selections.append("(" + "chain " + chain + " and " + "resi " + select + ")")
        selections = " or ".join(selections)

        pymol.cmd.show('sphere', selections)
        pymol.cmd.set("sphere_scale", '0.8')
        return

    def coloring (self, chain_resi, color):
        for chain in chain_resi:
            for resi in chain_resi[chain].keys():
                if resi < 0:
                    select = '\\' + str(resi)
                else:
                    select = str(resi)
                pymol.cmd.color(color, "(" + "chain " + chain + " and " + "resi " + select + ")")
        return
        

    def spectrum(self,
                 chain_resi_value,
                 color_list=[],
                 min=None,
                 max=None):

        selections = []
        for chain in chain_resi_value:
            select = []
            for resi in chain_resi_value[chain].keys():
                if resi < 0:
                    select.append('\\' + str(resi))
                else:
                    select.append(str(resi))
            select = ','.join(select)
            selections.append("(" + "chain " + chain + " and " + "resi " + select + ")")
        selections = " or ".join(selections)

        myspace = {"chain_resi_value":chain_resi_value}
        pymol.cmd.alter(selections, 'b=chain_resi_value[chain][resv]', space=myspace)
        
        if len(color_list) > 0:
            palette = " ".join(color_list)
        else:
            palette = "rainbow"
            
        pymol.cmd.spectrum('b', palette, selections, min, max)


    def save_session (self, fname):
        pymol.cmd.save("%s.pse" % (fname))
        return

    def save_pdb (self, note=None):
        if note != None:
            note = '_' + note
        pymol.cmd.save("%s%s.pdb" % (self.code, note)) 

    def clear_up (self):
        pymol.cmd.delete('all')
        return

    def done (self):
        pymol.cmd.quit()
        return

    def acetylate (self, chain_resi):
        selections = self.selection_string(chain_resi)
        pytms.acetylate(selection)
        return

    def methylate (self, chain_resi):
        selections = self.selection_string(chain_resi)
        pytms.methylate(selection)
        return

    def phosphorylate (self, chain_resi):
        selections = self.selection_string(chain_resi)
        pytms.phosphorylate(selection)
        return

    
# add PTM modifications
#NCP = Molecules("1kx5")
#pymol.cmd.remove('resn hoh')
#pymol.cmd.remove('resn mn')
#pymol.cmd.remove('resn cl')
#pymol.cmd.hide('all')
#pytms.acetylate(selection="chain F and resi 77")
#pymol.cmd.extend(pytms.acetylate)
#pymol.cmd.acetylate("chain F and resi 77")
#pymol.cmd.show('sphere', "chain F and resi 77")
#NCP.save_session("PTM_test")
#NCP.save_pdb("ac")

# load structure
#NCP = Molecules("1kx5")
#histone_chains = {'H2A':['C', 'G'], 'H2B':['D', 'H'], 'H3':['A', 'E'], 'H4':['B','F']}

"""
# load energy profile
path = "/home/spark159/Projects/slide-seq/"
fname = "mmlib_bubble_5_1rep_energy_wrt_601"
with open(path+fname + ".pickle", "rb") as f:
    size_dyad_shl_values = pickle.load(f,encoding='latin1')

# load structure
NCPandChd1 = Molecules("5O9G")
NCPandChd1.print_seq()
#NCPandChd1.spectrum({"A":{38:1, 39:2}}, ["blue", "white", "red"], min=1, max=2)

# assign value on Widom 601 sequence
Widom_seq = "ACAGGATGTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCTGCAGCACCGGGATTCTCCAG"
top_strand = NCPandChd1.chain_seq['I']
bott_strand = NCPandChd1.chain_seq['J']

size = 3
dyad = 92

chain_resi_value = {}
for shl in size_dyad_shl_values[size][dyad]:
    value = np.mean(size_dyad_shl_values[size][dyad][shl])

    if 'I' not in chain_resi_value:
        chain_resi_value['I'] = {}
    chain_resi_value['I'][shl] = value

    if 'J' not in chain_resi_value:
        chain_resi_value['J'] = {}
    chain_resi_value['J'][-shl] = value

# make spectrum plot
pymol.cmd.log_open("test.pml")
NCPandChd1.spectrum(chain_resi_value, color_list=['white', 'red'], min=-0.5, max=3)
#NCPandChd1.spectrum(chain_resi_value, min=-0.5, max=3)
pymol.cmd.log_close()

#pymol.cmd.save("test.pse")
pymol.cmd.quit()
#pymol.finish_launching()
#pymol.cmd.hide("all")
#pymol.cmd.show("")



# Load Structures
#pymol.cmd.fetch("5O9G")

# Display the information of structures
#for x in pymol.cmd.get_names():
#    print ('Protein', x)
#    for ch in pymol.cmd.get_chains(x):
#        print (x, "has chain ", ch)





#pymol.cmd.hide("all")
#pymol.cmd.show("")
#pymol.cmd.disable("all")
#pymol.cmd.enable(sname)
#pymol.cmd.png("my_image.png")

# Get out!
#pymol.cmd.quit()
"""
