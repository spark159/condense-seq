import math
import numpy as np
import matplotlib.pyplot as plt
import load_file


# load file
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file("data/hg19_chr1_171_everything_anot.cn")


def collision_prob(dist):
    
    return

def reaction_prob(score1, score2):
    prob = (1-math.log(-score1))*(1-math.log(-score2))
    return prob


def 
