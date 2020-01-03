import sys
import math
import random
import copy
import heapq
import numpy as np
import scipy.special as special
from scipy.integrate import quad
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import load_file

"""
# parameters
rmin, rmax = 0, 50 # minimum and maximum reaction distance (nm)
p_len = 50 # DNA persistence length (nm)
scale = 0.413 # total count ratio of test sample to control (UV data, sp9)
#scale = 0.154 # total count ratio of test sample to control (UV data, sp10)
occ_limit = 4 # occupancy number limit

# Random walk chain simulation
def random_walk (N, b):
    assert N >= 1 and b > 0
    R = [0.0, 0.0, float(b)]
    for i in range(N-1):
        u, v = random.uniform(-1, 1), random.random()
        x = b*np.sqrt(1-u**2)*np.cos(2*np.pi*v)
        y = b*np.sqrt(1-u**2)*np.sin(2*np.pi*v)
        z = b*u
        R[0] +=x; R[1] +=y; R[2] +=z    
    endToend = 0.0
    for k in range(3):
        endToend += R[k]**2
    return np.sqrt(endToend)

# Random walk chain radial function
def RW_radial (r, L, p_len):
    b = 2*p_len # kuhn length = 2 * persistence length        
    prob = ((3.0 / (2*np.pi*L*b))**(3.0/2))*np.exp(-3*(r**2)/ float(2*L*b))
    return 4*np.pi*r*r*prob

# Worm like chain radial function (Becker, et al, EPJ, 2010) 
def WLC_radial (r, L, p_len):
    rho = float(r)/L
    lamb = float(p_len)/L

    a, b = 14.054, 0.473
    c = 1.0 - (1 + (0.38*(lamb**(-0.95)))**(-5))**(-1.0/5)
    d = 1.0
    if lamb >= 1.0/8:
        den = 0.177/(lamb-0.111) + 6.4*((lamb-0.111)**(0.783))
        d += - 1.0 / den
    C = [[-3.0/4, 23.0/64, -7.0/64], [-1.0/2, 17.0/16, -9.0/16]]

    JSY = 1.0/((2.0*p_len)**3)
    #JSY = 1.0
    if lamb >= 1.0/8:
        JSY *= 896.32*(lamb**5)*np.exp(-14.054*lamb + 0.246/lamb)
    else:
        JSY *= ((3.0*lamb/np.pi)**(3.0/2))*(1.0 - 5.0*lamb/4.0 - 79.0*lamb*lamb/160.0)
        
    term1 = 0.0
    for i in [-1, 0]:
        for j in [1, 2, 3]:
            term1 += C[i+1][j-1]*(lamb**i)*(rho**(2*j))
    term2 = (-d*lamb*a*(1+b)*rho)/(1.0-b*b*rho*rho)
            
    prob = JSY*(((1.0-c*rho*rho)/(1.0-rho*rho))**(5.0/2))*np.exp(term1/(1.0-rho*rho))
    prob *= np.exp(term2*b*rho)
    prob *= special.iv(0, term2)

    return 4*np.pi*r*r*prob

def reaction_prob (score1, score2, metric, scale=scale):
    prob1 = 1.0 - scale*np.exp(-score1)
    if prob1 < 0:
        prob1 = 0
    if prob1 > 1:
        prob1 = 1
    #assert prob1 >=0 and prob1 <= 1
    prob2 = 1.0 - scale*np.exp(-score2)
    if prob2 < 0:
        prob2 = 0
    if prob2 > 1:
        prob2 = 1
    #assert prob2 >=0 and prob2 <= 1
    if metric == "product":
        prob = prob1*prob2
    elif metric == "GM":
        prob = np.sqrt(prob1*prob2)
    return prob

# load annotation file
path = "/home/spark159/../../media/spark159/sw/dataforcondense//hg19_chr1_171_everything_anot.cn"
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path)
ID_score = name_ID_value['work/condense_seq/sp9_hg19_chr1']
temp = {}
NCPst = 0
NCPnum = 100
for i in range(NCPst, NCPst+NCPnum):
    temp[i] = ID_score[i]
ID_score = temp
ID_AT = name_ID_value['ATcontent']

IDs = sorted(ID_score.keys())
ATproduct_matrix = np.zeros((len(IDs), len(IDs)))
for i in range(len(IDs)-1):
    for j in range(i+1, len(IDs)):
        ID1, ID2 = IDs[i], IDs[j]
        ATproduct = ID_AT[ID1]*ID_AT[ID2]
        ATproduct_matrix[i][j] = ATproduct
        ATproduct_matrix[j][i] = ATproduct

fig = plt.figure()
plt.imshow(ATproduct_matrix, cmap='Reds')
plt.colorbar()
plt.title("AT content x AT content")
#plt.show()
plt.close()

Prob_matrix = np.zeros((len(IDs), len(IDs)))
for i in range(len(IDs)-1):
    for j in range(i+1, len(IDs)):
        ID1, ID2 = IDs[i], IDs[j]
        Prob = reaction_prob (ID_score[ID1], ID_score[ID2], "product")
        Prob_matrix[i][j] = Prob
        Prob_matrix[j][i] = Prob

fig = plt.figure()
plt.imshow(Prob_matrix, cmap='Reds')
plt.colorbar()
plt.title("Reaction probability")
#plt.show()
plt.close()

Coll_matrix = np.zeros((len(IDs), len(IDs)))
for i in range(len(IDs)-1):
    for j in range(i+1, len(IDs)):
        ID1, ID2 = IDs[i], IDs[j]
        dist = abs(ID_pos[ID2] - ID_pos[ID1])
        p = quad(WLC_radial, min(rmin, dist), min(rmax, dist), args=(dist, p_len))[0]
        Coll_matrix[j][i] = p
        Coll_matrix[i][j] = p

fig = plt.figure()
plt.imshow(Coll_matrix, cmap='Reds')
plt.colorbar()
plt.title("Collision probability")
#plt.show()
plt.close()

Simple_matrix = np.zeros((len(IDs), len(IDs)))
for i in range(len(IDs)-1):
    for j in range(i+1, len(IDs)):
        #value = (Coll_matrix[i][j]*Prob_matrix[i][j]) / (Coll_matrix[i][j]*Prob_matrix[i][j]+1-Prob_matrix[i][j])
        value = (Coll_matrix[i][j]*Prob_matrix[i][j])
        Simple_matrix[i][j] = value
        Simple_matrix[j][i] = value

fig = plt.figure()
plt.imshow(Coll_matrix, cmap='Reds')
plt.colorbar()
plt.title("Collision x Reaction")
#plt.show()
plt.close()
"""

#fname = 'output2_bor'
fname = 'bor_N10_C9_bor'
mean_fname = fname + '_mean.txt'
ctnum_fname = fname + '_ctnum.txt'
trace_fname = fname + '_trace.txt'

def read_mean (mean_fname):
    mean_matrix = []
    for line in open(mean_fname):
        line = line.strip()
        if not line:
            continue
        if line.startswith('IDs'):
            _, IDs = line.split(':')
            IDs = [int(ID) for ID in IDs.split(',')]
            continue
        row = [float(value) for value in line.split(',')]
        assert len(row) == len(IDs)
        mean_matrix.append(row)
    assert len(mean_matrix) == len(IDs)
    return IDs, mean_matrix
#IDs, mean_matrix = read_mean (mean_fname)
#mean_matrix = np.zeros((1000,1000))
#for i in range(100):
#    mean_fname = "/home/spark159/scripts/condense-seq/tempdir/sp10_" + str(i+1) + "_rca_mean.txt"
#    #mean_fname = "check" + str(i+1) + '_rca_mean.txt'
#    _, temp_matrix = read_mean (mean_fname)
#    mean_matrix += temp_matrix
#mean_matrix = mean_matrix/100

def read_ctnum (ctnum_fname):
    ctnum_list = []
    for line in open(ctnum_fname):
        line = line.strip()
        if not line:
            continue
        if line.startswith('CycleNum'):
            cycle_num = int(line.split(':')[1])
            continue
        if line.startswith("IDs"):
            continue
        temp = [float(value) for value in line.split(',')]
        ctnum_list += temp
    #assert len(ctnum_list) == cycle_num
    return ctnum_list
#ctnum_list = read_ctnum (ctnum_fname)
#energy_list = read_ctnum (energy_fname)

def make_mean_matrix (trace_fname):
    #ctnum = 0
    #ctnum_list = []
    pre_cycle = 1
    for line in open(trace_fname):
        line = line.strip()
        if not line:
            continue
        if line.startswith('IDs'):
            IDs = [int(value) for value in line[4:].split(',')]
            ID_index = {}
            for i in range(len(IDs)):
                ID = IDs[i]
                ID_index[ID] = i
            mean_matrix = np.zeros((len(IDs), len(IDs)))
            D = np.zeros((len(IDs), len(IDs)))
            continue
        if line.startswith('CycleNum'):
            #cycle_num = int(line.split(':')[1])
            #cycle_num = 24264135
            #cycle_num = 88147315
            cycle_num = 106520234
            burn_in = int(cycle_num*0.1)
            continue
        if line.startswith('Initial'):
            continue
        if line.startswith('@'):
            cur_cycle = int(line[1:])
            repeat = cur_cycle - max(pre_cycle, burn_in+1)
            pre_cycle = cur_cycle
            if repeat > 0:
                mean_matrix += D*repeat
                #ctnum_list += [ctnum]*repeat
            continue
        type, contact = line.split(':')
        ID1, ID2 = contact.split(',')
        ID1, ID2 = int(ID1), int(ID2)
        idx1, idx2 = ID_index[ID1], ID_index[ID2]
        if type == '+':
            D[idx1][idx2] += 1
            D[idx2][idx1] += 1
            #ctnum +=1
        elif type == '-':
            D[idx1][idx2] -= 1
            D[idx2][idx1] -= 1
            #ctnum -=1
    repeat = cycle_num - max(pre_cycle, burn_in+1) + 1
    mean_matrix += D*repeat
    #ctnum_list += [ctnum]*repeat
    #assert len(ctnum_list) == cycle_num - burn_in
    mean_matrix = mean_matrix * (1.0/(cycle_num-burn_in))
    return IDs, mean_matrix
IDs, mean_matrix = make_mean_matrix(trace_fname)

def Video (trace_fname, num=None, step=None):
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='Movie Test', artist='Matplotlib', comment='Movie support!')
    writer = FFMpegWriter(fps=15, metadata=metadata)
    fig = plt.figure()
    pre_cycle = 1
    factor = 0
    with writer.saving(fig, "writer_test.mp4", 100):
        for line in open(trace_fname):
            line = line.strip()
            if not line:
                continue
            if line.startswith('IDs'):
                _, IDs = line.split(':')
                IDs = [int(ID) for ID in IDs.split(',')]
                ID_index = {}
                for i in range(len(IDs)):
                    ID = IDs[i]
                    ID_index[ID] = i
                D = np.zeros((len(IDs), len(IDs)))
                continue
            if line.startswith("CycleNum"):
                cycle_num = int(line.split(':')[1])
                #cycle_num = 10000000
                if not num:
                    num = 100
                if not step:
                    step = cycle_num / num
                continue
            if line.startswith('@'):
                cur_cycle = int(line[1:])
                while factor*step < cur_cycle:
                    plt.spy(D)
                    ax = plt.gca()
                    ax.xaxis.set_ticks_position('bottom')
                    plt.title(str(factor*step) + " th step")
                    writer.grab_frame()
                    plt.draw()
                    plt.pause(0.01)
                    factor += 1
                pre_cycle = cur_cycle
                continue
            type, contact = line.split(':')
            ID1, ID2 = contact.split(',')
            ID1, ID2 = int(ID1), int(ID2)
            idx1, idx2 = ID_index[ID1], ID_index[ID2]
            if type == '+':
                D[idx1][idx2] += 1
                D[idx2][idx1] += 1
            elif type == '-':
                D[idx1][idx2] -= 1
                D[idx2][idx1] -= 1
        while factor*step <= cycle_num:
            plt.spy(D)
            ax = plt.gca()
            ax.xaxis.set_ticks_position('bottom')
            plt.title(str(factor*step) + " th step")
            plt.draw()
            plt.pause(0.1)
            factor += 1
        plt.close()
        return
#Video("/home/spark159/scripts/condense-seq/tempdir/sp10_" + str(1) + "_rca_trace.txt")
#Video(trace_fname)
#Video("test10_rc_sim.txt")

fig = plt.figure()
plt.imshow(mean_matrix, cmap='Reds')
plt.colorbar()
plt.title("Mean Contact Probability")
plt.show()
plt.close()

"""
fig = plt.figure()
plt.plot(range(1, len(ctnum_list)+1), ctnum_list)
plt.xlabel("Steps")
plt.ylabel("Contact counts")
plt.show()
plt.close()


fig = plt.figure()
plt.plot(range(1, len(energy_list)+1), energy_list)
plt.xlabel("Steps")
plt.ylabel("Energy")
plt.show()
plt.close()


def ctprob_over_dist (IDs, mean_matrix, bin_size=None):
    dist_probs = {}
    for i in range(len(IDs)-1):
        for j in range(i+1, len(IDs)):
            ID1, ID2 = IDs[i], IDs[j]
            prob = mean_matrix[i][j]
            dist = ID_pos[ID2] - ID_pos[ID1]
            if bin_size:
                factor = dist // bin_size
                dist = factor*bin_size + bin_size*0.5
            if dist not in dist_probs:
                dist_probs[dist] = []
            dist_probs[dist].append(prob)
    return dist_probs
dist_contact_probs = ctprob_over_dist (IDs, mean_matrix, bin_size=100)
dist_collision_probs = ctprob_over_dist (IDs, Coll_matrix, bin_size=100)
dist_simple_probs = ctprob_over_dist (IDs, Simple_matrix, bin_size=100)

X1, Y1, Z1 = [], [], []
for dist in sorted(dist_contact_probs.keys()):
    X1.append(dist)
    Y1.append(np.mean(dist_contact_probs[dist]))
    Z1.append(np.std(dist_contact_probs[dist])/np.sqrt(len(dist_contact_probs[dist])))
    #Z1.append(np.std(dist_contact_probs[dist]))

X2, Y2, Z2 = [], [], []
for dist in sorted(dist_collision_probs.keys()):
    X2.append(dist)
    Y2.append(np.mean(dist_collision_probs[dist]))
    Z2.append(np.std(dist_collision_probs[dist])/np.sqrt(len(dist_collision_probs[dist])))
    #Z2.append(np.std(dist_collision_probs[dist]))

X3, Y3, Z3 = [], [], []
for dist in sorted(dist_simple_probs.keys()):
    X3.append(dist)
    Y3.append(np.mean(dist_simple_probs[dist]))
    Z3.append(np.std(dist_simple_probs[dist])/np.sqrt(len(dist_simple_probs[dist])))
    #Z3.append(np.std(dist_simple_probs[dist]))


    
fig = plt.figure()
plt.plot(X1, Y1, '.--', markersize=7, color='black', label='Contact probability')
plt.errorbar(X1, Y1, yerr=Z1, fmt='.', color='black', ecolor='green', alpha=0.5)
plt.plot(X2, Y2, '.--', markersize=7, color='red', label='Collision probability')
plt.errorbar(X3, Y3, yerr=Z3, fmt='.', color='red', ecolor='green', alpha=0.5)
plt.plot(X3, Y3, '.--', markersize=7, color='blue', label='Collision x Reaction')
plt.errorbar(X3, Y3, yerr=Z3, fmt='.', color='blue', ecolor='green', alpha=0.5)
plt.xscale("log")
#plt.yscale("log")
plt.legend()
plt.xlabel("Distance (bp)")
plt.ylabel("Probability")
plt.show()
plt.close()


            
#def Read_trace(trace_fname, step=100000):
step = 1000
cycle_num = 1
ctnum_list = []
ctnum = 0
factor = -1
#fig = plt.figure()
for line in open(trace_fname):
    line = line.strip()
    if not line:
        continue
    if line.startswith('IDs'):
        _, IDs = line.split(':')
        IDs = [int(ID) for ID in IDs.split(',')]
        ID_index = {}
        for i in range(len(IDs)):
            ID = IDs[i]
            ID_index[ID] = i
        mean_matrix = np.zeros((len(IDs), len(IDs)))
        D = np.zeros((len(IDs), len(IDs)))
        continue
    if line.startswith('@'):
        current_num = int(line[1:])
        repeat = current_num - cycle_num
        #print D
        mean_matrix += D*repeat
        cycle_num = current_num
        ctnum_list += [ctnum]*repeat
        #if count / step > factor:
        #    factor = count / step
        #    plt.spy(D)
        #    plt.draw()
        #    plt.pause(0.1)
        continue
    type, contact = line.split(':')
    ID1, ID2 = contact.split(',')
    ID1, ID2 = int(ID1), int(ID2)
    idx1, idx2 = ID_index[ID1], ID_index[ID2]
    if type == '+':
        D[idx1][idx2] += 1
        D[idx2][idx1] += 1
        ctnum +=1
    elif type == '-':
        D[idx1][idx2] -= 1
        D[idx2][idx1] -= 1
        ctnum -=1
#repeat = count - cycle_num
repeat = 1
mean_matrix += D*repeat
ctnum_list += [ctnum]*repeat
mean_matrix = mean_matrix / cycle_num
plt.close()

#IDs, mean_matrix, ctnum_list = Read_trace ("output_rc_sim.txt")

fig = plt.figure()
plt.imshow(mean_matrix, cmap='Reds')
plt.colorbar()
plt.title("Mean Contact Probability")
plt.show()
plt.close()

fig = plt.figure()
plt.plot(range(1, len(ctnum_list)+1), ctnum_list)
plt.xlabel("Steps")
plt.ylabel("Contact counts")
plt.show()
plt.close()


def count_contacts (Trace):
    count_list = []
    for i in range(len(Trace)):
        count = 0
        ID_contacts = Trace[i]
        for contacts in ID_contacts.values():
            count += len(contacts)
        assert count % 2 == 0
        count_list.append(count/2)
    return count_list
count_list = count_contacts(Trace)

fig = plt.figure()
plt.plot(range(1, len(count_list)+1), count_list)
plt.xlabel("Steps")
plt.ylabel("Contact counts")
plt.show()
plt.close()

def make_contact_matrix (IDs, ID_contacts):
    IDs, N = sorted(IDs), len(IDs)
    
    ID_index = {}
    for i in range(N):
        ID_index[IDs[i]] = i
        
    matrix = np.zeros((N, N))
    for ID1 in ID_contacts:
        for ID2 in ID_contacts[ID1]:
            if ID2 > ID1:
                i, j = ID_index[ID1], ID_index[ID2]
                matrix[i][j] = 1
                matrix[j][i] = 1

    return matrix

def Video (IDs, Trace, step):    
    fig = plt.figure()
    for i in range(0, len(Trace), step):
        matrix = make_contact_matrix(IDs, Trace[i])
        plt.spy(matrix)
        plt.draw()
        plt.pause(0.1)
    plt.close()
    return None
Video(IDs, Trace, 10000)

def mean_contact_matrix (IDs, Trace, burn_in=0):
    mean = np.zeros((len(IDs), len(IDs)))
    for i in range(burn_in, len(Trace)):
        matrix = make_contact_matrix(IDs, Trace[i])
        mean += matrix
    return mean/(len(Trace)-burn_in)
mean_matrix = mean_contact_matrix(IDs, Trace, burn_in=0)


fig = plt.figure()
plt.imshow(mean_matrix, cmap='Reds')
plt.colorbar()
plt.title("Mean Contact Probability")
plt.show()
plt.close()


IDs = sorted(ID_score.keys())
ATproduct_matrix = np.zeros((len(IDs), len(IDs)))
for i in range(len(IDs)-1):
    for j in range(i+1, len(IDs)):
        ID1, ID2 = IDs[i], IDs[j]
        ATproduct = ID_AT[ID1]*ID_AT[ID2]
        ATproduct_matrix[i][j] = ATproduct
        ATproduct_matrix[j][i] = ATproduct

fig = plt.figure()
plt.imshow(ATproduct_matrix)
plt.colorbar()
plt.title("AT content x AT content")
plt.show()
plt.close()

Prob_matrix = np.zeros((len(IDs), len(IDs)))
for i in range(len(IDs)-1):
    for j in range(i+1, len(IDs)):
        ID1, ID2 = IDs[i], IDs[j]
        Prob = reaction_prob (ID_score[ID1], ID_score[ID2], "product")
        Prob_matrix[i][j] = Prob
        Prob_matrix[j][i] = Prob

fig = plt.figure()
plt.imshow(Prob_matrix)
plt.colorbar()
plt.title("Reaction probability")
plt.show()
plt.close()

Contact_matrix = np.zeros((len(IDs), len(IDs)))
for i in range(len(IDs)-1):
    for j in range(i+1, len(IDs)):
        ID1, ID2 = IDs[i], IDs[j]
        dist = abs(ID_pos[ID2] - ID_pos[ID1])
        p = quad(WLC_radial, min(rmin, dist), min(rmax, dist), args=(dist, p_len))[0]
        Contact_matrix[j][i] = p
        Contact_matrix[i][j] = p

fig = plt.figure()
plt.imshow(Contact_matrix)
plt.colorbar()
plt.title("Collision probability")
plt.show()
plt.close()
"""
