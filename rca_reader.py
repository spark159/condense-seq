import sys
import math
import random
import copy
import heapq
import numpy as np
import scipy.special as special
from scipy.integrate import quad
from scipy.optimize import root_scalar
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import load_file
import matplotlib as mpl

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
#fname = 'bor_N10_C9_bor'

#fname = 'output_rca'
#mean_fname = fname + '_mean.txt'
#ctnum_fname = fname + '_ctnum.txt'
#trace_fname = fname + '_trace.txt'
#cluster_fname = fname + "_cluster.txt"

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
#cluster_list = read_ctnum (cluster_fname)
#energy_list = read_ctnum (energy_fname)

def read_cluster (cluster_fname):
    ID_size_dist = {}
    First = True
    for line in open(cluster_fname):
        line = line.strip()
        if First:
            First = False
            continue
        if not line:
            continue
        if line.startswith('@'):
            ID = line[1:]
            continue
        ID_size_dist[ID] = [float(value) for value in line.split(',')]
    return ID_size_dist

def read_anot (anot_fname):
    ID_energy = {}
    ID_valency = {}
    First = True
    for line in open(anot_fname):
        line = line.strip()
        if First:
            First = False
            continue
        if not line:
            continue
        cols = line.split()
        ID, energy, valency = cols[0], cols[3], cols[4]
        ID_energy[ID] = float(energy)
        ID_valency[ID] = int(valency)
    return ID_energy, ID_valency

def P_to_energy (ID_P, ID_valency, volume):
    # adjust P
    for ID, P in ID_P.items():
        if P <= 0:
            P = 10**-10
            ID_P[ID] = P
        elif P >= 1:
            P = 1.0 - 10**-10
            ID_P[ID] = P
    # find Pm and Em
    N = len(ID_P)
    fm = np.mean(ID_valency.values())
    Pm = 0.0
    for ID in ID_P:
        P, f = ID_P[ID], ID_valency[ID]
        Pm += np.log(P*f)
        #Pm += np.log(P*(f**2))
    Pm = (np.exp(Pm/float(N)))/float(fm)
    #Pm = (np.exp(Pm/float(N)))/float(fm**2)
    Em = np.log((fm*((fm-2)**2))/(fm-1)) + np.log(N) + np.log(Pm) - np.log(volume)
    # find E
    ID_energy = {}
    for ID in ID_P:
        P, f = ID_P[ID], ID_valency[ID]
        energy = Em + 2*np.log(float(f*P)/(fm*Pm))
        #energy = Em + 2*np.log(float((f**2)*P)/((fm**2)*Pm))
        ID_energy[ID] = energy
    return ID_energy

def P_to_valency (ID_P, ID_energy, volume):
    # adjust P
    for ID, P in ID_P.items():
        if P <= 0:
            P = 10**-10
            ID_P[ID] = P
        elif P >= 1:
            P = 1.0 - 10**-10
            ID_P[ID] = P
    # find Pm and fm
    N = len(ID_P)
    Em = np.mean(ID_energy.values())
    fm = 0.0
    for ID in ID_P:
        P, E = ID_P[ID], ID_energy[ID]
        fm += np.exp(0.5*(E-Em))/P
        #fm += np.exp(E) / P
    fm = np.sqrt(volume*np.exp(Em)*fm/float(N**2))
    #print fm
    #fm = np.sqrt(float(volume*fm) / N)
    #fm = 23.0
    # find f
    ID_valency = {}
    for ID in ID_P:
        P, E = ID_P[ID], ID_energy[ID]
        f = ((volume*np.exp(Em)*(1+(fm-2)*np.exp(0.5*(E-Em))))/float(N*((fm-2)**2)*P))
        #f = float(volume*np.exp(E))/(fm*P)
        #f = volume*np.exp(Em)*np.exp(0.5*(E-Em)) / float(N*P*fm)
        ID_valency[ID] = f
    return ID_valency

def residuals (ID_trues, ID_preds):
    output = []
    for ID in ID_trues:
        x, y = ID_trues[ID], ID_preds[ID]
        #total += float(abs(x-y))/math.sqrt(2)
        output.append(abs(x-y))
    return output

# Stockmayer cluster size distribution
def FS_size_dist (N, occ, energy, volume, z=None):
    # log of factorial
    def log_factorial (n):
        total = 0.0
        for i in range(1, n+1):
            total += np.log(i)
        return total
    
    # number of k-fold cluster
    def cluster_count (k,z):
        count = np.log(volume) + k*np.log(occ) + log_factorial(occ*k-k) - log_factorial(occ*k-2*k+2) + k*np.log(z) - energy*(k-1) - log_factorial(k)
        count = np.exp(count)
        return count
    
    # total particle number
    def total_num (z):
        total = 0.0
        for k in range(1, N+1):
            total += k*cluster_count(k, z)
        return total

    N_limit = volume * np.exp(energy) * float(occ-1) / occ*((occ-2)**2) # number limit
    if N > N_limit:
        print "number limit exceed:" + str(N_limit)
        N_s = N_limit - 1
        N_gel = N - N_s
    else:
        N_s = N
        N_gel = 0

    # function has to be solved
    def func (z):
        return total_num(z) - N_s

    if z == None:
        z = root_scalar(func, bracket=[10**-10, float(N)/volume]).root

    size_counts = [] # cluster size distribution    
    total = 0.0 # total particle number
    for k in range(1, N+1):
        count = cluster_count(k, z)
        size_counts.append(count)
        total += k*count
    #print z
    print total
    return size_counts

def number_count (cluster_list):
    total = 0.0
    for i in range(len(cluster_list)):
        total += (i+1)*cluster_list[i]
    return total

# pure case
N = 100
volume = 10**7

energy_list_list = []
P_list_list = []
pred_list_list = []
S_list_list = []
gel_point_list = []
for v in range(10):
    f = 3 + 10*v
    gel_point = math.log(N) + math.log(f) + 2*math.log(f-2) - math.log(volume) - math.log(f-1)
    gel_point_list.append(gel_point)
    energy_list = []
    P_list = []
    pred_list = []
    S_list = []
    for i in range(10):
        E = - 2*i
        fname = 'homo_' + str(N) + '_' + str(f) + '_' + str(E)
        cluster_list = read_ctnum (fname + "_rca_cluster.txt")
        P = number_count(cluster_list[:10])/float(N)
        Nc = volume*(f-1)*np.exp(E) / float(f*((f-2)**2))
        S = sum([cluster_list[k]*float(k+1) for k in range(len(cluster_list))])
        predict = 1.7*Nc/float(N)
        energy_list.append(E)
        P_list.append(P)
        pred_list.append(predict)
        S_list.append(S)
    energy_list_list.append(energy_list)
    P_list_list.append(P_list)
    pred_list_list.append(pred_list)
    S_list_list.append(S_list)

color_list = np.linspace(0.01, 1, num=len(energy_list_list))
cmap = mpl.cm.get_cmap("jet")
fig = plt.figure()
for i in range(len(energy_list_list)):
    f = 3 + 10*i
    #if f not in [23, 13]:
    #    continue
    gel_point = gel_point_list[i]
    energy_list = energy_list_list[i]
    P_list = P_list_list[i]
    pred_list = pred_list_list[i]
    S_list = S_list_list[i]
    plt.plot(energy_list, P_list, '.', color=cmap(color_list[i]), label='f=' + str(f))
    plt.plot(energy_list, pred_list, '--', color=cmap(color_list[i]), alpha=0.5)
    #plt.plot(energy_list, S_list, color=cmap(color_list[i]))
    plt.axvline(x=gel_point, color=cmap(color_list[i]), linestyle='-.', alpha=0.25)
plt.xlabel("Interaction energy (kT)")
plt.ylabel("Survival probability")
plt.title("Pure solution (F-S theory vs simulation)")
plt.yscale('log')
#plt.ylim([-0.02, 0.5])
#plt.ylim([-5, 60])
#plt.xlim([-18.5, -5])
plt.legend()
plt.show()
plt.close()



# almost pure case
N = 100
volume = 10**7
f = 23
E = -10
gel_point = math.log(N) + math.log(f) + 2*math.log(f-2) - math.log(volume) - math.log(f-1)
alpha = 1.0/float(f-1)
Nc = volume*(f-1)*np.exp(E) / float(f*((f-2)**2))
#Nc = volume*np.exp(E) / float((f)**2)
Pm = 1.7*Nc/float(N)
Sm = Pm + (N-Nc)*(1.0-Pm)

energy_list_list = []
P_list_list = []
pred_list_list = []
S_list_list = []
for v in range(10):
    f_single = 3 + 10*v
    energy_list = []
    P_list = []
    pred_list = []
    S_list = []
    for i in range(10):
        E_single = - 2*i
        fname = 'single_' + str(N) + '_' + str(f) + '_' + str(E) + ':' + str(f_single) + '_' + str(E_single)
        ID_size_dist = read_cluster (fname + "_rca_cluster.txt")
        ID_energy, ID_valency = read_anot(fname + "_anot.cn")
        dE = E_single - E
        Pg = number_count(ID_size_dist['total'][:10])/float(N)
        P = sum(ID_size_dist['0'][:10])
        S = sum([ID_size_dist['0'][k]*float(k+1) for k in range(len(ID_size_dist['0']))])
        P = P
        #predict = (float(f)/f_single)*np.exp(0.5*dE)
        predict = Pm*(float(f)/(f_single))*(np.exp(0.5*dE)*(1.0-alpha)+alpha)
        #predict = (float(f)/(f_single))*(np.exp(0.5*dE)*(1.0-alpha)+alpha)
        #predict = ((1.0 + f_single*(float(f-2)**(f-2))/((f-1)**(f-1))*np.exp(-0.5*dE))*((1.0+1.0/(f-2))**f))/(((1.0 + np.exp(-0.5*dE)/(f-2))**f_single)*(1.0 + f*(float(f-2)**(f-2))/((f-1)**(f-1))))
        energy_list.append(E_single)
        P_list.append(P)
        pred_list.append(predict)
        S_list.append(S)
    energy_list_list.append(energy_list)
    P_list_list.append(P_list)
    pred_list_list.append(pred_list)
    S_list_list.append(S_list)

color_list = np.linspace(0.01, 1, num=len(energy_list_list))
cmap = mpl.cm.get_cmap("jet")
fig = plt.figure()
for i in range(len(energy_list_list)):
    f_single = 3 + 10*i
    #if f not in [23, 13]:
    #    continue
    energy_list = energy_list_list[i]
    dE_list = [energy - E for energy in energy_list]
    df = f_single - f 
    P_list = P_list_list[i]
    pred_list = pred_list_list[i]
    S_list = S_list_list[i]
    plt.plot(dE_list, P_list, '.', color=cmap(color_list[i]), label='$\Delta$f=' + str(df))
    plt.plot(dE_list, pred_list, '--', color=cmap(color_list[i]), alpha=0.5)
    #plt.plot(dE_list, S_list, color=cmap(color_list[i]))
plt.axvline(x=gel_point-E, color='k', linestyle='-.', alpha=0.25)
plt.xlabel("$\Delta$Interaction energy (kT)")
plt.ylabel("Survival probability")
plt.title("Almost pure solution (F-S theory vs simulation)")
plt.yscale('log')
#plt.ylim([-0.02, 0.2])
#plt.ylim([-5, 60])
#plt.xlim([-8.5, 5])
plt.legend()
plt.show()
plt.close()

"""
ID_size_dist = read_cluster("check_rca_cluster.txt")
ID_energy, ID_valency = read_anot("hetero_100_23_-10_anot.cn")


N = len(ID_size_dist) - 1
volume = 10**7

f = ID_valency['1']
E = ID_energy['1']
z = ID_size_dist['total'][0]/volume

size_counts = FS_size_dist (N, f, E, volume, z=z)

temp = [0.0]*N
for ID in ID_size_dist:
    size_dist = ID_size_dist[ID]
    for i in range(len(temp)):
        temp[i] += size_dist[i]

fig = plt.figure()
plt.plot(ID_size_dist['0'], label='test', alpha=0.5)
plt.plot(ID_size_dist["1"], 'k--')
#plt.plot(temp)
plt.plot(ID_size_dist["total"], 'k--')
#plt.plot(size_counts)
plt.legend()
plt.show()
plt.close()


dE = ID_energy['0'] - ID_energy['1']
P = ID_size_dist['0'][0]
S = sum([ID_size_dist['0'][i]*float(i+1) for i in range(len(ID_size_dist['0']))])
alpha = 2*(N - sum(ID_size_dist['total'])) / float(f*N)
#alpha = 1.0 - (1.0/(2*(float(N)/volume)*f*np.exp(-E)))*(np.sqrt(1+4*N*f*np.exp(-E)/float(volume))-1)
#pred_P = (1-alpha)**f
#alpha = 1.0/(f-1)
pred_P = (1.0/(1.0 + (alpha/(1.0-alpha))*np.exp(-0.5*dE)))**ID_valency['0']
pred_S = 1 + float(alpha*ID_valency['0']) / ((1-(f-1)*alpha)*((1-alpha)*np.exp(0.5*dE)+alpha))
#pred_S = N - volume*(f-1)*np.exp(E) / float(f*((f-2)**2))

Pg = sum([(i+1)*ID_size_dist['total'][i] for i in range(10)]) / float(N)

pred_E = np.log((f*((f-2)**2))/(f-1)) + np.log(N) + np.log(Pg) - np.log(volume)
pred_f = np.sqrt(volume*np.exp(E)/(N*Pg))

Ps = sum(ID_size_dist['0'][:10])
Ps_mean = sum(ID_size_dist['1'][:10])
pred_Pm = volume*(f-1)*np.exp(E) / float(N*f*((f-2)**2))
pred_Ps = (f*pred_Pm*np.exp(0.5*dE) / float(ID_valency['0']))
#pred_Ps = (ID_valency['0']*pred_Pm*np.exp(-0.5*dE) / float(f)) #* ((f-1)**f/float((f-2+np.exp(-0.5*dE))**ID_valency['0']))
#pred_Ps = ((f-1)*pred_Pm*np.exp(0.5*dE) / float(ID_valency['0']))
pred_Ps2 = (f**2)*pred_Pm*np.exp(0.5*dE) / (float(ID_valency['0'])**2)

leftHS = np.exp(0.5*dE)
#rightHS = (float(f-1)/(f-2))*(float(ID_valency['0']*Ps)/float(f*pred_Pm) - 1.0/float(f-1))
rightHS = (ID_valency['0']/float(f))*(Ps/float(pred_Pm))
#rightHS = ((ID_valency['0']/float(f))**(0.5))*(Ps/float(pred_Pm))

print leftHS
print rightHS

#print abs(leftHS-rightHS) / min(leftHS, rightHS)

alpha = 1.0/float(f-1)
#alpha = 2*(N - sum(ID_size_dist['total'])) / float(f*N)
leftHS2 = np.exp(0.5*dE)
rightHS2 = (ID_valency['0']*(1-pred_Pm))/float(f*(1-Ps))

#print abs(leftHS2-rightHS2) / min(leftHS2, rightHS2)

#print 1-Ps
#print S / (N - volume*(f-1)*np.exp(E) / float(f*((f-2)**2)))

print 

leftHS = Ps / pred_Pm
rightHS = 2.0 / ((ID_valency['0']/f)*np.exp(-0.5*dE) + 1)
rightHS_p = (f/ID_valency['0'])*np.exp(0.5*dE)

print leftHS
print rightHS
print rightHS_p

leftHS2 = (1.0 - Ps) / (1.0 - pred_Pm)
rightHS2 = ((ID_valency['0']/f)*np.exp(-0.5*dE) + 1) / 2.0

print
print leftHS2
print rightHS2
"""


"""
fig = plt.figure()
for ID in ID_size_dist:
    if ID == "total":
        continue
    plt.plot(ID_size_dist[ID], label=ID, alpha=0.5)
plt.plot(ID_size_dist["total"], 'k--')
plt.legend()
plt.show()
plt.close()

volume = 10**7
ID_P = {}
ID_valency2 = {}
for ID in ID_size_dist:
    if ID == 'total':
        continue
    size_dist = ID_size_dist[ID]
    #P = sum([size_dist[i]/float(i+1) for i in range(len(size_dist))])
    #P = float(P) / sum(size_dist)
    #P = 1.0 / sum([size_dist[i]*float(i+1) for i in range(len(size_dist))])
    P = sum(ID_size_dist[ID][:10]) / float(sum(ID_size_dist[ID]))
    #print sum([size_dist[i]*float(i+1) for i in range(len(size_dist))])
    #print P
    ID_P[ID] = P
    ID_valency2[ID] = 23
ID_pred_energy = P_to_energy(ID_P, ID_valency, volume)
ID_pred_energy2 = P_to_energy(ID_P, ID_valency2, volume)
ID_pred_valency = P_to_valency (ID_P, ID_energy, volume)

print sum(residuals(ID_energy, ID_pred_energy))

IDs = ID_energy.keys()
X = [ID_valency[ID] for ID in IDs]
Y = [ID_pred_valency[ID] for ID in IDs]
Z = [ID_energy[ID] for ID in IDs]

fig = plt.figure()
plt.scatter(X, Y, c=Z, s=3)
#plt.plot(X, Y2, 'k.', alpha=0.5, label='control')
plt.plot([min(X)-2, max(X)+2], [min(X)-2, max(X)+2], 'k--')
plt.xlim([min(X)-2, max(X)+2])
plt.ylim([min(X)-2, max(X)+2])
plt.title("Mean-field theory validation")
plt.xlabel("Energy used for simulation")
plt.ylabel("F-S theory prediction")
plt.colorbar()
#plt.legend()
plt.show()
plt.close()


IDs = ID_energy.keys()
X = [ID_energy[ID] for ID in IDs]
Y = [ID_pred_energy[ID] for ID in IDs]
Y2 = [ID_pred_energy2[ID] for ID in IDs]
Z = [ID_valency[ID] for ID in IDs]
Z2 = [ID_valency2[ID] for ID in IDs]

fig = plt.figure()
plt.scatter(X, Y, c=Z, s=3)
plt.plot(X, Y2, 'k.', alpha=0.5, label='control')
plt.plot([min(X)-2, max(X)+2], [min(X)-2, max(X)+2], 'k--')
plt.xlim([min(X)-2, max(X)+2])
plt.ylim([min(X)-2, max(X)+2])
plt.title("Mean-field theory validation")
plt.xlabel("Energy used for simulation")
plt.ylabel("F-S theory prediction")
plt.colorbar()
#plt.legend()
plt.show()
plt.close()
"""


"""
N = 100
volume = 10**7
f = 23

gel_point = math.log(N) + math.log(f) + 2*math.log(f-2) - math.log(volume) - math.log(f-1)

ID_P = {}
Pg = 0.0
for ID in ID_size_dist:
    if ID == 'total':
        continue
    P = sum(ID_size_dist[ID][:10]) / float(sum(ID_size_dist[ID]))
    print P
    if P >= 1:
        P = 0.9999999999
    elif P <= 0:
        P = 0.0000000001
    Pg += np.log(P)
    ID_P[ID] = P
#Pg = np.mean(ID_P.values())
Pg = np.exp(Pg/float(N))

ID_pred_energy = {}
for ID, P in ID_P.items():    
    pred_energy = math.log((f*((f-2)**2))/(f-1)) + math.log(N) + math.log(Pg) - math.log(volume)
    #pred_energy += 2*(math.log(((P**(1.0/f))*(1-(Pg**(1.0/f))))/((Pg**(1.0/f))*(1-(P**(1.0/f))))))
    #pred_energy = np.log((f*(f-2)**2)/(f-1)) + np.log(N) + np.log(P) - np.log(volume)
    pred_energy += 2*(np.log(P) - np.log(Pg))
    ID_pred_energy[ID] = pred_energy

IDs = ID_energy.keys()
X = [ID_energy[ID] for ID in IDs]
Y = [ID_pred_energy[ID] for ID in IDs]

fig = plt.figure()
plt.plot(X, Y, '.')
plt.plot([min(X)-2, max(X)+2], [min(X)-2, max(X)+2], 'k--')
plt.xlim([min(X)-2, max(X)+2])
plt.ylim([min(X)-2, max(X)+2])
plt.xlabel("Energy used for simulation")
plt.ylabel("F-S theory prediction")
plt.show()
plt.close()
"""
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
#IDs, mean_matrix = make_mean_matrix(trace_fname)

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

# Stockmayer cluster size distribution
def FS_size_dist (N, occ, energy, volume, z=None):
    # log of factorial
    def log_factorial (n):
        total = 0.0
        for i in range(1, n+1):
            total += np.log(i)
        return total
    
    # number of k-fold cluster
    def cluster_count (k,z):
        count = np.log(volume) + k*np.log(occ) + log_factorial(occ*k-k) - log_factorial(occ*k-2*k+2) + k*np.log(z) - energy*(k-1) - log_factorial(k)
        count = np.exp(count)
        return count
    
    # total particle number
    def total_num (z):
        total = 0.0
        for k in range(1, N+1):
            total += k*cluster_count(k, z)
        return total

    N_limit = volume * np.exp(energy) * float(occ-1) / occ*((occ-2)**2) # number limit
    if N > N_limit:
        print "number limit exceed:" + str(N_limit)
        N_s = N_limit - 1
        N_gel = N - N_s
    else:
        N_s = N
        N_gel = 0

    # function has to be solved
    def func (z):
        return total_num(z) - N_s

    if z == None:
        z = root_scalar(func, bracket=[10**-10, float(N)/volume]).root

    size_counts = [] # cluster size distribution    
    total = 0.0 # total particle number
    for k in range(1, N+1):
        count = cluster_count(k, z)
        size_counts.append(count)
        total += k*count
    #print z
    print total
    return size_counts

def mean_cluster_size (size_counts):
    total_num = 0.0
    for i in range(len(size_counts)):
        counts = size_counts[i]
        total_num += (i+1)*counts
    return total_num/sum(size_counts)


"""
N = 100
volume = 10**7

f = 23
energy_list = []
P_list = []
pre_energy_list = []
mean_size_list = []
cluster_list_list = []
gel_point = math.log(N) + math.log(f) + 2*math.log(f-2) - math.log(volume) - math.log(f-1)
for i in range(10):
    energy = - 2*i
    fname = 'homo_' + str(N) + '_' + str(f) + '_' + str(energy) + '_rca'
    #mean_fname = fname + '_mean.txt'
    #ctnum_fname = fname + '_ctnum.txt'
    #trace_fname = fname + '_trace.txt'
    cluster_fname = fname + "_cluster.txt"
    cluster_list = read_ctnum (cluster_fname)
    P = number_count(cluster_list[:10])/float(N)
    P_list.append(P)
    cluster_num = float(sum(cluster_list))
    cluster_list = [value/cluster_num for value in cluster_list]
    pre_energy = np.log((f*(f-2)**2)/(f-1)) + np.log(N) + np.log(P) - np.log(volume)
    pre_energy_list.append(pre_energy)
    cluster_list_list.append(cluster_list)
    mean_size = mean_cluster_size(cluster_list)
    energy_list.append(energy)
    mean_size_list.append(mean_size/float(100))

fig = plt.figure()
plt.plot(energy_list, mean_size_list, '.-', color='k', label='simulation (f=' +str(f)+')'  )
plt.axvline(x=gel_point, color='r', linestyle='--', label='FS theory (f='+str(f)+')')
plt.xlabel("Interaction energy (kT)")
plt.ylabel("Average cluster size")
plt.title("Flory-Stockmayer simulations, " + "N=" + str(N))
plt.legend()
#plt.savefig("FS_simulation_free.png")
plt.show()
plt.close()

fig = plt.figure()
min_energy, max_energy = min(energy_list), max(energy_list)
plt.plot(energy_list, pre_energy_list, '.')
plt.plot([min_energy, max_energy], [min_energy, max_energy], '--')
plt.axvline(x=gel_point, color='red', linestyle='--', alpha=0.5)
plt.xlabel("Energy used for simulation")
plt.ylabel("F-S theory prediction")
plt.show()
plt.close()

color_list = np.linspace(0.01, 1, num=len(cluster_list_list))
cmap = mpl.cm.get_cmap("jet")
fig = plt.figure()
for i in range(10):
    energy = - 2*i
    cluster_list = cluster_list_list[i]
    plt.plot(range(1, len(cluster_list)+1), cluster_list, '.-', color=cmap(color_list[i]), label='E=' + str(energy) + 'kT')
plt.xscale("log")
plt.xlabel("Cluster size")
plt.ylabel("Fraction")
plt.title("Cluster size distribution, "+ "N=" + str(N))
plt.legend()
#plt.show()
plt.close()


energy_list_list = []
mean_size_list_list = []
gel_point_list = []
for v in range(10):
    f = 3 + 10*v
    energy_list = []
    P_list = []
    pre_energy_list = []
    mean_size_list = []
    gel_point = math.log(N) + math.log(f) + 2*math.log(f-2) - math.log(volume) - math.log(f-1)
    #print gel_point
    for i in range(10):
        energy = - 2*i
        fname = 'homo_' + str(N) + '_' + str(f) + '_' + str(energy) + '_rca'
        #mean_fname = fname + '_mean.txt'
        #ctnum_fname = fname + '_ctnum.txt'
        #trace_fname = fname + '_trace.txt'
        cluster_fname = fname + "_cluster.txt"
        cluster_list = read_ctnum (cluster_fname)
        P = number_count(cluster_list[:10])/float(N)
        P_list.append(P)
        cluster_num = float(sum(cluster_list))
        cluster_list = [value/cluster_num for value in cluster_list]
        pre_energy = np.log((f*(f-2)**2)/(f-1)) + np.log(N) + np.log(P) - np.log(volume)
        pre_energy_list.append(pre_energy)
        mean_size = mean_cluster_size(cluster_list)
        energy_list.append(energy)
        mean_size_list.append(mean_size/float(100))
        
    energy_list_list.append(energy_list)
    mean_size_list_list.append(mean_size_list)
    gel_point_list.append(gel_point)

    fig = plt.figure()
    min_energy, max_energy = min(energy_list), max(energy_list)
    plt.plot(energy_list, pre_energy_list, '.')
    plt.plot([min_energy, max_energy], [min_energy, max_energy], 'k--')
    plt.axvline(x=gel_point, color='red', linestyle='--', alpha=0.5)
    plt.xlabel("Energy used for simulation")
    plt.ylabel("F-S theory prediction")
    #plt.show()
    plt.close()

color_list = np.linspace(0.01, 1, num=len(energy_list_list))
cmap = mpl.cm.get_cmap("jet")
fig = plt.figure()
for i in range(len(energy_list_list)):
    f = 3 + 10*i
    energy_list = energy_list_list[i]
    mean_size_list = mean_size_list_list[i]
    gel_point = gel_point_list[i]
    plt.plot(energy_list, mean_size_list, '.-', color=cmap(color_list[i]), label='f=' + str(f))
    plt.axvline(x=gel_point, color=cmap(color_list[i]), linestyle='--', alpha=0.5)
    #plt.plot(energy_list, mean_size_list, '.-', color=cmap(color_list[i]), label='simulation (f=' +str(f)+')'  )
    #plt.axvline(x=gel_point, color=cmap(color_list[i]), linestyle='--', label='FS theory (f='+str(f)+')')
plt.xlabel("Interaction energy (kT)")
plt.ylabel("Average cluster size")
plt.title("Flory-Stockmayer simulations, " + "N=" + str(N))
plt.legend()
plt.savefig("FS_simulation_free_many.png")
#plt.show()
plt.close()



volume = 10**6
#FS = FS_size_dist(10, 30, -10, volume, z=None)
    
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

fig = plt.figure()
plt.plot(range(1, len(cluster_list)+1), cluster_list, alpha=0.5, label='simulation')
plt.plot(range(1, len(cluster_list)+1), cluster_list, '.')
#plt.plot(range(1, len(cluster_list)+1), FS[:len(cluster_list)], '--', label='FS theory')
plt.legend()
plt.xlabel("Cluster size")
plt.ylabel("Cluster counts")
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
