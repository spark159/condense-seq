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

def read_trace (trace_fname):
    def cluster_update (ID1, ID2, type):
        if type == '+':
            cID1, cID2 = ID_cID[ID1], ID_cID[ID2]
            if cID1 != cID2:
                new_cluster = cID_IDs[cID1] | cID_IDs[cID2]
                new_cID = str(len(new_cluster)) + '-' + str(min(new_cluster))
                del cID_IDs[cID1]; del cID_IDs[cID2]
                cID_IDs[new_cID] = new_cluster
                for ID in new_cluster:
                    ID_cID[ID] = new_cID
            ID_neighbors[ID1].add(ID2)
            ID_neighbors[ID2].add(ID1)
        else:
            ID_neighbors[ID1].remove(ID2)
            ID_neighbors[ID2].remove(ID1)
            cluster1 = BFS(ID_neighbors, ID1)
            if ID2 not in cluster1:
                cluster2 = BFS(ID_neighbors, ID2)
            else:
                cluster2 = cluster1
            if cluster1 != cluster2:
                new_cID1 = str(len(cluster1)) + '-' + str(min(cluster1))
                new_cID2 = str(len(cluster2)) + '-' + str(min(cluster2))
                cID1, cID2 = ID_cID[ID1], ID_cID[ID2]
                assert cID1 == cID2
                del cID_IDs[cID1]
                cID_IDs[new_cID1] = cluster1
                cID_IDs[new_cID2] = cluster2
                for ID in cluster1:
                    ID_cID[ID] = new_cID1
                for ID in cluster2:
                    ID_cID[ID] = new_cID2

    # reading trace file 
    pre_cycle = 1
    for line in open(trace_fname + '_rgs_trace.txt', 'r'):
        line = line.strip()
        if not line:
            continue
        if line.startswith('IDs'):
            IDs = line[4:].split(',')
            
            # mean contact matrix
            mean_matrix = np.zeros((len(IDs), len(IDs)))
            D = np.zeros((len(IDs), len(IDs)))

            # total contact number
            ctnum = 0
            ctnum_list = []

            # cluster information
            ID_neighbors = {ID:set([]) for ID in IDs} 
            cID_IDs = {'-'.join([str(1), str(IDs[i])]):set([IDs[i]]) for i in range(len(IDs))} # cluster ID to particle IDs
            ID_cID = {IDs[i]:'-'.join([str(1), str(IDs[i])]) for i in range(len(IDs))} # particle ID to cluster ID
            mean_size_num = [0]*len(IDs) # cluster size distribution
            ID_mean_size_num = {ID:[0]*len(IDs) for ID in IDs} # cluster size distribution by particle
            continue
        if line.startswith('CycleNum'):
            continue
        if line.startswith('@'):
            cur_cycle = int(line[1:])
            repeat = cur_cycle - max(pre_cycle, burn_in+1)
            pre_cycle = cur_cycle
            if repeat > 0:
                mean_matrix += D*repeat
                ctnum_list += [ctnum]*repeat
                for cluster in cID_IDs.values():
                    size = len(cluster)
                    mean_size_num[size-1] += repeat
                for ID, cID in ID_cID.items():
                    size = len(cID_IDs[cID])
                    ID_mean_size_num[ID][size-1] += repeat
                #for k in range(len(size_num)):
                #    mean_size_num[k] += repeat*size_num[k]
            continue
        type, contact = line.split(':')
        ID1, ID2 = contact.split(',')
        ID1, ID2 = int(ID1), int(ID2)
        idx1, idx2 = ID_index[ID1], ID_index[ID2]
        if type == '+':
            D[idx1][idx2] += 1
            D[idx2][idx1] += 1
            #cluster_update(ID_neighbors, size_num, ID1, ID2, type)
            cluster_update(ID1, ID2, type)
            ctnum +=1
        elif type == '-':
            D[idx1][idx2] -= 1
            D[idx2][idx1] -= 1
            #cluster_update(ID_neighbors, size_inum, ID1, ID2, type)
            cluster_update(ID1, ID2, type)
            ctnum -=1
    repeat = cycle_num - max(pre_cycle, burn_in+1) + 1
    mean_matrix += D*repeat
    ctnum_list += [ctnum]*repeat
    for cluster in cID_IDs.values():
        size = len(cluster)
        mean_size_num[size-1] += repeat
    for ID, cID in ID_cID.items():
        size = len(cID_IDs[cID])
        ID_mean_size_num[ID][size-1] += repeat
    #for k in range(len(size_num)):
    #    mean_size_num[k] += repeat*size_num[k]
    assert len(ctnum_list) == cycle_num - burn_in
    mean_matrix = mean_matrix * (1.0/(cycle_num-burn_in))
    for k in range(len(mean_size_num)):
        mean_size_num[k] = mean_size_num[k] * (1.0/(cycle_num-burn_in))
    for ID in ID_cID.keys():
        for k in range(len(ID_mean_size_num[ID])):
            ID_mean_size_num[ID][k] = ID_mean_size_num[ID][k] * (1.0/(cycle_num-burn_in))

    # writing the results
    f = open(trace_fname + '_rgs_ctprob.txt', 'w')
    print >> f, "IDs" + ":" + ",".join(str(ID) for ID in IDs)

    m, n = mean_matrix.shape
    for i in range(m):
        print >> f, ",".join(str(value) for value in mean_matrix[i])

    f.close()

    f = open(trace_fname + '_rgs_ctnum.txt', 'w')
    print >> f, "CycleNum:" + str(cycle_num)
    print >> f, ",".join(str(value) for value in ctnum_list)

    f.close()

    f = open(trace_fname + '_rgs_cluster.txt', 'w')
    print >> f, "IDs" + ":" + ",".join(str(ID) for ID in IDs)
    print >> f, "@total"
    print >> f, ",".join(str(value) for value in mean_size_num)
    for ID in IDs:
        print >> f, "@" + str(ID)
        print >> f, ",".join(str(value) for value in ID_mean_size_num[ID])

    f.close()

    
    print >> sys.stderr, "Done"
    return

    


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

def read_ctprob (ctprob_fname):
    for line in open(ctprob_fname):
        line = line.strip()
        if not line:
            continue
        if line.startswith('IDs'):
            _, IDs = line.split(':')
            IDs = [int(ID) for ID in IDs.split(',')]
            matrix = np.zeros((len(IDs), len(IDs)))
            matrix[:] = np.nan
            continue
        edge, value = line.split(':')
        value = float(value)
        node1, node2 = edge.split('-')
        node1, node2 = int(node1), int(node2)
        matrix[node1][node2] = value
        matrix[node2][node1] = value
    return IDs, matrix

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

def number_count (cluster_list):
    total = 0.0
    for i in range(len(cluster_list)):
        total += (i+1)*cluster_list[i]
    return total

# Stockmayer cluster size distribution
def FS_size_dist (N, occ, energy, volume, z=None):
    # log of factorial
    def log_factorial (n):
        total = 0.0
        for i in range(1, int(n+1)):
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

    """
    N_limit = volume * np.exp(energy) * float(occ-1) / occ*((occ-2)**2) # number limit
    if N > N_limit:
        print "number limit exceed:" + str(N_limit)
        N_s = N_limit - 1
        N_gel = N - N_s
    else:
        N_s = N
        N_gel = 0
    """
    
    # function has to be solved
    def func (z):
        #return total_num(z) - N_s
        return total_num(z) - N

    if z == None:
        z = root_scalar(func, bracket=[10**-12, float(N)/volume]).root

    size_counts = [] # cluster size distribution    
    total = 0.0 # total particle number
    for k in range(1, N+1):
        count = cluster_count(k, z)
        size_counts.append(count)
        total += k*count
    #print size_counts
    #print z
    #print total
    return size_counts

# FS theory for almost pure case
def FS_size_dist_almost (N, occ, energy, volume, f_single, E_single, z=None):
    # log of factorial
    def log_factorial (n):
        total = 0.0
        for i in range(1, int(n+1)):
            total += np.log(i)
        return total
    
    # number of k-fold mean-cluster
    def cluster_count (k,z):
        count = np.log(volume) + k*np.log(occ) + log_factorial(occ*k-k) - log_factorial(occ*k-2*k+2) + k*np.log(z) - energy*(k-1) - log_factorial(k)
        count = np.exp(count)
        return count
    
    # total number of mean-particles
    def total_num (z):
        total = 0.0
        for k in range(1, N+1):
            total += k*cluster_count(k, z)
        return total
    
    # test-cluster size proportion factor
    def q (s,z):
        q = 0.0
        for b in range(1, f_single+1):
            logA = log_factorial((occ-1)*s-occ) - log_factorial(f_single-b) - log_factorial(b-1) - log_factorial(s-b-1) - log_factorial((occ-2)*s-occ+2+b)
            q += np.exp(logA)*np.exp(-(E_single-energy)*b)
        q = q* np.exp(-energy*(s-1)) * (occ**(s-1)) * (z**(s-1))
        return q

    # mean test-cluster size
    def test_cluster_size (z):
        M = 0.0
        size = 0.0
        for s in range(1, N+1):
            M += q(s,z)
            size += (s-1)*q(s,z)
        size = float(size)/M + 1
        return size

    # function has to be solved
    def func (z):
        #return total_num(z) - N_s
        return total_num(z) + test_cluster_size(z) - N

    if z == None:
        z = root_scalar(func, bracket=[10**-12, float(N)/volume]).root

    test_size_prob = [] # test-cluster size distribution    
    M = 0.0 # normalizing factor
    for s in range(1, N+1):
        prob = q(s, z)
        test_size_prob.append(prob)
        M += prob
    test_size_prob = [ float(value)/M for value in test_size_prob ]
    return test_size_prob


def mean_cluster_size (size_counts):
    total_num = 0.0
    for i in range(len(size_counts)):
        counts = size_counts[i]
        total_num += (i+1)*counts
    return total_num/sum(size_counts)

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

def mean_contact_matrix (IDs, Trace, burn_in=0):
    mean = np.zeros((len(IDs), len(IDs)))
    for i in range(burn_in, len(Trace)):
        matrix = make_contact_matrix(IDs, Trace[i])
        mean += matrix
    return mean/(len(Trace)-burn_in)

"""
fname = "single_100_23_-10:33_-4_rgs"
#fname = 'homo_100_23_-4_rca'
#fname = 'check_rgs'
ctprob_fname = fname + '_ctprob.txt'
ctnum_fname = fname + '_ctnum.txt'
#trace_fname = fname + '_trace.txt'
cluster_fname = fname + "_cluster.txt"

ctnum_list = read_ctnum(ctnum_fname)
fig = plt.figure()
plt.plot(ctnum_list)
plt.show()
plt.close()

ID_size_dist = read_cluster (cluster_fname)
#ID_size_dist = read_ctnum(cluster_fname)
fig = plt.figure()
plt.plot(range(1, len(ID_size_dist['0'])+1), ID_size_dist["0"], 'k--')
#plt.plot(range(1, len(ID_size_dist)+1), ID_size_dist, 'k--')
plt.show()
plt.close()

IDs, matrix = read_ctprob (ctprob_fname)
fig = plt.figure()
current_cmap = matplotlib.cm.get_cmap()
current_cmap.set_bad(color='white')
plt.imshow(matrix)
plt.colorbar()
plt.show()
plt.close()
"""

# pure case
N = 100
volume = 10**7

energy_list_list = []
P_list_list = []
pred_list_list = []
S_list_list = []
gel_point_list = []
for v in [2]:
    f = 3 + 10*v
    #f = 23
    gel_point = math.log(N) + math.log(f) + 2*math.log(f-2) - math.log(volume) - math.log(f-1)
    gel_point_list.append(gel_point)
    energy_list = []
    P_list = []
    pred_list = []
    S_list = []
    for i in range(10):
        E = - 2*i
        fname = 'homo_' + str(N) + '_' + str(f) + '_' + str(E)
        ID_size_dist = read_cluster (fname + "_rgs_cluster.txt")
        cluster_list = ID_size_dist['total']

        #fig = plt.figure()
        #plt.plot(cluster_list, '--')
        #plt.plot(FS_size_dist(N, f, E, volume), '.', markersize=3)
        #plt.title(str(E))
        #plt.ylim([-2,102])
        #plt.show()
        #plt.close()
        
        #cluster_list = read_ctnum (fname + "_rca_cluster.txt")
        P = number_count(cluster_list[:10])/float(N)
        #P = sum(ID_size_dist['0'][:10])
        Nc = volume*(f-1)*np.exp(E) / float(f*((f-2)**2))
        S = sum([cluster_list[k]*float(k+1) for k in range(len(cluster_list))])
        predict = 1.7*Nc/float(N)
        #predict = number_count(FS_size_dist(N, f, E, volume)[:10])/float(N)
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
#Po = float(volume)*np.exp(E-1) / float(N*f*(f-1))
Po = (volume/N) * (np.exp(E)/f) * ((f-2)**(f-2)) / ((f-1)**(f-1))
Sm = Pm + (N-Nc)*(1.0-Pm)

energy_list_list = []
P_list_list = []
pred_list_list = []
S_list_list = []
for v in range(10):
    f_single = 3 + 10*v
    #f_single = 23
    energy_list = []
    P_list = []
    pred_list = []
    S_list = []
    for i in range(10):
        E_single = - 2*i
        fname = 'single_' + str(N) + '_' + str(f) + '_' + str(E) + ':' + str(f_single) + '_' + str(E_single)
        ID_size_dist = read_cluster (fname + "_rca_cluster.txt")
        ID_energy, ID_valency = read_anot(fname + "_anot.cn")
        #ID_size_dist = read_cluster (fname + "_rgs_cluster.txt")
        #ID_energy, ID_valency = read_anot(fname + "_para.cn")
        dE = E_single - E -4.5
        Pg = number_count(ID_size_dist['total'][:10])/float(N)

        #fig = plt.figure()
        #plt.plot(range(1, len(ID_size_dist['0'])+1), ID_size_dist["0"], 'k--', label=str(E_single) + ":" + str(f_single))
        #plt.plot(range(1, len(ID_size_dist)+1), ID_size_dist, 'k--')
        #plt.legend()
        #plt.show()
        #plt.close()

        P = sum(ID_size_dist['0'][:10])
        #P = ID_size_dist['0'][0]
        S = sum([ID_size_dist['0'][k]*float(k+1) for k in range(len(ID_size_dist['0']))])
        P = (P)
        #C = ((float(f)/f_single))*((1+1.0/float(f-1))**(f-1))/((1+np.exp(-0.5*dE)/float(f-1))**(f_single-1))
        C = ((float(f)/f_single)*(float(f-1)/(f-2))**(f-1))/((1+np.exp(-0.5*dE)/float(f-2))**(f_single-1))
        #C = ((float(f)/f_single)*np.exp(0.5*dE+2))/((1+np.exp(-0.5*dE-1)/float(f-1))**(f_single-1))
        #Co = np.exp(2)/((1+np.exp(-1)/float(f-1))**(f-1))
        #predict = Pm*(Po*(1-C) + np.sqrt((Po*(1-C))**2 + 4*Po*C))/(Po*(1-Co) + np.sqrt((Po*(1-Co))**2 + 4*Po*Co))
        predict = 0.5*(Po*(1-C) + np.sqrt((Po*(1-C))**2 + 4*Po*C))
        #predict = np.exp(np.exp(0.5*dE*(f_single/float(f))))
        #predict = 1.0 / (1.0 + (float(alpha)/(1-alpha))*np.exp(-0.5*dE))**f_single
        #predict = Pm*(float(f)/f_single)*np.exp(0.5*dE)
        #predict = np.log((float(f)/f_single)*np.exp(0.5*dE))
        #predict = Pm*(float(f)/(f_single))*(np.exp(0.5*dE)*(1.0-alpha)+alpha)
        #predict = (float(f)/(f_single))*(np.exp(0.5*dE)*(1.0-alpha)+alpha)
        #predict = Pm*((float(f-1)/(f-2))**f)*((float(f-2)/(f-2 + np.exp(-0.5*dE)))**f_single) + 10**-10
        #predict = Pm*(1+f_single*((float(f-2)**(f-2))/(float(f-1)**(f-1)))*np.exp(-0.5*dE))/(1+np.exp(-0.5*dE)/float(f-2))**f_single
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

