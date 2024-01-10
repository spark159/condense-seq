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
import load_file


# load file
path = "/home/spark159/../../media/spark159/sw/dataforcondense//hg19_chr1_171_everything_anot.cn"
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path)
ID_score = name_ID_value['work/condense_seq/sp9_hg19_chr1']
temp = {}
for i in range(0, 10):
    temp[i] = ID_score[i]
ID_score = temp
ID_AT = name_ID_value['ATcontent']

# parameters
rmin, rmax = 0, 50 # minimum and maximum reaction distance (nm)
p_len = 50 # DNA persistence length (nm)
scale = 0.413 # total count ratio of test sample to control (UV data, sp9)
#scale = 0.154 # total count ratio of test sample to control (UV data, sp10)
#occ_limit = sys.maxint # occupancy number limit
occ_limit = 4

def tuple_cmp (a, b):
    if a[0] < b[0]:
        return -1
    elif a[0] > b[0]:
        return 1
    else:
        if a[1] < b[1]:
            return -1
        return 1

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

def collision_test (dist, model, rmin=rmin, rmax=rmax, p_len=p_len):
    if dist <= 0:
        return dist >= rmin and dist <= rmax
    
    if model == "RW":
        if float(dist) / p_len < 10:
            b = 2*p_len
            N = int(round(float(dist)/b))
            if N < 1:
                endToend = dist
            else:
                endToend = random_walk(N, b)
            if endToend >= rmin and endToend <= rmax:
                return True
            return False
        else:
            p = quad(RW_radial, min(rmin, dist), min(rmax, dist), args=(dist, p_len))[0]
            if random.random() < p:
                return True
            return False
        
    elif model == "WLC":
        p = quad(WLC_radial, min(rmin, dist), min(rmax, dist), args=(dist, p_len))[0]
        if random.random() < p:
            return True
        return False
    
def reaction_prob (score1, score2, metric, scale=scale):
    prob1 = 1.0 - scale*np.exp(-score1)
    if prob1 <= 0:
        prob1 = 0.001
    if prob1 >= 1:
        prob1 = 0.999
    #assert prob1 >=0 and prob1 <= 1
    prob2 = 1.0 - scale*np.exp(-score2)
    if prob2 <= 0:
        prob2 = 0.001
    if prob2 >= 1:
        prob2 = 0.999
    #assert prob2 >=0 and prob2 <= 1
    if metric == "product":
        prob = prob1*prob2
    elif metric == "GM":
        prob = np.sqrt(prob1*prob2)
    return prob

def Dijkstra(Vertexes, v1_v2_weight, source, target=None):
    dist, prev = {source:0}, {}
    Q = []
    for v in Vertexes:
        if v != source:
            dist[v] = np.inf
        prev[v] = None
        heapq.heappush(Q, (dist[v], v))
    while len(Q) > 0:
        _, u = heapq.heappop(Q)
        if target and u == target:
            return dist, prev
        for v in v1_v2_weight[u]:
            alt = dist[u] + v1_v2_weight[u][v]
            if alt < dist[v]:
                Q.remove((dist[v],v))
                heapq.heapify(Q)
                heapq.heappush(Q, (alt, v))
                dist[v] = alt
                prev[v] = u
    return dist, prev
        
def random_contact_algorithm (ID_score, cycle_num, model, metric):
    Trace = []
    i = 0
    order = -1
    while i < cycle_num:
        if int(np.log10(i+1)) > order:
            order = int(np.log10(i+1))
            print "cycle num", i+1
        
        ID_contacts = {}
        if len(Trace) > 0:
            ID_contacts = copy.deepcopy(Trace[-1])
            
        ID1, ID2 = sorted(random.sample(ID_score.keys(), 2))

        contacts1, contacts2 = [], []
        if ID1 in ID_contacts:
            contacts1 = ID_contacts[ID1]
        if ID2 in ID_contacts:
            contacts2 = ID_contacts[ID2]
            
        # the case of already bonded pair
        if ID1 in contacts2 and ID2 in contacts1:
            P_rec = reaction_prob (ID_score[ID1], ID_score[ID2], metric=metric)
            if random.random() > P_rec:
                contacts1.remove(ID2)
                contacts2.remove(ID1)
                if len(contacts1) <= 0:
                    del ID_contacts[ID1]
                else:
                    ID_contacts[ID1] = contacts1
                if len(contacts2) <= 0:
                    del ID_contacts[ID2]
                else:
                    ID_contacts[ID2] = contacts2
            Trace.append(ID_contacts)
            i += 1
            continue

        # the case of occupancy limt violation
        if len(contacts1) >= occ_limit or len(contacts2) >= occ_limit:
            Trace.append(ID_contacts)
            i += 1
            continue

        # get the shortest distance of the chosen pair
        v1_v2_weight = {}
        for ID in range(ID1, ID2):
            v1 = ID
            if v1 not in v1_v2_weight:
                v1_v2_weight[v1] = {}
            v2 = ID + 1
            v1_v2_weight[v1][v2] = ID_pos[v2] - ID_pos[v1]
            try:
                for v2 in ID_contacts[v1]:
                    if v2 > v1 and v2 <= ID2:
                        v1_v2_weight[v1][v2] = 0
                        #v1_v2_weight[v1][v2] = min(ID_pos[v2] - ID_pos[v1], rmax/0.34)
            except:
                None
        dist, _ = Dijkstra(range(ID1, ID2+1), v1_v2_weight, ID1, target=ID2)
        dist = dist[ID2]*0.34
        #dist = 0.34*(ID_pos[ID2] - ID_pos[ID1])
        
        # try collision
        if not collision_test (dist, model=model):
            Trace.append(ID_contacts)
            i += 1
            continue        

        # try contact reaction
        P_rec = reaction_prob (ID_score[ID1], ID_score[ID2], metric=metric)
        if random.random() > P_rec:
            Trace.append(ID_contacts)
            i += 1
            continue
        contacts1.append(ID2)
        contacts2.append(ID1)
        ID_contacts[ID1], ID_contacts[ID2] = sorted(contacts1), sorted(contacts2)
        Trace.append(ID_contacts)
        i += 1
        
    assert len(Trace) == cycle_num
    print "Done"
    return Trace
Trace = random_contact_algorithm (ID_score, cycle_num=1000000, model="WLC", metric="product")
        
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
#Video(ID_score.keys(), Trace, 1000)

def mean_contact_matrix (IDs, Trace, burn_in=0):
    mean = np.zeros((len(IDs), len(IDs)))
    for i in range(burn_in, len(Trace)):
        matrix = make_contact_matrix(IDs, Trace[i])
        mean += matrix
    return mean/(len(Trace)-burn_in)
mean_matrix = mean_contact_matrix(ID_score.keys(), Trace, burn_in=0)

fig = plt.figure()
plt.imshow(mean_matrix, cmap='Reds')
plt.colorbar()
plt.title("Mean Contact Probability")
plt.show()
plt.close()

"""
def Read_trace(trace_fname):
    Trace = []
    for line in open(trace_fname):
        line = line.strip()
        if not line:
            continue
        if line.startswith('@'):
            Trace.append(ID_contacts)
            ID_contacts = {}
            cycle_num = int(line[1:])
            continue
        ID, contacts = line.split(':')
        ID = int(ID)
        contacts = [int(e) for e in contacts.split(',')]
        ID_contacts[ID] = contacts
    assert len(Trace) == cycle_num
    return Trace

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
