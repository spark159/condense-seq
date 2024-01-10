import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import random
import copy
import heapq
import scipy.special as special
from scipy.integrate import quad

def tuple_cmp (a, b):
    if a[0] < b[0]:
        return -1
    elif a[0] > b[0]:
        return 1
    else:
        if a[1] < b[1]:
            return -1
        return 1

# bisection search
def bisect (sL, target):
    st, ed = 0, len(sL)
    while st != ed:
        mid = (st + ed) / 2
        if target < sL[mid]:
            ed = mid
        elif target > sL[mid]:
            st = mid + 1
        elif target == sL[mid]:
            return mid
    return None

def half_right_normal (mu, sigma):
    value = random.gauss(mu, sigma)
    if value < mu:
        return 2*mu - value
    return value

class Matrix (object):
    def __init__ (self, m, n, value=0.0):
        self.m=m; self.n=n;
        self.matrix=[]
        for i in range(m):
            self.matrix.append([value]*n)
    
    def __getitem__ (self, idx):
        return self.matrix[idx]
    
    def __setitem__ (self, idx, value):
        self.matrix[idx]=value
    
    def getrank (self):
        return (self.m,self.n)
        
    def __add__ (self, other):
        if type(other) == Matrix and self.getrank() == other.getrank():
            new_one = Matrix(self.m,self.n)
            for i in range(self.m):
                for j in range(self.n):
                    new_one.matrix[i][j] = self.matrix[i][j] + other.matrix[i][j]
        elif type(other) == int or type(other) == float:
            new_one = Matrix(self.m,self.n)
            for i in range(self.m):
                for j in range(self.n):
                    new_one[i][j] = self.matrix[i][j] + other
        else:
            print "Not possible"
            return None
        return new_one
    
    def __sub__ (self, other):
        if type(other) == Matrix and self.getrank() == other.getrank():
            new_one = Matrix(self.m,self.n)
            for i in range(self.m):
                for j in range(self.n):
                    new_one.matrix[i][j] = self.matrix[i][j] - other.matrix[i][j]
        elif type(other) == int or type(other) == float:
            new_one = Matrix(self.m,self.n)
            for i in range(self.m):
                for j in range(self.n):
                    new_one[i][j] = self.matrix[i][j] - other
        else:
            print "Not possible"
            return None
        return new_one

    def __mul__ (self, other):
        if type(other) == Matrix and self.n == other.m:
            new_one = Matrix(self.m, other.n)
            for i in range(self.m):
                for j in range(other.n):
                    value = 0
                    for k in range(self.n):
                        value += self.matrix[i][k]*other.matrix[k][j]
                    new_one.matrix[i][j] = value
        elif type(other) == int or type(other) == float:
            new_one = Matrix(self.m,self.n)
            for i in range(self.m):
                for j in range(self.n):
                    new_one[i][j] = self.matrix[i][j]*other
        else:
            print "Not possible"
            return None
        return new_one
    
    def __str__ (self):
        s=''
        for i in range(self.m):
            s += str(self.matrix[i]) + '\n'
        return s[0:len(s)-1]

def random_graph_simulation (para_fname,
                             field,
                             cycle_num,
                             burn_in,
                             PTL_st,
                             PTL_num,
                             model,
                             metric,
                             rmin,
                             rmax,
                             p_len,
                             volume,
                             expN,
                             scale,
                             short_cut,
                             circular,
                             bethe,
                             indis,
                             out_fname):

    # load parameter file
    def read_para_file(fname, target_names=None, num_max=sys.maxint):
        ID_chr, ID_pos = {}, {}
        name_ID_value = {}
        First = True
        count = 0
        for line in open(fname):
            if count > num_max:
                break
            cols = line.strip().split()
            if First:
                names = cols[3:]
                First = False
                continue
            ID, chr, pos = int(cols[0]), cols[1], int(cols[2])
            ID_chr[ID] = chr
            ID_pos[ID] = pos
            cols = cols[3:]
            for i in range(len(cols)):
                name = names[i]
                if target_names and name not in target_names:
                    continue
                if name not in name_ID_value:
                    name_ID_value[name] = {}
                assert ID not in name_ID_value[name]
                try:
                    value = float(cols[i])
                except:
                    value = cols[i]
                name_ID_value[name][ID] = value
            count += 1
        return ID_chr, ID_pos, name_ID_value
    ID_chr, ID_pos, name_ID_value = read_para_file(para_fname)
    ID_score = name_ID_value[field]
    ID_valency = name_ID_value["Valency"]

    # sub-sample the data
    temp = {}
    for i in range(PTL_st, PTL_st + PTL_num):
        temp[i] = ID_score[i]
    ID_score = temp

    # sort by position
    pos_ID = []
    for ID in ID_score:
        pos = ID_pos[ID]
        pos_ID.append([pos, ID])
    pos_ID = sorted(pos_ID, cmp=tuple_cmp)

    IDs, ID_index = [], {}
    for i in range(len(pos_ID)):
        pos, ID = pos_ID[i]
        IDs.append(ID)
        ID_index[ID] = i

    # Random walk chain simulation
    def random_walk (N, b):
        assert N >= 1 and b > 0
        R = [0.0, 0.0, float(b)]
        for i in range(N-1):
            u, v = random.uniform(-1, 1), random.random()
            x = b*math.sqrt(1-u**2)*math.cos(2*math.pi*v)
            y = b*math.sqrt(1-u**2)*math.sin(2*math.pi*v)
            z = b*u
            R[0] +=x; R[1] +=y; R[2] +=z    
        endToend = 0.0
        for k in range(3):
            endToend += R[k]**2
        return math.sqrt(endToend)

    # Random walk chain radial function
    def RW_radial (r, L, p_len):
        b = 2*p_len # kuhn length = 2 * persistence length        
        prob = ((3.0 / (2*math.pi*L*b))**(3.0/2))*math.exp(-3*(r**2)/ float(2*L*b))
        return 4*math.pi*r*r*prob

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
            JSY *= 896.32*(lamb**5)*math.exp(-14.054*lamb + 0.246/lamb)
        else:
            JSY *= ((3.0*lamb/math.pi)**(3.0/2))*(1.0 - 5.0*lamb/4.0 - 79.0*lamb*lamb/160.0)

        term1 = 0.0
        for i in [-1, 0]:
            for j in [1, 2, 3]:
                term1 += C[i+1][j-1]*(lamb**i)*(rho**(2*j))
        term2 = (-d*lamb*a*(1+b)*rho)/(1.0-b*b*rho*rho)

        prob = JSY*(((1.0-c*rho*rho)/(1.0-rho*rho))**(5.0/2))*math.exp(term1/(1.0-rho*rho))
        prob *= math.exp(term2*b*rho)
        prob *= special.iv(0, term2)

        return 4*math.pi*r*r*prob

    def collision_prob (dist, model, rmin=rmin, rmax=rmax, p_len=p_len):
        if model == 'Free':
            if dist == "not connected":
                return 1.0 / float(volume)
            else:
                return 1.0
            
        if dist <= 0:
            if dist >= rmin and dist <= rmax:
                return 1.0
            return 0.0

        if model == "RW":
            if float(dist) / p_len < 10:
                b = 2*p_len
                N = int(round(float(dist)/b))
                if N < 1:
                    endToend = dist
                else:
                    endToend = random_walk(N, b)
                if endToend >= rmin and endToend <= rmax:
                    return 1.0   # Need to change (not probability)
                return 0.0 # Need to change (not probability)
            else:
                return quad(RW_radial, min(rmin, dist), min(rmax, dist), args=(dist, p_len))[0] 

        elif model == "WLC":
            return quad(WLC_radial, min(rmin, dist), min(rmax, dist), args=(dist, p_len))[0]

    # convert scores to the interaction energy
    def scores_to_energy (score1, score2, metric, scale=scale):
        if metric == "raw":
            return 0.5*(score1+score2)
        prob1 = scale*math.exp(-score1)
        if prob1 <= 0:
            prob1 = 0.00001
        if prob1 >= 1:
            prob1 = 0.99999
        prob2 = scale*math.exp(-score2)
        if prob2 <= 0:
            prob2 = 0.00001
        if prob2 >= 1:
            prob2 = 0.99999
        if metric == "product":
            return 0.5*math.log(prob1*prob2)
        if metric == "F-S": # need better formula for inhomogenous particles case
            structure = (occ_limit*(occ_limit-2)**2) / float(occ_limit-1)
            expN1 = expN * (10**-7) #Need to know frac1
            expN2 = expN * (10**-7) #Need to know frac2
            factor1 = structure * expN1 / float(volume)
            factor2 = structure * expN2 / float(volume)
            return 0.5*math.log(factor1*prob1*factor2*prob2)

    def reaction_prob (score1, score2, metric, scale=scale):
        energy = scores_to_energy(score1, score2, metric=metric, scale=scale)
        return math.exp(-energy) / float(1.0 + math.exp(-energy))

    def Dijkstra (v1_v2_weights, source, target=None):
        dist, prev = {source:0}, {}
        Q = []
        V = v1_v2_weights.keys()
        for v in V:
            if v != source:
                dist[v] = float("inf")
            prev[v] = None
            heapq.heappush(Q, (dist[v], v))
        while len(Q) > 0:
            _, u = heapq.heappop(Q)
            if target and u == target:
                return dist, prev
            for v in v1_v2_weights[u]:
                alt = dist[u] + min(v1_v2_weights[u][v].values())
                if alt < dist[v]:
                    Q.remove((dist[v],v))
                    heapq.heapify(Q)
                    heapq.heappush(Q, (alt, v))
                    dist[v] = alt
                    prev[v] = u
        return dist, prev

    def BFS (ID_neighbors, source, target=None):
        visited = {source:True}
        queue = [source]
        while len(queue) > 0 :
            v = queue.pop(0)
            if target and v == target:
                return True
            try:
                neighbors = list(ID_neighbors[v])
            except:
                continue
            for w in neighbors:
                try:
                    visited[w]
                    continue
                except:
                    visited[w] = True
                    queue.append(w)
        if target:
            return False
        return set(visited.keys())

    def clustering (ID_neighbors):
        size_clusters = {}
        IDs = set(ID_neighbors.keys())
        while len(IDs) > 0:
            ID = IDs.pop()
            cluster = BFS(ID_neighbors, ID)
            size = len(cluster)
            if size not in size_clusters:
                size_clusters[size] = []
            size_clusters[size].append(cluster)
            IDs = IDs - cluster
        return size_clusters
    

    # random graph simulation
    print >> sys.stderr, "Random graph simulation"
    
    f = open(out_fname + '_rca_trace.txt', 'w')

    print >> f, "IDs" + ":" + ",".join(str(ID) for ID in IDs)
    print >> f, "CycleNum:" + str(cycle_num)

    v1_v2_weights = {} # genome + contact distance graph
    for k in range(len(IDs)):
        v1 = IDs[k]
        v1_v2_weights[v1] = {}
        v2_list = []
        if circular and k >= 1:
            v2_list.append(IDs[k-1])
        if k < len(IDs) - 1:
            v2_list.append(IDs[k+1])
        for v2 in v2_list:
            v1_v2_weights[v1][v2] = {}
            v1_v2_weights[v1][v2]['genome'] = abs(ID_pos[v1] - ID_pos[v2])

    ID_neighbors = {ID:set([]) for ID in IDs} # contact neighbors
        
    i = 0
    order = -1
    while i < cycle_num:
        if int(math.log10(i+1)) > order:
            order = int(math.log10(i+1))
            print >> sys.stderr, "cycle num", i+1
        
        ## choose a random pair
        #if len(IDs) >= 100:
        #    step = sys.maxint
        #    while step > len(IDs) - 1:
        #        step = int(half_right_normal (mu=1, sigma=len(IDs)/10))
        #    idx1 = random.randint(0, len(IDs) -1 - step)
        #    idx2 = idx1 + step
        #    ID1, ID2 = IDs[idx1], IDs[idx2]
        #else:
        #    ID1, ID2 = sorted(random.sample(IDs, 2))
        
        
        idx1, idx2 = sorted(random.sample(range(len(IDs)), 2))
        ID1, ID2 = IDs[idx1], IDs[idx2]

        try:
            v1_v2_weights[ID1][ID2]['contact']
            already = True
        except:
            already = False

        # the case of already bonded pair
        if already:
            # break the contact for a moment
            del v1_v2_weights[ID1][ID2]['contact']
            if circular:
                del v1_v2_weights[ID2][ID1]['contact']
            ID_neighbors[ID1].remove(ID2)
            ID_neighbors[ID2].remove(ID1)

            # check the connectivity (Free) or the shortest distance (Polymer) of the chosen pair
            if model == 'Free':
                if bethe or len(ID_neighbors[ID1]) <= 0 or len(ID_neighbors[ID2]) <= 0 or not BFS(ID_neighbors, ID1, target=ID2):
                    dist = "not connected"
                else:
                    dist = "connected"
            else:
                dist, _ = Dijkstra(v1_v2_weights, ID1, target=ID2)
                dist = dist[ID2]*0.34

            # try break the contact
            P_col = collision_prob (dist, model=model)
            energy = scores_to_energy (ID_score[ID1], ID_score[ID2], metric=metric)
            if indis: # indistingushable 
                conf = 1
            else: # distingushable
                #conf = (occ_limit - len(ID_neighbors[ID1])) * (occ_limit - len(ID_neighbors[ID2]))
                conf = (ID_valency[ID1] - len(ID_neighbors[ID1])) * (ID_valency[ID2] - len(ID_neighbors[ID2]))
            if random.random() > min(1, math.exp(energy)/float(P_col*conf)):
                v1_v2_weights[ID1][ID2]['contact'] = short_cut/0.34
                if circular:
                    v1_v2_weights[ID2][ID1]['contact'] = short_cut/0.34
                ID_neighbors[ID1].add(ID2)
                ID_neighbors[ID2].add(ID1)
                i += 1
                print i
                continue
            
            print >> f, "\n@" + str(i+1)
            print >> f, "-:" + str(ID1) + ',' + str(ID2)
            i += 1
            print i
            continue

        # occupancy limt violation
        #if len(ID_neighbors[ID1]) >= occ_limit or len(ID_neighbors[ID2]) >= occ_limit:
        if len(ID_neighbors[ID1]) >= ID_valency[ID1] or len(ID_neighbors[ID2]) >= ID_valency[ID2]:
            continue

        # Bethe lattice violation
        if bethe and  len(ID_neighbors[ID1]) > 0 and len(ID_neighbors[ID2]) > 0:
            if BFS (ID_neighbors, ID1, target=ID2):
                continue

        # check the connectivity (Free) or the shortest distance (Polymer) of the chosen pair
        if model == 'Free':
            if bethe or len(ID_neighbors[ID1]) <= 0 or len(ID_neighbors[ID2]) <= 0 or not BFS(ID_neighbors, ID1, target=ID2):
                dist = "not connected"
            else:
                dist = "connected"
        else:
            dist, _ = Dijkstra(v1_v2_weights, ID1, target=ID2)
            dist = dist[ID2]*0.34
            
        # try contact reaction
        P_col = collision_prob (dist, model=model)
        energy = scores_to_energy (ID_score[ID1], ID_score[ID2], metric=metric)
        if indis: # indistingushable 
            conf = 1
        else: # distingushable
            #conf = (occ_limit - len(ID_neighbors[ID1])) * (occ_limit - len(ID_neighbors[ID2]))
            conf = (ID_valency[ID1] - len(ID_neighbors[ID1])) * (ID_valency[ID2] - len(ID_neighbors[ID2]))
        if random.random() > min(1, P_col * conf * math.exp(-energy)):
            i += 1
            print i
            continue

        if ID2 not in v1_v2_weights[ID1]:
            v1_v2_weights[ID1][ID2] = {}
        v1_v2_weights[ID1][ID2]['contact'] = short_cut/0.34
        if circular:
            if ID1 not in v1_v2_weights[ID2]:
                v1_v2_weights[ID2][ID1] = {}
            v1_v2_weights[ID2][ID1]['contact'] = short_cut/0.34
        ID_neighbors[ID1].add(ID2)
        ID_neighbors[ID2].add(ID1)
        print >> f, "\n@" + str(i+1)
        print >> f, "+:" + str(ID1) + ',' + str(ID2)
        i += 1
        print i

    print >> sys.stderr, "Done"
    f.close()

    print >> sys.stderr, "\nAveraging the data"

    # mean contact matrix
    mean_matrix = Matrix(len(IDs), len(IDs))
    D = Matrix(len(IDs), len(IDs))

    # total contact number
    ctnum = 0
    ctnum_list = []

    # cluster information
    ID_neighbors = {ID:set([]) for ID in IDs} 
    cID_IDs = {'-'.join([str(1), str(IDs[i])]):set([IDs[i]]) for i in range(len(IDs))} # cluster ID to particle IDs
    ID_cID = {IDs[i]:'-'.join([str(1), str(IDs[i])]) for i in range(len(IDs))} # particle ID to cluster ID
    mean_size_num = [0]*len(IDs) # cluster size distribution
    ID_mean_size_num = {ID:[0]*len(IDs) for ID in IDs} # cluster size distribution by particle
    
    ## mean cluster size distribution
    #mean_size_num = [0]*len(IDs)
    #size_num = [len(IDs)] + [0]*(len(IDs)-1)

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
                    
    
    #def cluster_update (ID_neighbors, size_num, ID1, ID2, type):
    #    if type == "+":
    #        cluster1 = BFS(ID_neighbors, ID1)
    #        if ID2 not in cluster1:
    #            cluster2 = BFS(ID_neighbors, ID2)
    #        else:
    #            cluster2 = cluster1
    #        size1, size2 = len(cluster1), len(cluster2)
    #        if cluster1 != cluster2:
    #            size_num[size1-1] -=1
    #            size_num[size2-1] -=1
    #            size_num[size1+size2-1] +=1
    #        ID_neighbors[ID1].add(ID2)
    #        ID_neighbors[ID2].add(ID1)
    #    else:
    #        ID_neighbors[ID1].remove(ID2)
    #        ID_neighbors[ID2].remove(ID1)
    #        cluster1 = BFS(ID_neighbors, ID1)
    #        if ID2 not in cluster1:
    #            cluster2 = BFS(ID_neighbors, ID2)
    #        else:
    #            cluster2 = cluster1
    #        size1, size2 = len(cluster1), len(cluster2)
    #        if cluster1 != cluster2:
    #            size_num[size1+size2-1] -=1
    #            size_num[size1-1] +=1
    #            size_num[size2-1] +=1
    #    return

    # reading trace file 
    pre_cycle = 1
    for line in open(out_fname + '_rca_trace.txt', 'r'):
        line = line.strip()
        if not line:
            continue
        if line.startswith('IDs'):
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
    f = open(out_fname + '_rca_mean.txt', 'w')
    print >> f, "IDs" + ":" + ",".join(str(ID) for ID in IDs)

    m, n = mean_matrix.getrank()
    for i in range(m):
        print >> f, ",".join(str(value) for value in mean_matrix[i])

    f.close()

    f = open(out_fname + '_rca_ctnum.txt', 'w')
    print >> f, "CycleNum:" + str(cycle_num)
    print >> f, ",".join(str(value) for value in ctnum_list)

    f.close()

    f = open(out_fname + '_rca_cluster.txt', 'w')
    print >> f, "IDs" + ":" + ",".join(str(ID) for ID in IDs)
    print >> f, "@total"
    print >> f, ",".join(str(value) for value in mean_size_num)
    for ID in IDs:
        print >> f, "@" + str(ID)
        print >> f, ",".join(str(value) for value in ID_mean_size_num[ID])

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

    parser = ArgumentParser(description='Random contact simulation')
    parser.add_argument(metavar='--para',
                        dest="para_fname",
                        type=str,
                        help='parameter file')
    parser.add_argument('--field',
                        dest="field",
                        type=str,
                        default='Score',
                        help='data field name')
    parser.add_argument('--cycle',
                        dest="cycle_num",
                        type=int,
                        default=10000,
                        help='iteration cycle number')
    parser.add_argument('--burn',
                        dest="burn_in",
                        type=int,
                        help='burn-in cycle number')
    parser.add_argument('--st',
                        dest="PTL_st",
                        type=int,
                        default=0,
                        help='start particle ID')
    parser.add_argument('--num',
                        dest="PTL_num",
                        type=int,
                        default=100,
                        help='simulation particle number')
    parser.add_argument('--model',
                        dest="model",
                        type=str,
                        default='WLC',
                        help='free solution(Free) or polymer model(WLC/RW)')
    parser.add_argument('--metric',
                        dest="metric",
                        type=str,
                        default='F-S',
                        help='interaction energy metric')
    parser.add_argument('--rmin',
                        dest="rmin",
                        type=float,
                        default=0,
                        help='minimum collision distance (nm)')
    parser.add_argument('--rmax',
                        dest="rmax",
                        type=float,
                        default=50,
                        help='maximum collision distance (nm)')
    parser.add_argument('--plen',
                        dest="p_len",
                        type=float,
                        default=50,
                        help='persistence length of dsDNA (nm)')
    parser.add_argument('--vfrac',
                        dest="vfrac",
                        type=float,
                        default=10**-5,
                        help='volume fraction (volume of particle / volume of solution)')
    parser.add_argument('--expN',
                        dest="expN",
                        type=float,
                        default=10**11,
                        help='total particle number in the experiment')
    parser.add_argument('--scale',
                        dest="scale",
                        type=float,
                        default=0.413,
                        help='global survivial probabilty')
    parser.add_argument('--short-cut',
                        dest="short_cut",
                        type=float,
                        help='short cut length by contact (nm)')
    parser.add_argument('--cir',
                        dest="circular",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='circular graph option')
    parser.add_argument('--bethe',
                        dest="bethe",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='Bethe-lattice constraint')
    parser.add_argument('--indis-site',
                        dest="indis",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='indistinguishable binding sites')    
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()

    if args.burn_in == None:
        burn_in = int(args.cycle_num*0.1) # trash first 10% of total cycle
    else:
        burn_in = args.burn_in
    
    if args.short_cut == None:
        short_cut = 0.5*(args.rmin + args.rmax) # short cut length is mean collision distance 
    else:
        short_cut = args.short_cut

    if args.model == 'Free':
        volume = args.PTL_num / float(args.vfrac) # simulation solution volume (in the unit of particle volume)
    else:
        volume = args.expN / float(args.vfrac) # experiment solution volume (in the unit of particle volume)

    print >> sys.stderr, "simulation parameters"
    print >> sys.stderr, "field: " + args.field
    print >> sys.stderr, "particles: " + str(args.PTL_st) + '-' + str(args.PTL_st + args.PTL_num - 1)
    print >> sys.stderr, "cycle number: " + str(args.cycle_num)
    print >> sys.stderr, "scale: " + str(args.scale)
    print >> sys.stderr, "collision model : " + str(args.model)
    if args.model == 'Free':
        print >> sys.stderr, "simulation volume: " + str(volume)
    else:
        print >> sys.stderr, "persistence length: " + str(args.p_len)
        print >> sys.stderr, "short cut length: " + str(short_cut)
        print >> sys.stderr, "circular choice: " + str(args.circular)
    print >> sys.stderr, "Bethe lattice constraint: " + str(args.bethe)
    print >> sys.stderr, "indistinguishable binding sites: " + str(args.indis)
    print >> sys.stderr

    random_graph_simulation (args.para_fname,
                             args.field,
                             args.cycle_num,
                             burn_in,
                             args.PTL_st,
                             args.PTL_num,
                             args.model,
                             args.metric,
                             args.rmin,
                             args.rmax,
                             args.p_len,
                             volume,
                             args.expN,
                             args.scale,
                             short_cut,
                             args.circular,
                             args.bethe,
                             args.indis,
                             args.out_fname)
