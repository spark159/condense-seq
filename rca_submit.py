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

def random_contact_simulation (anot_fname,
                               field,
                               cycle_num,
                               burn_in,
                               NCP_st,
                               NCP_num,
                               model,
                               metric,
                               rmin,
                               rmax,
                               p_len,
                               scale,
                               short_cut,
                               occ_limit,
                               circular,
                               out_fname):

    # load anotation file
    def read_anot_file(fname, target_names=None, num_max=sys.maxint):
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
    ID_chr, ID_pos, name_ID_value = read_anot_file(anot_fname)
    ID_score = name_ID_value[field]

    # sub-sample the data
    temp = {}
    for i in range(NCP_st, NCP_st + NCP_num):
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
        prob1 = 1.0 - scale*math.exp(-score1)
        if prob1 <= 0:
            prob1 = 0.001
        if prob1 >= 1:
            prob1 = 0.999
        #assert prob1 >=0 and prob1 <= 1
        prob2 = 1.0 - scale*math.exp(-score2)
        if prob2 <= 0:
            prob2 = 0.001
        if prob2 >= 1:
            prob2 = 0.999
        #assert prob2 >=0 and prob2 <= 1
        if metric == "product":
            prob = prob1*prob2
        elif metric == "GM":
            prob = math.sqrt(prob1*prob2)
        return prob

    def Dijkstra (v1_v2_weight, source, target=None):
        dist, prev = {source:0}, {}
        Q = []
        V = v1_v2_weight.keys()
        for v in V:
            if v != source:
                dist[v] = float("inf")
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

    # random contact simulation
    print >> sys.stderr, "Random contact simulation"
    
    f = open(out_fname + '_rca_trace.txt', 'w')

    print >> f, "IDs" + ":" + ",".join(str(ID) for ID in IDs)
    print >> f, "CycleNum:" + str(cycle_num)
    
    v1_v2_weight = {} # linear genome + contact graph
    v1_v2_boolean = {} # contact graph
    for k in range(len(IDs)):
        v1 = IDs[k]
        v1_v2_weight[v1] = {}
        v1_v2_boolean[v1] = {}
        v2_list = []
        if circular and k >= 1:
            v2_list.append(IDs[k-1])
        if k < len(IDs) - 1:
            v2_list.append(IDs[k+1])
        for v2 in v2_list:
            weight = abs(ID_pos[v1] - ID_pos[v2])
            v1_v2_weight[v1][v2] = weight
        
    ID_occ = {} # occupancy number for each NCP
    for ID in IDs:
        ID_occ[ID] = 0

    i = 0
    order = -1
    while i < cycle_num:
        if int(math.log10(i+1)) > order:
            order = int(math.log10(i+1))
            print >> sys.stderr, "cycle num", i+1

        """
        # choose a random pair
        if len(IDs) >= 100:
            step = sys.maxint
            while step > len(IDs) - 1:
                step = int(half_right_normal (mu=1, sigma=len(IDs)/10))
            idx1 = random.randint(0, len(IDs) -1 - step)
            idx2 = idx1 + step
            ID1, ID2 = IDs[idx1], IDs[idx2]
        else:
            ID1, ID2 = sorted(random.sample(IDs, 2))
        """
        
        idx1, idx2 = sorted(random.sample(range(len(IDs)), 2))
        ID1, ID2 = IDs[idx1], IDs[idx2]
        
        try:
            already = v1_v2_boolean[ID1][ID2]
        except:
            already = False

        # the case of already bonded pair
        if already:
            P_rec = reaction_prob (ID_score[ID1], ID_score[ID2], metric=metric)
            if random.random() > P_rec:
                del v1_v2_weight[ID1][ID2]
                if circular:
                    del v1_v2_weight[ID2][ID1]
                if ID_index[ID2] - ID_index[ID1] <= 1: # neighbor in linear genome
                    v1_v2_weight[ID1][ID2] = ID_pos[ID2] - ID_pos[ID1]
                    if circular:
                        v1_v2_weight[ID2][ID1] = ID_pos[ID2] - ID_pos[ID1]
                del v1_v2_boolean[ID1][ID2]
                ID_occ[ID1] -= 1
                ID_occ[ID2] -= 1
                print >> f, "\n@" + str(i+1)
                print >> f, "-:" + str(ID1) + ',' + str(ID2)
            i += 1
            continue

        # the case of occupancy limt violation
        if ID_occ[ID1] >= occ_limit or ID_occ[ID2] >= occ_limit:
            i += 1
            continue

        # get the shortest distance of the chosen pair        
        dist, _ = Dijkstra(v1_v2_weight, ID1, target=ID2)
        dist = dist[ID2]*0.34

        # try collision
        if not collision_test (dist, model=model):
            i += 1
            continue        

        # try contact reaction
        P_rec = reaction_prob (ID_score[ID1], ID_score[ID2], metric=metric)
        if random.random() > P_rec:
            i += 1
            continue
        v1_v2_weight[ID1][ID2] = min(ID_pos[ID2]-ID_pos[ID1], short_cut/0.34)
        if circular:
            v1_v2_weight[ID2][ID1] = min(ID_pos[ID2]-ID_pos[ID1], short_cut/0.34)
        v1_v2_boolean[ID1][ID2] = True
        ID_occ[ID1] += 1
        ID_occ[ID2] += 1
        print >> f, "\n@" + str(i+1)
        print >> f, "+:" + str(ID1) + ',' + str(ID2)
        i += 1

    print >> sys.stderr, "Done"

    f.close()

    print >> sys.stderr, "\nCalculating mean contact matrix"

    mean_matrix = Matrix(len(IDs), len(IDs))
    D = Matrix(len(IDs), len(IDs))

    ctnum = 0
    ctnum_list = []

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
    repeat = cycle_num - max(pre_cycle, burn_in+1) + 1
    mean_matrix += D*repeat
    ctnum_list += [ctnum]*repeat
    assert len(ctnum_list) == cycle_num - burn_in
    mean_matrix = mean_matrix * (1.0/(cycle_num-burn_in))

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
    parser.add_argument(metavar='--anot',
                        dest="anot_fname",
                        type=str,
                        help='annotation file')
    parser.add_argument('--field',
                        dest="field",
                        type=str,
                        default='work/condense_seq/sp9_hg19_chr1',
                        help='condensabiltiy data field name')
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
                        dest="NCP_st",
                        type=int,
                        default=0,
                        help='start nucleosome ID')
    parser.add_argument('--num',
                        dest="NCP_num",
                        type=int,
                        default=100,
                        help='target nucleosome number')
    parser.add_argument('--model',
                        dest="model",
                        type=str,
                        default='WLC',
                        help='Polymer model (WLC or RW)')
    parser.add_argument('--metric',
                        dest="metric",
                        type=str,
                        default='product',
                        help='probability metric')
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
    parser.add_argument('--scale',
                        dest="scale",
                        type=float,
                        default=0.413,
                        help='global survivial probabilty')
    parser.add_argument('--short-cut',
                        dest="short_cut",
                        type=float,
                        help='short cut length by contact (nm)')
    parser.add_argument('--occ',
                        dest="occ_limit",
                        type=int,
                        default=4,
                        help='occupancy number limit')
    parser.add_argument('--cir',
                        dest="circular",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='circular graph option')    
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

    print >> sys.stderr, "simulation parameters"
    print >> sys.stderr, "field: " + args.field
    print >> sys.stderr, "NCP: " + str(args.NCP_st) + '-' + str(args.NCP_st + args.NCP_num - 1)
    print >> sys.stderr, "cycle number: " + str(args.cycle_num)
    print >> sys.stderr, "scale: " + str(args.scale)
    print >> sys.stderr, "short cut length: " + str(short_cut)
    print >> sys.stderr, "occupancy number: " + str(args.occ_limit)
    print >> sys.stderr, "circular choice: " + str(args.circular)
    print >> sys.stderr

    random_contact_simulation (args.anot_fname,
                               args.field,
                               args.cycle_num,
                               burn_in,
                               args.NCP_st,
                               args.NCP_num,
                               args.model,
                               args.metric,
                               args.rmin,
                               args.rmax,
                               args.p_len,
                               args.scale,
                               short_cut,
                               args.occ_limit,
                               args.circular,
                               args.out_fname)
