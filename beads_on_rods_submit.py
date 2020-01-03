import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import random
import copy

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

def beads_on_rods (anot_fname,
                   field,
                   cycle_num,
                   burn_in,
                   sigma,
                   NCP_st,
                   NCP_num,
                   model,
                   metric,
                   rmin,
                   rmax,
                   p_len,
                   scale,
                   occ_limit,
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

    def ID_to_index (IDs, b):
        Pos_list = [ID_pos[ID] for ID in IDs]
        minPos, maxPos = min(Pos_list), max(Pos_list)
        ID_index, index_ID = {}, {}
        for ID in IDs:
            pos = ID_pos[ID]
            index = int(round(0.34*(pos - minPos)/b))
            assert index >= 0  and index < N
            ID_index[ID] = index
            index_ID[index] = ID
        return ID_index, index_ID

    def norm (v):
        total = 0.0
        for i in range(3):
            total += v[i]**2
        return math.sqrt(total)

    def normalize (v):
        new_v = [v[0]/norm(v), v[1]/norm(v), v[2]/norm(v)]
        return new_v

    def cross (v1, v2):
        x = v1[1]*v2[2] - v1[2]*v2[1]
        y = v1[2]*v2[0] - v1[0]*v2[2]
        z = v1[0]*v2[1] - v1[1]*v2[0]
        return [x, y, z]

    def distance (v1, v2):
        total = 0.0
        for i in range(3):
            total += (v1[i]-v2[i])**2
        return math.sqrt(total)

    def get_angle (v1, v2):
        dot_product = 0.0
        for i in range(3):
            dot_product += v1[i]*v2[i]
        cos = dot_product/(norm(v1)*norm(v2))
        return math.acos(cos)

    def polar_to_cartesian (theta, phi):
        x = math.sin(phi)*math.cos(theta)
        y = math.sin(phi)*math.sin(theta)
        z = math.cos(phi)
        return [x, y, z]

    def cartesian_to_polar (x, y, z):
        rho = math.sqrt(x**2 + y**2 + z**2)
        r = math.sqrt(x**2 + y**2)
        theta = math.arccos(float(z)/rho)
        phi = math.arccos(float(x)/r)
        return theta, phi

    def rotation (u, angle, p):
        def h_product(q1, q2):
            x = q1[0]*q2[0]-q1[1]*q2[1]-q1[2]*q2[2]-q1[3]*q2[3]
            y = q1[0]*q2[1]+q1[1]*q2[0]+q1[2]*q2[3]-q1[3]*q2[2]
            z = q1[0]*q2[2]-q1[1]*q2[3]+q1[2]*q2[0]+q1[3]*q2[1]
            w = q1[0]*q2[3]+q1[1]*q2[2]-q1[2]*q2[1]+q1[3]*q2[0]
            return [x, y, z, w]
        q = [math.cos(float(angle)/2), u[0]*math.sin(float(angle)/2), u[1]*math.sin(float(angle)/2), u[2]*math.sin(float(angle)/2)]
        q_inv = [math.cos(float(angle)/2), -u[0]*math.sin(float(angle)/2), -u[1]*math.sin(float(angle)/2), -u[2]*math.sin(float(angle)/2)]
        p = [0, p[0], p[1], p[2]]
        new_p = h_product(h_product(q, p), q_inv) 
        return new_p[1:]

    def pick_on_sphere ():
        u, v = random.uniform(-1, 1), random.random()
        x = math.sqrt(1-u**2)*math.cos(2*math.pi*v)
        y = math.sqrt(1-u**2)*math.sin(2*math.pi*v)
        z = u
        return [x, y, z]

    def update (old_u, sigma):
        theta = random.gauss(0, sigma)
        dum = copy.deepcopy(old_u)
        dum[0] += 1
        r = cross(dum, old_u)
        r[0] *= math.tan(theta)/norm(r); r[1] *= math.tan(theta)/norm(r); r[2] *= math.tan(theta)/norm(r)
        v = [old_u[0]+r[0], old_u[1]+r[1], old_u[2]+r[2]]
        phi = random.uniform(0, 2*math.pi)
        new_v = rotation(old_u, phi, v)
        new_u = [ float(value)/norm(new_v) for value in new_v ]
        return new_u

    def go_through (U, NCP_index):
        ete = [0.0, 0.0, 0.0]
        R = []
        ete_list = []
        NCP_index = sorted(NCP_index, reverse=True)
        for i in range(len(U)):
            u = U[i]
            if NCP_index and i == NCP_index[-1]:
                NCP_index.pop()
                R.append(copy.deepcopy(ete))
            ete[0] += b*u[0]; ete[1] += b*u[1]; ete[2] += b*u[2]
            ete_list.append(copy.deepcopy(ete))
        assert len(NCP_index) == 0
        return ete_list, R

    def bind_energy (score1, score2, metric, scale=scale):
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
        Ebind = math.log(prob/(1-prob)) # unit of KT
        return Ebind

    def energy (U, index_ID, model, metric, contact_info=False):
        total = 0.0
        ete = [0.0, 0.0, 0.0]
        temp_index = sorted(index_ID.keys(), reverse=True)
        R = []
        for i in range(len(U)):
            if model == 'WLC' and i < len(U)-1:
                u, u_next = U[i], U[i+1]
                angle = get_angle(u, u_next)
                total += 0.5*k*angle*angle
            if temp_index and i == temp_index[-1]:
                temp_index.pop()
                R.append(copy.deepcopy(ete))
            ete[0] += b*u[0]; ete[1] += b*u[1]; ete[2] += b*u[2]
        assert len(temp_index) == 0
        NCP_index = sorted(index_ID.keys())
        contacts = []
        for i in range(len(R)-1):
            for j in range(i+1, len(R)):
                r1, r2 = R[i], R[j]
                dist = distance(r1, r2)
                if dist >= rmin and dist <= rmax:
                    index1, index2 = NCP_index[i], NCP_index[j]
                    ID1, ID2 = index_ID[index1], index_ID[index2]
                    contacts.append((ID1, ID2))
                    score1, score2 = ID_score[ID1], ID_score[ID2]
                    Ebind = bind_energy (score1, score2, metric=metric)
                    total -= Ebind
        if contact_info:
            return total, set(contacts)
        return total

    def make_contact_matrix (U, NCP_index):
        ete_list, R = go_through (U, NCP_index)
        N = len(R)
        matrix = Matrix(N, N)
        for i in range(N-1):
            for j in range(i+1, N):
                r1, r2 = R[i], R[j]
                dist = distance(r1, r2)
                if dist >= rmin and dist <= rmax:
                    matrix[i][j] = 1
                    matrix[j][i] = 1
        return matrix

    # set parameters for beads-on-rods model
    b = 0.1*p_len # one rod length
    k = float(p_len) / b # stiffness
    Pos_list = [ID_pos[ID] for ID in ID_score]
    L = 0.34*(max(Pos_list) - min(Pos_list) + 1) # full length
    N = int(math.ceil(float(L)/b)) + 1 # number of monomers
    
    # initialization
    IDs = sorted(ID_score.keys())
    ID_index, index_ID = ID_to_index (IDs, b)
    U = [[0.0, 0.0, 1.0]]*N
    E, Ct = energy(U, index_ID, model=model, metric=metric, contact_info=True) 
    
    # Beads-on-rods simulation
    print >> sys.stderr, "Beads-on-rods simulation"

    f = open(out_fname + '_bor_trace.txt', 'w')
    print >> f, "IDs" + ":" + ",".join(str(ID) for ID in IDs)
    print >> f, "CycleNum:" + str(cycle_num)
    print >> f, "Initial Contacts:"
    for ID1, ID2 in list(Ct):
        print >> f, "+:" + str(ID1) + ',' + str(ID2)

    g = open(out_fname + '_bor_energy.txt', 'w')
    print >> g, "IDs" + ":" + ",".join(str(ID) for ID in IDs)
    g.write(str(E))
    
    # MC cycle
    Trace = [U]
    E_list = [E]
    Ct_list = [Ct]
    i = 0
    order = -1
    while i < cycle_num:
        if int(math.log10(i+1)) > order:
            order = int(math.log10(i+1))
            print >> sys.stderr, "cycle num", i+1

        U = Trace.pop()
        E = E_list.pop()
        Ct = Ct_list.pop()

        new_U = [copy.deepcopy(U[0])]
        for j in range(1, len(U)):
            old_u = U[j]
            new_u = update(old_u, sigma)
            new_U.append(new_u)
        new_E, new_Ct = energy(new_U, index_ID, model=model, metric=metric, contact_info=True)

        if random.random() < min(1, math.exp(-new_E+E)):
            Trace.append(new_U)
            E_list.append(new_E)
            Ct_list.append(new_Ct)

            plus, minus = list(new_Ct - Ct), list(Ct - new_Ct)
            if len(plus + minus) > 0:
                print >> f, "\n@" + str(i+1)
                for ID1, ID2 in plus:
                    print >> f, "+:" + str(ID1) + ',' + str(ID2)
                for ID1, ID2 in minus:
                    print >> f, "-:" + str(ID1) + ',' + str(ID2)
            
            g.write(',' + str(new_E))
            i +=1
            continue

        Trace.append(U)
        E_list.append(E)
        Ct_list.append(Ct)
        g.write(',' + str(E))
        i +=1

    print >> sys.stderr, "Done"

    f.close()
    g.close()

    print >> sys.stderr, "\nCalculating mean contact matrix"

    ID_index = {}
    for i in range(len(IDs)):
        ID = IDs[i]
        ID_index[ID] = i

    mean_matrix = Matrix(len(IDs), len(IDs))
    D = Matrix(len(IDs), len(IDs))

    ctnum = 0
    ctnum_list = []

    pre_cycle = 1
    
    for line in open(out_fname + '_bor_trace.txt', 'r'):
        line = line.strip()
        if not line:
            continue
        if line.startswith('IDs'):
            continue
        if line.startswith('CycleNum'):
            continue
        if line.startswith('Initial'):
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

    f = open(out_fname + '_bor_mean.txt', 'w')
    print >> f, "IDs" + ":" + ",".join(str(ID) for ID in IDs)

    m, n = mean_matrix.getrank()
    for i in range(m):
        print >> f, ",".join(str(value) for value in mean_matrix[i])

    f.close()

    f = open(out_fname + '_bor_ctnum.txt', 'w')
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

    parser = ArgumentParser(description='Beads-on-rods simulation')
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
    parser.add_argument('--sigma',
                        dest="sigma",
                        type=int,
                        default=0.01,
                        help='update noise')
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
    parser.add_argument('--occ',
                        dest="occ_limit",
                        type=int,
                        default=4,
                        help='occupancy number limit')
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()

    if not args.burn_in:
        burn_in = int(args.cycle_num*0.1) # trash first 10% of total cycle
    else:
        burn_in = args.burn_in

    beads_on_rods (args.anot_fname,
                   args.field,
                   args.cycle_num,
                   burn_in,
                   args.sigma,
                   args.NCP_st,
                   args.NCP_num,
                   args.model,
                   args.metric,
                   args.rmin,
                   args.rmax,
                   args.p_len,
                   args.scale,
                   args.occ_limit,
                   args.out_fname)
