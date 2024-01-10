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
from mpl_toolkits.mplot3d import Axes3D

# load file
path = "/home/spark159/../../media/spark159/sw/dataforcondense//hg19_chr1_171_everything_anot.cn"
ID_chr, ID_pos, name_ID_value = load_file.read_anot_file(path)
ID_score = name_ID_value['work/condense_seq/sp10_hg19_chr1']
temp = {}
for i in range(0, 10):
    temp[i] = ID_score[i]
ID_score = temp

# parameters
rmin, rmax = 0, 50 # minimum and maximum reaction distance (nm)
p_len = 50 # DNA persistence length (nm)
#scale = 0.413 # total count ratio of test sample to control (UV data, sp9)
scale = 0.154 # total count ratio of test sample to control (UV data, sp10)
occ_limit = 4 # occupancy number limit

# set parameters for beads-on-rods model
b = 0.1*p_len # one rod length
k = float(p_len) / b # stiffness
Pos_list = [ID_pos[ID] for ID in ID_score]
L = 0.34*(max(Pos_list) - min(Pos_list) + 1) # full length
N = int(math.ceil(float(L)/b)) + 1 # number of monomers

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
    return np.sqrt(total)

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
    return np.sqrt(total)

def get_angle (v1, v2):
    dot_product = 0.0
    for i in range(3):
        dot_product += v1[i]*v2[i]
    cos = dot_product/(norm(v1)*norm(v2))
    return np.arccos(cos)

def polar_to_cartesian (theta, phi):
    x = np.sin(phi)*np.cos(theta)
    y = np.sin(phi)*np.sin(theta)
    z = np.cos(phi)
    return [x, y, z]

def cartesian_to_polar (x, y, z):
    rho = np.sqrt(x**2 + y**2 + z**2)
    r = np.sqrt(x**2 + y**2)
    theta = np.arccos(float(z)/rho)
    phi = np.arccos(float(x)/r)
    return theta, phi

def rotation (u, angle, p):
    def h_product(q1, q2):
        x = q1[0]*q2[0]-q1[1]*q2[1]-q1[2]*q2[2]-q1[3]*q2[3]
        y = q1[0]*q2[1]+q1[1]*q2[0]+q1[2]*q2[3]-q1[3]*q2[2]
        z = q1[0]*q2[2]-q1[1]*q2[3]+q1[2]*q2[0]+q1[3]*q2[1]
        w = q1[0]*q2[3]+q1[1]*q2[2]-q1[2]*q2[1]+q1[3]*q2[0]
        return [x, y, z, w]
    q = [np.cos(float(angle)/2), u[0]*np.sin(float(angle)/2), u[1]*np.sin(float(angle)/2), u[2]*np.sin(float(angle)/2)]
    q_inv = [np.cos(float(angle)/2), -u[0]*np.sin(float(angle)/2), -u[1]*np.sin(float(angle)/2), -u[2]*np.sin(float(angle)/2)]
    p = [0, p[0], p[1], p[2]]
    new_p = h_product(h_product(q, p), q_inv) 
    return new_p[1:]

def pick_on_sphere ():
    u, v = random.uniform(-1, 1), random.random()
    x = np.sqrt(1-u**2)*np.cos(2*np.pi*v)
    y = np.sqrt(1-u**2)*np.sin(2*np.pi*v)
    z = u
    return [x, y, z]

def update (old_u, sigma):
    theta = random.gauss(0, sigma)
    dum = copy.deepcopy(old_u)
    dum[0] += 1
    r = cross(dum, old_u)
    r[0] *= np.tan(theta)/norm(r); r[1] *= np.tan(theta)/norm(r); r[2] *= np.tan(theta)/norm(r)
    v = [old_u[0]+r[0], old_u[1]+r[1], old_u[2]+r[2]]
    phi = random.uniform(0, 2*np.pi)
    new_v = rotation(old_u, phi, v)
    new_u = new_v/norm(new_v)
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
    prob1 = 1.0 - scale*np.exp(-score1)
    assert prob1 >=0 and prob1 <= 1
    prob2 = 1.0 - scale*np.exp(-score2)
    assert prob2 >=0 and prob2 <= 1
    if metric == "product":
        prob = prob1*prob2
    elif metric == "GM":
        prob = np.sqrt(prob1*prob2)
    Ebind = np.log(prob/(1-prob)) # unit of KT
    return Ebind

def energy (U, index_ID, model, metric):
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
    for i in range(len(R)-1):
        for j in range(i+1, len(R)):
            r1, r2 = R[i], R[j]
            dist = distance(r1, r2)
            if dist >= rmin and dist <= rmax:
                index1, index2 = NCP_index[i], NCP_index[j]
                ID1, ID2 = index_ID[index1], index_ID[index2]
                score1, score2 = ID_score[ID1], ID_score[ID2]
                Ebind = bind_energy (score1, score2, metric=metric)
                total -= Ebind
    return total

def plot (U, NCP_index):
    ete_list, R = go_through (U, NCP_index)
    P = [[0.0, 0.0, 0.0]] + ete_list[:-1]
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for p, u in zip(P, U):
        ax.quiver(p[0], p[1], p[2], u[0], u[1], u[2], length=b, color='k')
    X, Y, Z = [], [], []
    for r in R:
        X.append(r[0])
        Y.append(r[1])
        Z.append(r[2])
    ax.scatter(X, Y, Z, color='b', s=50)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
    plt.close()
    return

def make_contact_matrix (U, NCP_index):
    ete_list, R = go_through (U, NCP_index)
    N = len(R)
    matrix = np.zeros((N, N))
    for i in range(N-1):
        for j in range(i+1, N):
            r1, r2 = R[i], R[j]
            dist = distance(r1, r2)
            if dist >= rmin and dist <= rmax:
                matrix[i][j] = 1
                matrix[j][i] = 1
    return matrix

def beads_on_rods_simulation (ID_score, cycle_num, sigma, model, metric, graph=False):
    # initialization
    #U = []
    #for i in range(N):
    #    U.append(pick_on_sphere())
    U = [[0.0, 0.0, 1.0]]*N
    ID_index, index_ID = ID_to_index (ID_score.keys(), b)

    if graph:
        plot(U, index_ID.keys())

    # MC cycle
    Trace = [U]
    E_list = [energy(U, index_ID, model=model, metric=metric)]
    k = 0
    order = -1
    while k < cycle_num:
        if int(np.log10(k+1)) > order:
            order = int(np.log10(k+1))
            print "cycle num", k+1
        U = copy.deepcopy(Trace[-1])
        E = copy.deepcopy(E_list[-1])
        new_U = [copy.deepcopy(U[0])]
        for i in range(1, len(U)):
            old_u = U[i]
            new_u = update(old_u, sigma)
            new_U.append(new_u)
        new_E = energy(new_U, index_ID, model=model, metric=metric)
        if random.random() < min(1, np.exp(-new_E+E)):
            Trace.append(new_U)
            E_list.append(new_E)
            k +=1
            continue
        Trace.append(U)
        E_list.append(E)
        k +=1
    print "Done"

    if graph:
        plot(U, index_ID.keys())

    return Trace, E_list
Trace, E_list = beads_on_rods_simulation (ID_score, cycle_num=100000, sigma=0.01, model='WLC', metric='product', graph=True)

fig = plt.figure()
plt.plot(range(len(E_list)), E_list)
plt.xlabel("Steps")
plt.ylabel("Energy (kT)")
plt.show()
plt.close()

ID_index, index_ID = ID_to_index (ID_score.keys(), b)

def Video (Trace, NCP_index, step):    
    fig = plt.figure()
    for i in range(0, len(Trace), step):
        matrix = make_contact_matrix(Trace[i], NCP_index)
        #matrix = M_list[i]
        plt.spy(matrix)
        plt.draw()
        plt.pause(0.1)
    plt.close()
    return None
Video (Trace, index_ID.keys(), step=1000)

def mean_contact_matrix (Trace, NCP_index, burn_in=None):
    if burn_in == None:
        burn_in = int(len(Trace)*0.1)
    for i in range(burn_in, len(Trace)):
        matrix = make_contact_matrix(Trace[i], NCP_index)
        if i == burn_in:
            mean = matrix
            continue
        mean += matrix
    return mean/(len(Trace)-burn_in)
mean_matrix = mean_contact_matrix(Trace, index_ID.keys())

fig = plt.figure()
plt.imshow(mean_matrix, cmap='Reds')
plt.colorbar()
plt.title("Mean Contact Probability")
plt.show()
plt.close()
