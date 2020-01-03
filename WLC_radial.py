import os, sys, subprocess, re
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

# parameters
L = 100 # full contour length (nm) 
p_len = 50 # DNA persistence length (nm)

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

def biased_pick_on_sphere (u, sigma):
    acc_prob = -1
    while random.random() > acc_prob:
        new_u = pick_on_sphere ()
        angle = get_angle (u, new_u)
        acc_prob = 1.0/(np.sqrt(2.0*np.pi)*sigma)*np.exp(-0.5*angle*angle/(sigma*sigma))
    return new_u

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

def go_through (U, b):
    ete = [0.0, 0.0, 0.0]
    ete_list = []
    for i in range(len(U)):
        u = U[i]
        ete[0] += b*u[0]; ete[1] += b*u[1]; ete[2] += b*u[2]
        ete_list.append(copy.deepcopy(ete))
    return ete_list

def energy (U, k):
    total = 0.0
    for i in range(len(U)-1):
        u, u_next = U[i], U[i+1]
        angle = get_angle(u, u_next)
        total += 0.5*k*angle*angle
    return total

def plot (U, b):
    ete_list = go_through (U, b)
    P = [[0.0, 0.0, 0.0]] + ete_list[:-1]
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for p, u in zip(P, U):
        ax.quiver(p[0], p[1], p[2], u[0], u[1], u[2], length=b, color='k')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
    plt.close()
    return

# Worm like chain radial function (Becker, et al, EPJ, 2010) 
def WLC_radial (r, L=L, p_len=p_len):
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

# walking type WLC simulation
def WLC_walk (sample_num, L=L, p_len=p_len):
    
    # set parameters for beads-on-rods model
    b = 0.1*p_len # one rod length
    k = float(p_len) / b # stiffness
    N = int(math.ceil(float(L)/b)) + 1 # number of monomers

    print "persistence length:", p_len
    print "Full length:", L
    print "rod length:", b
    print "stiffness:", k
    print "monomer #:", N
    print "\nWLC walk simulation"
    
    Trace = []
    ete_dists = []
    i = 0
    order = -1
    while i < sample_num:
        if int(np.log10(i+1)) > order:
            order = int(np.log10(i+1))
            print "sample num", i+1

        U = [[0.0, 0.0, 1.0]]
        for j in range(N-1):
            prev_u = U[-1]
            new_u = biased_pick_on_sphere (prev_u, 1.0/np.sqrt(k))
            U.append (new_u)

        Trace.append(U)
        ete = go_through(U, b) [-1]
        ete_dists.append(norm(ete))
        i += 1
        
    assert len(Trace) == sample_num
    return Trace, ete_dists
#Trace, ete_dists = WLC_walk (1000)

# MC type WLC simulation
def WLC_MC (cycle_num, sigma, burn_in=None, L=L, p_len=p_len, fname='WLC_MC', erase=False, full_config=False):
    # set parameters for beads-on-rods model
    b = 0.1*p_len # one rod length
    k = float(p_len) / b # stiffness
    N = int(math.ceil(float(L)/b)) + 1 # number of monomers

    print "persistence length:", p_len
    print "Full length:", L
    print "rod length:", b
    print "stiffness:", k
    print "monomer #:", N
    print "\nWLC MC simulation"

    #initialization
    U = [[0.0, 0.0, 0.1]]*N
    #for i in range(N):
    #    U.append(pick_on_sphere())
    E = energy(U, k)
    eted = norm(go_through(U, b)[-1])

    if burn_in == None:
        burn_in  = int(0.1*cycle_num)

        
    if full_config:
        f = open(fname + "_trace.txt", "w")
        print >> f, "Monomer #:" + str(N)
        print >> f, "Cycle num:" + str(cycle_num)
        print >> f, "Burn in:" + str(burn_in)
        print >> f, "Initial config:"
        for l in range(len(U)):
            print >> f, str(l) + ":" + ",".join(str(value) for value in U[l])

    g = open(fname + "_energy.txt", "w")
    print >> g, "Monomer #:" + str(N)
    print >> g, "Cycle num:" + str(cycle_num)
    print >> g, "Burn in:" + str(burn_in)
    g.write(str(round(E,3)))

    h = open(fname + "_eted.txt", "w")
    print >> h, "Monomer #:" + str(N)
    print >> h, "Cycle num:" + str(cycle_num)
    print >> h, "Burn in:" + str(burn_in)
    h.write(str(round(eted,3)))
    
    Trace = [U]
    E_list = [E]
    eted_list = [eted]
    i = 0
    order = -1
    while i < cycle_num:
        if int(np.log10(i+1)) > order:
            order = int(np.log10(i+1))
            print "cycle num", i+1

        if full_config:
            print >> f, "\n@" + str(i+1)

        U = Trace.pop()
        E = E_list.pop()
        eted = eted_list.pop()
        
        new_U = [copy.deepcopy(U[0])]
        for j in range(1, len(U)):
            old_u = U[j]
            new_u = update(old_u, sigma)
            #new_u = biased_pick_on_sphere (old_u, sigma)
            new_U.append(new_u)
        new_E = energy(new_U, k)

        if random.random() < min(1, np.exp(-new_E+E)):
            Trace.append(new_U)
            E_list.append(new_E)
            new_eted = norm(go_through(new_U, b)[-1])
            eted_list.append(new_eted)
            if full_config:
                for l in range(len(new_U)):
                    print >> f, str(l) + ":" + ",".join(str(value) for value in new_U[l])
            g.write(',' + str(round(new_E,3)))
            h.write(',' + str(round(new_eted,3)))
            i +=1
            continue

        Trace.append(U)
        E_list.append(E)
        eted_list.append(eted)
        if full_config:
            for l in range(len(U)):
                print >> f, str(l) + ":" + ",".join(str(value) for value in U[l])
        g.write(',' + str(round(E,3)))
        h.write(',' + str(round(eted,3)))
        i +=1

    if full_config:
        f.close()

    g.close()
    h.close()

    if full_config:
        Trace = []
        for line in open(fname + "_trace.txt", 'r'):
            line = line.strip()
            if not line:
                continue
            if line.startswith("Monomer #:"):
                continue
            if line.startswith("Cycle num:"):
                continue
            if line.startswith("Burn in:"):
                continue
            if line.startswith("Initial config:"):
                temp = []
                continue
            if line.startswith("@"):
                cur_cycle = int(line[1:])
                Trace.append(temp)
                temp = []
                continue
            index, coordin = line.split(':')
            u = [float(value) for value in coordin.split(',')]
            temp.append(u)
        Trace.append(temp)
        assert len(Trace) == cycle_num + 1

    E_list = []
    for line in open(fname + "_energy.txt", "r"):
        line = line.strip()
        if not line:
            continue
        if line.startswith("Monomer #:"):
            continue
        if line.startswith("Cycle num:"):
            continue
        if line.startswith("Burn in:"):
            continue
        E_list = line.split(',')
    E_list = [float(value) for value in E_list]
    assert len(E_list) == cycle_num + 1

    eted_list = []
    for line in open(fname + "_eted.txt", "r"):
        line = line.strip()
        if not line:
            continue
        if line.startswith("Monomer #:"):
            continue
        if line.startswith("Cycle num:"):
            continue
        if line.startswith("Burn in:"):
            continue
        eted_list = line.split(',')
    eted_list = [float(value) for value in eted_list]
    assert len(eted_list) == cycle_num + 1

    if erase:
        subprocess.call(['rm', fname + '_trace.txt'])
        subprocess.call(['rm', fname + '_energy.txt'])
        subprocess.call(['rm', fname + '_eted.txt'])

    eted_list = eted_list[burn_in+1:]

    if full_config:
        return Trace, E_list, eted_list

    return E_list, eted_list
#E_list, ete_dists = WLC_MC (100000000, 0.01, erase=False)
#E_list, ete_dists = WLC_MC (10000, 0.01)

def read_temp_list (fname):
    for line in open(fname, "r"):
        line = line.strip()
        if not line:
            continue
        if line.startswith("Monomer #:"):
            continue
        if line.startswith("Cycle num:"):
            continue
        if line.startswith("Burn in:"):
            continue
        data_list = line.split(',')
    data_list = [float(value) for value in data_list]
    #assert len(eted_list) == cycle_num + 1
    return data_list
ete_dists = read_temp_list ("WLC_MC_eted.txt")
#energy_list = read_temp_list ("temp_energy.txt")

fig = plt.figure()
X = np.linspace(0.01, L-0.01)
Y = [ WLC_radial (x) for x in X]
plt.plot(X, Y, '--', label='analytic')
plt.hist(ete_dists, 50, density=True, label='simulation')
plt.title("Worm-like-Chain Polymer")
plt.xlabel("End to End distance (nm)")
plt.ylabel("Probability density")
plt.legend()
plt.show()
