# MC_looping, By Sangwoo Park, 2016
# MC simulation to calculate J factor, 2016
# Used CEHS system, ring-closure boundary condition (R, gamma, tau)

import math
import random
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

# ______________define functions
# define matrix class
class Matrix (object):
    def __init__ (self, m,n):
        self.m=m; self.n=n;
        self.matrix=[]
        for i in range(m):
            self.matrix.append([0]*n)
    
    def __getitem__ (self, idx):
        return self.matrix[idx]
    
    def __setitem__ (self, idx, value):
        self.matrix[idx]=value
    
    def getrank (self):
        return (self.m,self.n)
        
    def __add__ (self, other):
        if self.getrank() != other.getrank():
            print "Not possible"
        else:
            new_one=Matrix(self.m,self.n)
            for i in range(self.m):
                for j in range(self.n):
                    new_one.matrix[i][j]=self.matrix[i][j]+other.matrix[i][j]
        return new_one
    
    def __sub__ (self, other):
        if self.getrank() != other.getrank():
            print "Not possible"
        else:
            new_one=Matrix(self.m,self.n)
            for i in range(self.m):
                for j in range(self.n):
                    new_one.matrix[i][j]=self.matrix[i][j]-other.matrix[i][j]
        return new_one    

    def __mul__ (self, other):
        if self.n != other.m:
            print "Not possible"
        else:
            new_one=Matrix(self.m, other.n)
            for i in range(self.m):
                for j in range(other.n):
                    value=0
                    for k in range(self.n):
                        value +=self.matrix[i][k]*other.matrix[k][j]
                    new_one.matrix[i][j]=value
        return new_one
    
    def __str__ (self):
        s=''
        for i in range(self.m):
            s += str(self.matrix[i]) + '\n'
        return s[0:len(s)-1]

def unit_rot (direct, a):
    A=Matrix(3,3)
    
    if direct=='x':
        A[0][0]=1.0
        A[1][1]=math.cos(a)
        A[1][2]=-math.sin(a)
        A[2][1]=math.sin(a)
        A[2][2]=math.cos(a)
    
    elif direct=='y':
        A[0][0]=math.cos(a)
        A[0][2]=math.sin(a)
        A[1][1]=1.0
        A[2][0]=-math.sin(a)
        A[2][2]=math.cos(a)
    
    elif direct=='z':
        A[0][0]=math.cos(a)
        A[0][1]=-math.sin(a)
        A[1][0]=math.sin(a)
        A[1][1]=math.cos(a)
        A[2][2]=1.0

    return A

# TiltRoll approximation, change (Tilt, Roll) to (TiltRoll, phi)
def TiltRoll_app (Tilt, Roll):
    TiltRoll=math.sqrt(Tilt**2 + Roll**2)
    phi=math.atan2(float(Tilt),float(Roll))
    return TiltRoll, phi

# inverse of TiltRoll approximation
def inv_TiltRoll (TiltRoll, phi):
    Tilt=TiltRoll*math.sin(phi)
    Roll=TiltRoll*math.cos(phi)
    return Tilt, Roll

# define rotation matrix: transform i+1 frame to i frame
def rot_matrix (Tilt, Roll, Twist):
    TR, phi = TiltRoll_app (Tilt, Roll)
    return unit_rot('z', Twist*0.5-phi) * unit_rot('y', TR) * unit_rot('z', Twist*0.5+phi)

# get dinstance between two list
def dist (V,W):
    return math.sqrt(math.pow(V[0]-W[0],2) + math.pow(V[1]-W[1],2)  + math.pow(V[2]-W[2],2) )

# get Basis coordinate in reference frame
def get_basis (A):
    N=A.getrank()[0]+1
    B=Matrix(3*N,3)
    B[0]=[1,0,0]; B[1]=[0,1,0]; B[2]=[0,0,1];
    
    for i in range(1,N):
        R=rot_matrix(A[i-1][0],A[i-1][1],A[i-1][2])
        if i==1:
            product=R
        else:
            product=product*R
        for k in range(3):
            B[3*(i)+k][0]=product[0][k];B[3*(i)+k][1]=product[1][k];B[3*(i)+k][2]=product[2][k]
    return B

# get end-to-end distance vector
def get_ete (A):
    N=A.getrank()[0]
    B=get_basis(A)
    ete=Matrix(1,3)
    for i in range(N):
        add=Matrix(1,3); add[0]=B[3*(i+1)-1]
        ete=ete+add
    return dist(ete[0],[0,0,0])

#get dot-product of end-base-pair unit vector and z-axis unit vector in reference frame
def get_gamma (A):
    B=get_basis(A)
    return B[B.getrank()[0]-1][2]

#get torsional alignment of end-basepair and z-axis in reference frame
def get_tau (A):
    B=get_basis(A)
    return B[B.getrank()[0]-3][0]

# get configurational energy of DNA
def get_energy (A):
    global A0; global K
    energy=0
    diff=A-A0
    [m,n]=diff.getrank()
    for i in range(m):
        for j in range(n):
            energy += K[i][j]*math.pow(diff[i][j],2)
    return energy

# generate a random DNA chain
def random_chain ():
    global A0; global K;
    N=A0.getrank()[0]
    A=Matrix(N,3)
    for i in range(N):
        for j in range(3):
            A[i][j]=random.gauss(A0[i][j],K[i][j])
    return A

# collect data through selection on random chains
def Looping (num, bound1, bound2, bound3):
    hit=0; ete_list=[]; gam_list=[]; tau_list=[];
    for i in range(num):
        A=random_chain()
        ete=get_ete(A); gam=get_gamma(A); tau=get_tau(A)
        #DNA_plot(A)
        ete_list.append(ete)
        if ete < bound1:
            gam_list.append(gam)
            if gam > bound2:
                tau_list.append(tau)
                if tau > bound3:
                    hit +=1; DNA_plot(A)
    return hit/float(num), ete_list, gam_list, tau_list

# draw DNA chian
def DNA_plot (A):
    global a; global r0;
    N=A.getrank()[0]
    B=get_basis(A)
    
    X=[0,0]; Y=[0,0]; Z=[0,a]
    haxis0=Matrix(3,1); haxis0[0][0]=r0*math.cos(0); haxis0[1][0]=r0*math.sin(0);
    ihaxis0=Matrix(3,1); ihaxis0[0][0]=r0*math.cos(0+ 4*math.pi/3); ihaxis0[1][0]=r0*math.sin(0+ 4*math.pi/3);
    hX=[haxis0[0][0]]; hY=[haxis0[1][0]]; hZ=[haxis0[2][0]]
    ihX=[ihaxis0[0][0]]; ihY=[ihaxis0[1][0]]; ihZ=[ihaxis0[2][0]]
    for i in range(1,N+1):
        R=rot_matrix(A[i-1][0],A[i-1][1],A[i-1][2])
        if i==1:
            product=R
        else:
            product=product*R
        haxis=product*haxis0; hX.append(X[i]+haxis[0][0]); hY.append(Y[i]+haxis[1][0]); hZ.append(Z[i]+haxis[2][0])
        ihaxis=product*ihaxis0; ihX.append(X[i]+ihaxis[0][0]); ihY.append(Y[i]+ihaxis[1][0]); ihZ.append(Z[i]+ihaxis[2][0])        
        zaxis=Matrix(3,1); zaxis[2][0]=a
        zaxis=product*zaxis; X.append(X[i]+zaxis[0][0]); Y.append(Y[i]+zaxis[1][0]); Z.append(Z[i]+zaxis[2][0])
    #hX.append(hX[0]); hY.append(hY[0]); hZ.append(hZ[0]); 
    #ihX.append(ihX[0]); ihY.append(ihY[0]); ihZ.append(ihZ[0]);
    
    # plot DNA loop
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')
    ax.plot(hX,hY,hZ)
    ax.plot(X, Y, Z)
    ax.plot(ihX,ihY,ihZ)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')    
    ax.set_zlabel('Z')
    
    # make fake box to equalize axis units
    # Create cubic bounding box to simulate equal aspect ratio
    max_range = np.array([max(hX)-min(hX), max(hY)-min(hY), max(hZ)-min(hZ)]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(max(hX)+min(hX))
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(max(hY)+min(hY))
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(max(hZ)+min(hZ))
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')

    plt.grid()
    plt.show()
    
def get_radial_dist(ete_list,n):
    L=N*a; d=L/n; Total=float(len(ete_list))
    prob_list=Matrix(1,n)
    for ete in ete_list:
        ete=ete*a
        prob_list[0][int(ete//d)] +=1
    prob_den=[]; x_line=[]
    for i in range(n):
        if prob_list[0][i] !=0:
            prob_den.append(math.log10((prob_list[0][i]/Total)/((4*math.pi/3)*((d*(i+1))**3-(d*i)**3))))
            x_line.append(d*(i+0.5))
    return x_line, prob_den

def get_orient_dist(gam_list,n):
    d=math.pi/n; Total=len(gam_list)
    prob_list=Matrix(1,n)
    for gam in gam_list:
        ang=math.acos(gam)
        prob_list[0][int(ang//d)] +=1
    prob_den=[]; x_line=[]
    for i in range(n):
        if prob_list[0][i] !=0:
            prob_den.append(math.log10((prob_list[0][i]/Total)/d))
            x_line.append(d*(i+0.5))            
    return x_line, prob_den

def graph_fit (x,y):
    z=np.polyfit(x, y, 3)
    f=np.poly1d(z)
    x_new=np.linspace(x[0],x[-1],50)
    y_new=f(x_new)
    plt.figure()
    plt.plot(x,y,'o',x_new,y_new)
    return f(0)

def get_J (iter_num, bound1, bound2=-2, bound3=-2):
    [prob, ete_list, gam_list, tau_list]=Looping(iter_num, bound1, bound2, bound3)
    [x_line1, prob_den1]=get_radial_dist (ete_list, 20)
    radial0=10**graph_fit(x_line1,prob_den1)
    J= radial0/(6.02*(10**23))
    if bound2 !=-2:
        [x_line2, prob_den2]=get_orient_dist (gam_list, 20)
        gamma0=10**graph_fit(x_line2,prob_den2)
        J= (2*radial0*gamma0)/(6.02*(10**23))
    return J

# ______________physical parameters of DNA chain model
# number of base-pair in the loop
N=100;
# each base-pair step length/radius of DNA/
a=0.34; r0=1;
# (tilt,roll,twist) angles in natural configuration
A0=Matrix(N-1,3);
for i in range(N-1):
    A0[i][2]=math.pi/5
# (tilt,roll,twist) interaction coefficeints
K=Matrix(N-1,3);
for i in range(N-1):
    #K[i][0]=0.084; K[i][1]=0.084; K[i][2]=0.064
    K[i][0]=0.15; K[i][1]=0.15; K[i][2]=0.064
    
# error and energy list for generations
error_list=[]; energy_list=[];

# ______________main codes
#[prob, ete_list, gam_list, tau_list]=Looping(100,6,0.9,-2)
#[prob, ete_list, gam_list, tau_list]=Looping(1000,6,-2,-2)
#[x_line, prob_den]=get_radial_dist(ete_list,20)
#print graph_fit(x_line, prob_den)
#[x_line, prob_den]=get_orient_dist(gam_list,20)
print get_J (1000,6)