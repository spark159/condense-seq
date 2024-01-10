# random_loop By Sangwoo Park, 2014.
# generate random DNA chain loop to minimize some specified condition

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
    
# define rotation matrix: transform i+1 frame to i frame
def rot_matrix(tilt,roll,twist):
    a=tilt; b=roll; c=twist;
    
    #tilt matrix
    tilt_M=Matrix(3,3)
    tilt_M[0][0]=1
    tilt_M[1][1]=math.cos(a)
    tilt_M[1][2]=-math.sin(a)
    tilt_M[2][1]=math.sin(a)
    tilt_M[2][2]=math.cos(a)
    
    #roll matrix
    roll_M=Matrix(3,3)
    roll_M[0][0]=math.cos(b)
    roll_M[0][2]=math.sin(b)
    roll_M[1][1]=1
    roll_M[2][0]=-math.sin(b)
    roll_M[2][2]=math.cos(b)
    
    #twist matrix
    twist_M=Matrix(3,3)
    twist_M[0][0]=math.cos(c)
    twist_M[0][1]=-math.sin(c)
    twist_M[1][0]=math.sin(c)
    twist_M[1][1]=math.cos(c)
    twist_M[2][2]=1
    
    #order: tilt->roll->twist
    matrix=tilt_M*roll_M*twist_M

    # order: twist->tilt->roll
    #matrix=twist_M*tilt_M*roll_M

    return matrix

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
    return ete

# set dye coordinate in base-pair frame
def set_dye_bco (D,t0):
    D_bco=Matrix(len(D),3);
    for i in range(len(D)):
        if i < len(D)/2:
            D_bco[i][0]=r*math.cos((math.pi/5)*D[i]+t0); D_bco[i][1]=r*math.sin((math.pi/5)*D[i]+t0);
        else:
            D_bco[i][0]=r*math.cos((math.pi/5)*D[i]+t0+math.pi); D_bco[i][1]=r*math.sin((math.pi/5)*D[i]+t0+math.pi);    
    return D_bco

# get dye coordinate in reference frame
def get_dye_co (A):
    global a; global D; global D_bco;
    D_co=Matrix(len(D),3);
    for i in range(len(D)):
        if D[i] == 0:
            D_co[i]=D_bco[i] 
        
    for i in range(1,max(D)+1):
        R=rot_matrix(A[i-1][0],A[i-1][1],A[i-1][2])
        if i==1:
            product=R
            increase=Matrix(3,1); increase[2][0]=a;
        else:
            z_unit=Matrix(3,1); z_unit[2][0]=a;
            increase=increase+product*z_unit            
            product=product*R
        if i in D:
            idx=D.index(i)
            dye_loc=Matrix(3,1); dye_loc[0][0]=D_bco[idx][0]; dye_loc[1][0]=D_bco[idx][1]; dye_loc[2][0]=D_bco[idx][2];
            dye_co=product*dye_loc+increase
            D_co[idx][0]=dye_co[0][0]; D_co[idx][1]=dye_co[1][0]; D_co[idx][2]=dye_co[2][0]
    return D_co

# get calcuated FRET values
def get_E (A):
    global R0
    Cal_E=[]
    D_co=get_dye_co(A)
    n=D_co.getrank()[0]
    for i in range(0,n/2):
        for j in range(n/2,n):
            d=dist(D_co[i],D_co[j])
            Cal_E.append((1.0/(1.0+math.pow((d/R0),6))))
    return Cal_E

# get error between calculated FRET and experimental FRET
def get_error (A):
    global Exp_E
    Cal_E=get_E (A)
    error=0
    for i in range(len(Exp_E)):
        error+=math.pow((Cal_E[i]-Exp_E[i]),2)
    return error

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

#get dot-product of fake end-base-pair vector and z-axis vector in reference frame
def get_ebc (A):
    B=get_basis(A)
    ebsz=Matrix(1,3); ebsz[0]=B[B.getrank()[0]-1]
    zaxis=Matrix(3,1); zaxis[2][0]=1
    ebsy=Matrix(1,3); ebsy[0]=B[B.getrank()[0]-2]
    yaxis=Matrix(3,1); yaxis[1][0]=1    
    valuez=ebsz*zaxis; valuey=ebsy*yaxis
    return [valuez[0][0],valuey[0][0]]

# generate random loops from initial configuration
def random_loop (A, num, bound1, bound2, ebound):
    L=[]; n=0; N=A.getrank()[0]
    
    ete=get_ete(A); ebc=get_ebc(A)
    error=get_error(A)
    energy=get_energy(A)
    global error_list; global energy_list
    error_list.append(error); energy_list.append(energy);
    print str(n) +': ' 'Eted '+ str(dist(ete[0],[0,0,0])) +', Ebc ' + str(ebc) +', Error ' + str(error) +', Energy ' + str(energy) + '\n'
    L.append(A)
    
    while n<num:
        nA=Matrix(N,3);
        for i in range(N):
            for j in range(3):
                nA[i][j]=A[i][j]+random.uniform(-0.005,0.005)
        ete=get_ete(nA); ebc=get_ebc(nA); energy=get_energy(nA)
        if (dist(ete[0],[0,0,0]) < bound1) & (ebc[0] > bound2) & (ebc[1] > bound2) & (energy < ebound):
            n+=1
            error=get_error(nA)
            print str(n) +': ' 'Eted '+ str(dist(ete[0],[0,0,0])) +', Ebc ' + str(ebc) +', Error ' + str(error) +', Energy ' + str(energy) + '\n'            
            L.append(nA)
        else:
            continue
    return L

# draw DNA loop
def loop_plot (A):
    global a; global r0;
    N=A.getrank()[0]
    B=get_basis(A)
    
    X=[0,0]; Y=[0,0]; Z=[0,a]
    hX=[r0]; hY=[0]; hZ=[0]
    ihX=[r0*math.cos(4*math.pi/3)]; ihY=[r0*math.sin(4*math.pi/3)]; ihZ=[0]
    for i in range(1,N):
        R=rot_matrix(A[i-1][0],A[i-1][1],A[i-1][2])
        if i==1:
            product=R
        else:
            product=product*R
        haxis=Matrix(3,1); haxis[0][0]=r0*math.cos((math.pi/5)*i); haxis[1][0]=r0*math.sin((math.pi/5)*i);
        haxis=product*haxis; hX.append(X[i]+haxis[0][0]); hY.append(Y[i]+haxis[1][0]); hZ.append(Z[i]+haxis[2][0])
        ihaxis=Matrix(3,1); ihaxis[0][0]=r0*math.cos((math.pi/5)*i+4*math.pi/3); ihaxis[1][0]=r0*math.sin((math.pi/5)*i+4*math.pi/3);
        ihaxis=product*ihaxis; ihX.append(X[i]+ihaxis[0][0]); ihY.append(Y[i]+ihaxis[1][0]); ihZ.append(Z[i]+ihaxis[2][0])        
        zaxis=Matrix(3,1); zaxis[2][0]=a
        zaxis=product*zaxis; X.append(X[i]+zaxis[0][0]); Y.append(Y[i]+zaxis[1][0]); Z.append(Z[i]+zaxis[2][0])
    hX.append(hX[0]); hY.append(hY[0]); hZ.append(hZ[0]); 
    ihX.append(ihX[0]); ihY.append(ihY[0]); ihZ.append(ihZ[0]);
    
    D_co=get_dye_co(A)
    dX=[];dY=[];dZ=[]
    for i in range(D_co.getrank()[0]):
        dX.append(D_co[i][0]); dY.append(D_co[i][1]); dZ.append(D_co[i][2])
    
    # plot DNA loop
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')
    ax.plot(hX,hY,hZ)
    ax.plot(X, Y, Z)
    ax.plot(ihX,ihY,ihZ)
    ax.scatter(dX,dY,dZ,color='black')

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
    
# draw bar graph to compare Cal_E and Exp_E
def bar_plot(A):
    Cal_E = get_E(A)
    global Exp_E
    N=len(Exp_E)
    
    ind = np.arange(1,N+1)  # the x locations for the groups
    width = 0.35       # the width of the bars
    
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, Cal_E, width, color='b')
    rects2 = ax.bar(ind+width, Exp_E, width, color='r')
    
    ax.set_ylabel('FRET value')
    #ax.set_title('FRET value comparison')
    ax.set_xticks(ind+width)
    ax.set_xticklabels( ('1-1', '1-2', '1-3', '1-4', '2-1', '2-2', '2-3', '2-4','3-1', '3-2', '3-3', '3-4','4-1', '4-2', '4-3', '4-4') )
    
    ax.legend( (rects1[0], rects2[0]), ('Model', 'Experiment') )    

# do loop evolution through selection on random loops
def evolution (A, evo_num, ran_num, bound1, bound2, ebound):
    loop_plot(A); bar_plot(A)
    n=0
    while n < evo_num:
        L=random_loop(A, ran_num, bound1, bound2, ebound)
        eL=[]
        for i in range(len(L)):
            eL.append(get_error(L[i]))
        idx=eL.index(min(eL))
        A=L[idx]
        #loop_plot(A)
        n+=1
        print str(n) + '\n'
    loop_plot(A); bar_plot(A)
    ete=get_ete(A); ebc=get_ebc(A)
    error=get_error(A)
    energy=get_energy(A)
    global error_list; global energy_list
    error_list.append(error); energy_list.append(energy); 
    ee_plot(error_list, energy_list);  
    print 'Final' +': ' 'Eted '+ str(dist(ete[0],[0,0,0])) +', Ebc ' + str(ebc) +', Error ' + str(error) +', Energy ' + str(energy) + '\n'
    return A

def ee_plot(Er,En):
  fig, ax1 = plt.subplots()
  t = range(len(Er))  
  s1 = Er
  ax1.plot(t, s1, 'b-')
  ax1.set_xlabel('Iteration Num')
  ax1.set_ylabel('Error', color='b')
  for tl in ax1.get_yticklabels():
      tl.set_color('b')
      
  ax2 = ax1.twinx()
  s2 = En
  ax2.plot(t, s2, 'r-')
  ax2.set_ylabel('Energy', color='r')
  for tl in ax2.get_yticklabels():
      tl.set_color('r')
  plt.show()


# ______________physical parameters of DNA chain model
# number of base-pair in the loop
N=90;
# each base-pair step length/ dye-linker length/initial angle/radius of DNA/Foster distance
a=0.34; r=1.3; t0=1.7*math.pi; r0=1; R0=5;
# base-pair position(step, coordinate) of Dyes/ experimental FRET value
DT = [9, 21, 39, 48]; DB = [18, 63, 84, 0]; D=DT+DB;
D_bco=set_dye_bco(D,t0)
Exp_E=[0.79112, 0.19345, 0.49401, 0.77861, 0.87837, 0.13387, 0.15293, 0.31435, 0.21743, 0.31719, 0.18315, 0.28295, 0.1641, 0.68687, 0.22236, 0.35428]
# (tilt,roll,twist) angles in each dinucleotide steps  
A=Matrix(N,3);    
for i in range(N):
    A[i][1]=-math.pi+(N-2.0)*math.pi/N      # make perfect circle for initial condition
# (tilt,roll,twist) angles in natural configuration
A0=Matrix(N,3);
# (tilt,roll,twist) interaction coefficeints
K=Matrix(N,3);
for i in range(N):
    for j in range(3):
        K[i][j]=1

# error and energy list for generations
error_list=[]; energy_list=[];




# ______________main codes
#print get_ete(A)
#print ete, dist(ete[0],[0,0,0])
#print get_basis(A)
#print get_dye_co(A)
#Af=evolution(A,100,50,0.2)
#loop_plot(A); bar_plot(A)
#evolution(A,1000,20,0.15,0.9,1)  #input/generation/progeny/closure/continuity/energy
#evolution(A,1000,20,0.15,0.9,0.6)
evolution(A,500,20,0.15,0.9,5)
#loop_plot(A); bar_plot(A)
