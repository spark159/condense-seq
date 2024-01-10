import math
import random

def comb(n, b):
    return float(math.factorial(n)) / (math.factorial(n-b)*math.factorial(b))

def w(n,f,l,g):
    return (f**n)*(g**l)*math.factorial((f*n-n+g*l-l)) / float(math.factorial((f*n-2*n+g*l-2*l+2)))


def re_w(n,f,g):
    output = 0.0
    for b in range(1, min(n+1, g+1)):
        output += comb(n, b)*comb(g,b)*math.factorial(b)*(f**b)*w(n-b,f,1,b*(f-1))
    return output

N = 100
Ps = []
fs = [10.0]*N
for i in range(N):
    Ps.append(0.5 + (0.1) * random.gauss(0,1))
    #Ps.append(0.5)
Pmean = sum(Ps)/float(len(Ps))
fmean = sum(fs)/float(len(fs))

product = 1.0
for i in range(len(Ps)):
    product = product * ((Ps[i]**(1.0/fs[i]))*(1-(Pmean**(1.0/fmean)))) / ((Pmean**(1.0/fmean))*(1-(Ps[i]**(1.0/fs[i]))))

print product



