import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, sys, struct, os, operator
import matplotlib.pyplot as plt
import solver, data_processing, writer, measure

N = 3
L = 12
s = 3*L
g = 1.
eps = 0.1
c0 = -(eps**(N/2))
gamma = 1.
bc = 0

def m_element_P(i, j, s, bc):
    d = (j - i)
    if bc == 1: d = d % s
    if bc == 0: 
        if d == -1: return 1.
    else: 
        if d == s-1: return 1.
    if d == 2: return gamma
    return 0

R = range(0,s)
M_P = numpy.zeros((s,s))
for i in R:
    for j in R:
        M_P[i][j] = m_element_P(i, j, s, bc)
print M_P
eigs = numpy.linalg.eigvals(M_P)
res = 0.
print eigs
for e in eigs:
    if abs(e.imag) < 1E-8:
        print e

def omega(j):
    return numpy.exp(2.*math.pi*1.J*j/L)

def eigen(j):
    res = 0.
    return c0 + omega(j) + omega(j) ** (L-1)

#for i in R:
    #print eigen(i).real

def eps_per(k, L, N):
    return pow(2. * numpy.cos(math.pi * k / 2./ L), 2./N)
def exact_per(L, N):
    res = 0.
    for j in range(L):
        res -= eps_per(j+1, L, N)
    return res/L

#print '____'
print solver.exact_eig2(L, N, 1.)
#print exact_per(L, N)
xdata = []
ydata = []
for e in eigs:
    xdata.append(e.real)
    ydata.append(e.imag)
plt.plot(xdata, ydata, linestyle='', marker='o')
#plt.show()

res = 0.
for e in eigs:
    if e.real > 0:
        res += e**(1./N)

print res / L
