import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, sys, struct, os, operator
import matplotlib.pyplot as plt
import solver, data_processing, writer, measure

N = 3
L = 8
R = range(0,L)
g = 1.
eps = 0.1
c0 = -(eps**(N/2))

def m_element_P(i, j):
    d = (j - i) % L
    if d == 0: return -(eps**(N/2))
    if d == (-1 % L): return 1.
    if d == 1: return 1.
    return 0

M_P = numpy.zeros((L,L))
for i in R:
    for j in R:
        M_P[i][j] = m_element_P(i, j)
print M_P
print numpy.linalg.eigvals(M_P)

def omega(j):
    return numpy.exp(2.*math.pi*1.J*j/L)

def eigen(j):
    res = 0.
    return c0 + omega(j) + omega(j) ** (L-1)

for i in R:
    print eigen(i).real

def eps_per(k, L, N):
    return pow(2. * numpy.cos(math.pi * k / 2./ L), 2./N)
def exact_per(L, N):
    res = 0.
    for j in range(L):
        res -= eps_per(j+1, L, N)
    return res/L

print '____'
print solver.exact_eigenvalue(L, N, 1.)
print exact_per(L, N)
