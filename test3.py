import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, sys, struct, os, operator
import matplotlib.pyplot as plt
import solver, data_processing, writer, measure

N = 3
L = 200
R = range(0,L)
g = 1.

def omega(j):
    return numpy.exp(math.pi*1.J*j/L)

def eps_per(k, L, N):
    z = omega(k) + numpy.power(omega(k), 2.*L-1)
    if z < 0.: return 0
    return numpy.power(omega(k) + numpy.power(omega(k), 2.*L-1), 2./N)

def exact_per(L, N):
    res = 0.
    for j in range(2*L):
        res -= eps_per(j, L, N)
    return res/L

print '____'
print solver.exact_eigenvalue(L, N, 1.)
print exact_per(L, N)
