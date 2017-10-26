import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, sys, struct, os, operator
import matplotlib.pyplot as plt
import solver, data_processing, writer, measure

N = 3
L = 12
R = range(0,L)
g = 1.

def m_element_P(i, j):
    d = (j - i) % L
    if d == N-1: return g
    if d == (-1 % L): return 1.
    return 0

def m_element_O(i, j):
    d = (j - i)
    if d == N-1: return g
    if d == -1: return 1.
    return 0

M_P = numpy.zeros((L,L))
for i in R:
    for j in R:
        M_P[i][j] = m_element_P(i, j)
print M_P
M_O = numpy.zeros((L,L))
for i in R:
    for j in R:
        M_O[i][j] = m_element_O(i, j)
print M_O

print numpy.linalg.eigvals(M_P)
print numpy.linalg.eigvals(M_O)
