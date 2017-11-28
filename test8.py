import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, sys, struct, os, operator
import matplotlib.pyplot as plt
import solver, data_processing, writer, measure

N = 3
L = 3
Q = N**L

columns = []

for total_charge in range(N):
    columns.append([])

for i in range(Q):
    for total_charge in range(N):
        x = 0
        for j in range(L):
            x += solver.sigma(i, j, N)
        x = x % N
        if total_charge == x:
            columns[total_charge].append(i)

X_mapping = []
for j in range(L):
    X_mapping.append([])
for i in range(N ** L):
    for j in range(L):
        solver.sigma_j = solver.sigma(i, j, N)
        if (solver.sigma_j < 1):
            X_mapping[j].append(solver.index(i, j, solver.sigma_j+1, N))
        elif (solver.sigma_j < N - 1):
            X_mapping[j].append(solver.index(i, j, solver.sigma_j+1, N))
        else:
            X_mapping[j].append(solver.index(i, j, 0, N))


def links(i):
    res = []
    for j in range(L):
        jf = (solver.sigma(i, j, N) + 1) % N
        res.append(solver.index(i, j, jf, N))
    return res

def matrix_element(x, y):
    if x in links(y): return 1
    return 0

import sympy

Z = 1

def matrix_element2(x, y, Z):
    x = columns[Z][x]
    y = columns[Z][y]
    if x in links(y): return 1
    return 0

V = len(columns[0])
M = sympy.Matrix(V, V, lambda i,j: 0)

for x in range(V):
    ls = links(columns[(Z-1)%N][x])
    for y in range(V):
        if columns[Z][y] in ls:
            M[x,y] = 1

rhs = sympy.Matrix(V, 1, lambda i,j: 0)


x, y, z, a, b, c = sympy.symbols('x, y, z, a, b, c')

print sympy.linsolve((M, rhs), sympy.symbols('a0:%i' % V))
        







        

