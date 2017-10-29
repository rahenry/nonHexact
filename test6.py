import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, sys, struct, os, operator
import matplotlib.pyplot as plt
import solver, data_processing, writer, measure

def omega(j):
    return numpy.exp(2.*math.pi*1.J*j/N)

N = 3
L = 15
beta = (1-omega(1))**(N-1) / (1-omega(N-1))
beta = 1.
shift = 1
s = N*L
g = 1.
eps = 1.0
c0 = -(eps**(N/2))
gamma = 1.0
bc = 1

eps_list = []
for k in range(L):
    eps_list.append(solver.eps(k, L, N, gamma))

sol_list = [0]

for e in eps_list:
    new_list = []
    for s in sol_list:
        for i in range(N):
            x = s + omega(i) * e
            d = 1.0
            q = 1.
            #if not i==0: x -= d + q * 1.J * i
            new_list.append(x)
        #for i in range(N):
            #new_list.append(s + e * numpy.exp(1.9*math.pi*1.J / 2. / N * (i-1)))
    sol_list = new_list

print len(sol_list)

xdata = []
ydata = []
for e in sol_list:
    xdata.append(-e.real)
    ydata.append(-e.imag)
#print sol_list
plt.plot(xdata, ydata, linestyle='', marker='o', ms=3)
plt.axes().set_aspect('equal', 'datalim')
plt.show()

