import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, sys, struct, os, operator
import matplotlib.pyplot as plt
import solver, data_processing, writer, measure

def omega(j):
    return numpy.exp(2.*math.pi*1.J*j/N)

N = 3
<<<<<<< HEAD
L = 10
=======
L = 7
>>>>>>> 7a34989283f953bed4a0212702f22fb7bb197b22
beta = (1-omega(1))**(N-1) / (1-omega(N-1))
beta = 1.
beta = omega(1)
beta1 = omega(1)*0
beta2 = -omega(2)
shift = 1
s = N*L
g = 1.
eps = 1.0
c0 = -(eps**(N/2))
gamma = 1.0
bc = 1

def m_element_P(i, j, s, bc):
    d = (j - i)
    if d == -1: return 1.
    if bc == 1: 
        d = d % s
    if d == s-1:
        if j == s-1: return beta1
    if d == N-1: 
        z = gamma
        if j == 0: z *= beta2
        if j % N == 0: return z
        if j % N == N-1: return z
        #if j % N == 1 or i % N == 1: return gamma
        return 0
    return 0

R = range(0,s)
M_P = numpy.zeros((s,s), dtype=numpy.complex_)
#M_P = numpy.zeros((s,s))
for i in R:
    for j in R:
        M_P[i, j] = m_element_P(i, j, s, bc)
print M_P
eigs = numpy.linalg.eigvals(M_P)
res = 0.
print eigs
counter = 0
for e in eigs:
    if abs(e.imag) < 1E-5:
        counter +=1
        print e
print "total = " + str(counter)


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

res = 0.
eps_list = []
for e in eigs:
    if e.real > 0:
        eps_list.append(e)
        res += e

print "e = " +  str(res / L)

xdata = []
ydata = []
for e in eigs:
    xdata.append(e.real)
    ydata.append(e.imag)
#plt.plot(xdata, ydata, linestyle='', marker='o')
#plt.show()

sol_list = [(0,0)]

for e in eps_list:
    new_list = []
    for s in sol_list:
        for i in range(N):
            charge = s[1] + i
            x = s[0] + omega(i) * e
            d = 1.0
            q = 1.
            #if not i==0: x -= d + q * 1.J * i
            new_list.append((x, charge))
        #for i in range(N):
            #new_list.append(s + e * numpy.exp(1.9*math.pi*1.J / 2. / N * (i-1)))
    sol_list = new_list

new_list = []
Q = -0.333
for s in sol_list:
    charge = s[1] % N
    charge = s[1]
    x = s[0] + Q*charge
    new_list.append(x)

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

