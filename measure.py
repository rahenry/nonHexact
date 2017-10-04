import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, sys, struct, os, operator
import matplotlib.pyplot as plt
import solver, data_processing, writer

def current_1(i, sol, k):
    res = 0.
    L = sol['L']
    for j in range(1, sol['L']):
        s1 = solver.sigma(i, j, sol['N'])
        s2 = solver.sigma(i, j+1, sol['N'])
        #res += numpy.sin(2*math.pi/sol['N'] * (s2-s1))
        #res += numpy.exp(u * k * numpy.exp(2.*math.pi*(s2-s1)/sol['N']))
    if sol['bc'] == 1:
        s1 = solver.sigma(i, sol['L'], sol['N'])
        s2 = solver.sigma(i, 1, sol['N'])
        #res += numpy.sin(2*math.pi/sol['N'] * (s2-s1))
        #res += numpy.exp(u * k * numpy.exp(2.*math.pi*(s2-s1)/sol['N']))
    s1 = solver.sigma(i, 0, sol['N'])
    s2 = solver.sigma(i, L-1, sol['N'])
    res += numpy.sin(2.*math.pi/sol['N'] * (s2-s1))
    div = sol['L']
    #if sol['bc'] == 1:
        #div += 1
    res = res.real
    return res / div
def momentum_1(i, sol, k):
    res = 0. + 0.J
    u = 1.J
    L = sol['L']
    for j in range(1, sol['L']):
        s1 = solver.sigma(i, j, sol['N'])
        s2 = solver.sigma(i, j+1, sol['N'])
        res += numpy.exp(u * k * numpy.exp(2.*math.pi*(s2-s1)/sol['N']))
    if sol['bc'] == 1:
        s1 = solver.sigma(i, 1, sol['N'])
        s2 = solver.sigma(i, sol['L'], sol['N'])
        #res += numpy.exp(u * k * numpy.exp(2.*math.pi*(s2-s1)/sol['N']))
    div = sol['L']
    #if sol['bc'] == 1:
        #div += 1
    return res / div

def momentum_3(i, sol, k):
    res = 0. + 0.J
    u = 1.J
    L = sol['L']
    for j in range(sol['L']):
        s1 = solver.sigma(i, j, sol['N'])
        res += numpy.exp(u * k * j) * numpy.exp(2.J*math.pi*(s1)/sol['N'])
    if sol['bc'] == 1:
        s1 = solver.sigma(i, 1, sol['N'])
        s2 = solver.sigma(i, sol['L'], sol['N'])
        #res += numpy.exp(u * k * numpy.exp(2.*math.pi*(s2-s1)/sol['N']))
    div = sol['L']
    #if sol['bc'] == 1:
        #div += 1
    return res / div

def momentum_2(i, sol, k):
    res = 0. + 0.J
    u = 1.J
    L = sol['L']
    for j in range(2, sol['L']):
        s1 = solver.sigma(i, j-1, sol['N'])
        s2 = solver.sigma(i, j, sol['N'])
        s3 = solver.sigma(i, j+1, sol['N'])
        d1 = numpy.exp(u * k * numpy.exp(2.*math.pi*(s2-s1)/sol['N']))
        d2 = numpy.exp(u * k * numpy.exp(2.*math.pi*(s3-s2)/sol['N']))
        res += 0.5 * (d1+d2)
    if sol['bc'] == 1:
        s1 = solver.sigma(i, L, sol['N'])
        s2 = solver.sigma(i, 1, sol['N'])
        s3 = solver.sigma(i, 2, sol['N'])
        d1 = numpy.exp(u * k * numpy.exp(2.*math.pi*(s2-s1)/sol['N']))
        d2 = numpy.exp(u * k * numpy.exp(2.*math.pi*(s3-s2)/sol['N']))
        res += 0.5 * (d1+d2)
        s1 = solver.sigma(i, L-1, sol['N'])
        s2 = solver.sigma(i, L, sol['N'])
        s3 = solver.sigma(i, 1, sol['N'])
        d1 = numpy.exp(u * k * numpy.exp(2.*math.pi*(s2-s1)/sol['N']))
        d2 = numpy.exp(u * k * numpy.exp(2.*math.pi*(s3-s2)/sol['N']))
        res += 0.5 * (d1+d2)
    else:
        s1 = solver.sigma(i, 1, sol['N'])
        s2 = solver.sigma(i, 2, sol['N'])
        d1 = numpy.exp(u * k * numpy.exp(2.*math.pi*(s2-s1)/sol['N']))
        res += d1
        s1 = solver.sigma(i, L-1, sol['N'])
        s2 = solver.sigma(i, L, sol['N'])
        d1 = numpy.exp(u * k * numpy.exp(2.*math.pi*(s2-s1)/sol['N']))
        res += d1

    div = sol['L']
    #if sol['bc'] == 1:
        #div += 1
    return res / div



def measure_scalar(sol, f, *args):

    res = 0. + 0.J
    sol['spin_expectations'] = []
    for j in range(0,sol['L']):
        se = 0. + 0.J
        for i, e in enumerate(sol['evec']):
            s1 = solver.sigma(i, j, sol['N'])
            angle = 2.*math.pi/sol['N']*s1
            se += numpy.exp(1.J*angle)*e*numpy.conj(e)
            #se += angle*e*numpy.conj(e)

        #se /= len(sol['evec'])
        se = se.real
        sol['spin_expectations'].append(se)

    print sol['spin_expectations']
    for i, e in enumerate(sol['evec']):
        res += numpy.conj(e) * f(i, sol, *args) * e


    return res
