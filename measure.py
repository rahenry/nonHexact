import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, sys, struct, os, operator
import matplotlib.pyplot as plt
import solver, data_processing, writer

def momentum_1(i, sol, k):
    res = 0. + 0.J
    u = 1.J
    for j in range(1, sol['L']):
        s1 = solver.sigma(i, j, sol['N'])
        s2 = solver.sigma(i, j+1, sol['N'])
        res += numpy.exp(u * k * numpy.exp(2.*math.pi*(s2-s1)/sol['N']))
    if sol['bc'] == 1:
        s1 = solver.sigma(i, 1, sol['N'])
        s2 = solver.sigma(i, sol['L'], sol['N'])
        #res += numpy.exp(u * k * numpy.exp(2.*math.pi*(s2-s1)/sol['N']))
    return res


def measure_scalar(sol, f, *args):

    res = 0. + 0.J
    for i, e in enumerate(sol['evec']):
        res += numpy.conj(e) * f(i, sol, *args) * e


    return res / sol['L']
