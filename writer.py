import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, sys, struct, os, operator
import matplotlib.pyplot as plt
import solver, data_processing

INDEX_KEYS = ['bc', 'N', 'L', 'lambda', 'index']
SORTING_SCHEME_DEFAULT = ['bc', 'L', 'N', 'gamma']
SORTING_SCHEME_DEFAULT = ['L', 'lambda', 'bc', 'q1']
OUTPUT_SCHEME_DEFAULT = ['bc', 'lambda', 'L', 'cX_real', 'cX_expected', 'e', 'e_expected']
OUTPUT_SCHEME_DEFAULT = ['N', 'bc', 'lambda', 'L', 'cX_real', 'e', 'e_infinity']
OUTPUT_SCHEME_DEFAULT = ['L', 'e']
FLOAT_PRECISION_DEFAULT = 16
FLOAT_WIDTH_DEFAULT = 26
INTEGER_WIDTH_DEFAULT = 6
OUTPUT_SETTINGS= {
        'lambda' : {'p1' : 8, 'p2' : 2},
        'bc' : {'values' : {0 : 'open', 1 : 'periodic'},
            'p1' : 9},
        'N' : {'p1' : 3},
        'L' : {'p1' : 4},
        'index' : {'p1' : 4},
        }

def write_any(x, name, p1=None, p2=None):
    if data_processing.is_iterable(x): x = x[0]
    if name in OUTPUT_SETTINGS:
        ost = OUTPUT_SETTINGS[name]
        if 'values' in ost:
            if x in ost['values']: x = ost['values'][x]
        if 'p1' in ost: p1 = ost['p1']
        if 'p2' in ost: p2 = ost['p2']

    dtype = 'i'
    just = ''
    if isinstance(x, complex):
        x = x.real
    if isinstance(x, int):
        if not p1: p1 = INTEGER_WIDTH_DEFAULT
        dtype = 'i'
    if isinstance(x, float):
        if not p1: p1 = FLOAT_WIDTH_DEFAULT
        if not p2: p2 = FLOAT_PRECISION_DEFAULT
        dtype = 'f'
    if not p1: p1 = FLOAT_WIDTH_DEFAULT
    if isinstance(x, str):
        dtype = 's'
        just = '+'
    z = '%' + just + str(p1) + dtype
    if isinstance(x, float):
        z = '%' + str(p1) + '.' + str(p2) + dtype
    #print z, n
    return z % x

def write_name(x, n, p1=None):
    if n in OUTPUT_SETTINGS:
        ost = OUTPUT_SETTINGS[n]
        if 'values' in ost:
            if x in ost['values']: x = ost['values'][x]
        if 'p1' in ost: p1 = ost['p1']
        if 'p2' in ost: p2 = ost['p2']
    if data_processing.is_iterable(x): x = x[0]
    if isinstance(x, int):
        if not p1: p1 = INTEGER_WIDTH_DEFAULT
    if isinstance(x, float):
        if not p1: p1 = FLOAT_WIDTH_DEFAULT
    if not p1: p1 = FLOAT_WIDTH_DEFAULT
    z = '%' + str(p1) + 's'
    #print z, n
    return z % n

def generate_output(s, scheme):
    res = ''
    for x in scheme:
        if not data_processing.is_iterable(x): x = [x]
        res += write_any(s[x[0]], x[0], *x[1:])
    return res

def generate_header(s, scheme):
    res = ''
    for x in scheme:
        if not data_processing.is_iterable(x): x = [x]
        res += write_name(s[x[0]], x[0], *x[1:])
    return res

def sort_sols(sol_list, scheme):
    for x in scheme:
        if x in sol_list[0]:
            sol_list = sorted(sol_list, key=lambda y: y[x])
    return sol_list

def make_ident(c):
    keys = INDEX_KEYS
    res = []
    for k in keys: 
        if k in c:
            res.append(float(c[k]))
    return tuple(res)
