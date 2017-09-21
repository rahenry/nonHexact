import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, sys, struct, os
import matplotlib.pyplot as plt
import solver

DATA_MARKER_START = '***'
DATA_MARKER_END = '___'

def is_iterable(x):
    return hasattr(x, '__iter__')

def read_input_file(input_file):
    f = open(input_file)
    res = {}
    for l in f.readlines():
        l = l.split()
        res[l[0]] = []
        for x in l[1:]:
            if x == '=': continue
            if ':' in x:
                y = x.split(':')
                res[l[0]] += range(int(y[0]), int(y[1])+1)
            else:
                try: 
                    res[l[0]].append(int(x))
                    continue
                except ValueError: pass
                try: 
                    res[l[0]].append(float(x))
                    continue
                except ValueError: res[l[0]].append(x)
    return res

def encode_data(sol, key):
    res = ''
    val = sol[key]
    if not is_iterable(val):
        val = [val]
    for v in val:
        if isinstance(v, complex):
            return ''
    res += key + ' '
    if isinstance(val[0], int): res += 'int'
    res += DATA_MARKER_START
    for x in val:
        if isinstance(x, int):
            res += struct.pack('q', x)
        else:
            res += struct.pack('d', x)
    res += DATA_MARKER_END
    return res

def recode_complex(sol, key):
    val = sol[key]
    if not is_iterable(val):
        val = [val]
    res = {key + '_real': [],
            key + '_imag' : [],
            }
    for v in val:
        res[key+'_real'].append(v.real)
        res[key+'_imag'].append(v.imag)
    return res

def encode_solution(sol):
    res = ''

    new_data = {}
    for key, val in sol.iteritems():
        if not is_iterable(val):
            val = [val]
        is_complex = False
        for v in val:
            if isinstance(v, complex):
                is_complex = True
        if is_complex: new_data.update(recode_complex(sol, key))
    sol.update(new_data)
    for key in sol:
        res += encode_data(sol, key)
    return res

def decode_solution(d):
    res = {}
    for l in d.split(DATA_MARKER_END):
        l = l.split(DATA_MARKER_START)
        if 'int' in l[0]:
            dtype = 'q'
        else: dtype = 'd'
        
        n_data = len(l[-1]) / 8
        if n_data == 0: continue
        offset = 0
        name = l[0].split()[0]
        if (n_data == 1):
            res[name] = struct.unpack(dtype, l[-1])[0]
            continue
        if (n_data > 0):
            res[name] = struct.unpack(dtype * n_data, l[-1])

    new_data = {}
    for key in res:
        if '_real' in key:
            basekey = key.replace('_real', '')
            if is_iterable(res[key]):
                new_data[basekey] = []
                for (x, y) in zip(res[key], res[basekey+'_imag']):
                    new_data[basekey].append(x + y*1.J)
            else:
                new_data[basekey] = res[key] + res[basekey+'_imag']*1.J
    res.update(new_data)

    return res

def process_combinations(old, new):
    res = []
    for x in old:
        for y in new[1]:
            res.append(dict(x))
            res[-1][new[0]] = y
    return res

def read_inputs():
    defaults = read_input_file('defaults')
    input_files = sys.argv[1:]
    sols = {}
    for input_file in input_files:
        params = dict(defaults)
        params.update(read_input_file(input_file))
        sols[input_file] = {}
        combs = [[]]
        for key, val in params.iteritems():
            combs = process_combinations(combs, (key, val))

        ident = 0
        for c in combs:
            sols[input_file][ident] = c
            sols[input_file][ident]['ident'] = ident
            ident += 1

        data_file = input_file + '.data'
        new_data = {}
        if os.path.isfile(data_file):
            f = open(data_file)
            for d in f.read().split('\n\n'):
                sol = decode_solution(d)
                try:
                    new_data[sol['ident']] = sol
                except KeyError:
                    pass
            f.close()
        sols[input_file].update(new_data)

        #for ident, s in sols[input_file].iteritems():
            #for x in s: print x
        for ident, s in sols[input_file].iteritems():
            if 'e' not in s: 
                print 55
                solver.solve(s)
            enc = encode_solution(s)
            dec = decode_solution(enc)

        f = open(data_file, 'w+')
        for ident, s in sols[input_file].iteritems():
            f.write(encode_solution(s) + '\n\n')
        f.close()





