import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, sys, struct, os
import matplotlib.pyplot as plt
import solver

DATA_MARKER_START = '*******'
DATA_MARKER_END = '___#*$'
ENCODE_EVECS = False

def is_iterable(x):
    return hasattr(x, '__iter__')

def read_input_file(input_file):
    f = open(input_file)
    res = {}
    for l in f.readlines():
        l = l.split()
        if len(l) == 0: continue
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
    if 'evec' in key: return ''
    res = ''
    val = sol[key]
    if not is_iterable(val):
        val = [val]
    for v in val:
        if isinstance(v, complex):
            return ''
    res += struct.pack('q', len(key))
    res += key
    if isinstance(val[0], int): res += '0'
    else: res+= '1'

    dat = ''
    for x in val:
        if isinstance(x, int):
            dat += struct.pack('q', x)
        else:
            dat += struct.pack('d', x)
    res += struct.pack('q', len(dat))
    res += dat
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
    res = struct.pack('q', len(res)) + res
    return res

def decode_solutions(data_file):
    if not os.path.isfile(data_file): return {}
    f = open(data_file)
    d = f.read()
    index = 0
    res = {}
    while index < len(d):
        if len(d[index:]) < 5: break
        size = struct.unpack('q', d[index:index+8])[0]
        new_sol = decode_solution(d[index:index+8+size])
        index += 8+size
        res[new_sol['ident']] = new_sol
    return res

def decode_solution(d):
    res = {}
    index = 8
    while index < len(d):
        if len(d[index:]) < 5: break
        subdata = d[index:]

        name_len = struct.unpack('q', subdata[0:8])[0]
        name = struct.unpack(str(name_len)+'s', subdata[8:8+name_len])[0]

        flag = d[index+8+name_len]
        if flag == '0': dtype = 'q'
        else: dtype = 'd'

        data_len = struct.unpack('q', subdata[9+name_len:17+name_len])[0]
        n_data = data_len / 8
        res[name] = struct.unpack(dtype * n_data, subdata[17+name_len:17+name_len+data_len])

        if n_data == 1: 
            res[name] = res[name][0]

        index += 17+name_len+data_len

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

        for c in combs:
            ident = solver.make_ident(c)
            sols[input_file][ident] = c
            sols[input_file][ident]['ident'] = ident

        data_file = input_file + '.data'
        sols[input_file].update(decode_solutions(data_file))
        #print decode_solutions(data_file)

        #for ident, s in sols[input_file].iteritems():
            #for x in s: print x
        for ident, s in sols[input_file].iteritems():
            if 'e' not in s: 
                print "Solving " + str(ident)
                solver.solve(s)
            #enc = encode_solution(s)
            #dec = decode_solution(enc)
            s['e_infinity'] = solver.exact_eigenvalue(s['L'], s['N'], s['lambda'])

        f = open(data_file, 'w+')
        for ident, s in sols[input_file].iteritems():
            f.write(encode_solution(s))
        f.close()

    
    return sols




