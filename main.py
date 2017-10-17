import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, sys, struct, os, operator
import matplotlib.pyplot as plt
import solver, data_processing, writer, measure

OUTPUT_SCHEME_DEFAULT = writer.OUTPUT_SCHEME_DEFAULT
SORTING_SCHEME_DEFAULT = writer.SORTING_SCHEME_DEFAULT

data = data_processing.read_inputs().iteritems()
header = ''
sol_list = []
for input_file, d in data:
    for ident, s in d.iteritems():
        if header == '': header = writer.generate_header(s, OUTPUT_SCHEME_DEFAULT)
        s['output'] = writer.generate_output(s, OUTPUT_SCHEME_DEFAULT)
        sol_list.append(s)

sol_list = writer.sort_sols(sol_list, SORTING_SCHEME_DEFAULT)

output = header
for s in sol_list:
    output += '\n' + s['output']
output += '\n'

output_file = 'latest_output'
f = open(output_file, 'w+')
f.write(output)
f.close()


k = 11.111
#for s in sol_list:
    #if s['L'] != 6: continue
    #print s['ident'], measure.measure_scalar(s, measure.momentum_3, k)
    #print s['ident'], measure.measure_scalar(s, measure.momentum_1, k)

print '-------'

#for s in sol_list:
    #print s['ident'], measure.measure_scalar(s, measure.current_1, k)

for sol in sol_list:
    L = sol['L']
    N = sol['N']
    res = ''
    total = 0.
    for j in range(0, L-1):
        se = 0. +0.J
        for i, e in enumerate(sol['evec']):
            s1 = solver.sigma(i, j, N)
            s2 = solver.sigma(i, j+1, N)
            diff = (s2 - s1) % N
            se += diff * e * numpy.conj(e)
        total += se
        #res += '%9.2E %9.2E' % (se.real, se.imag) + ' '
        res += '%9.2E' % se.real + ' '
    if sol['bc'] == 1:
        se = 0. +0.J
        for i, e in enumerate(sol['evec']):
            s1 = solver.sigma(i, L-1, N)
            s2 = solver.sigma(i, 0, N)
            diff = (s2 - s1) % N
            se += diff * e * numpy.conj(e)
        #res += '%9.2E %9.2E' % (se.real, se.imag) + ' '
        res += '%9.2E' % se.real + ' '
        total += se

    div = L
    if sol['bc'] == 0 : div += 1
    print L, sol['bc'], sol['lambda'], ' | ', total.real / div
    #print res
    #print ''
        

