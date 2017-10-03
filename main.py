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

output_file = 'latest_output'
f = open(output_file, 'w+')
f.write(output)
f.close()


k = 9.9
for s in sol_list:
    if s['L'] != 9: continue
    print s['ident'], measure.measure_scalar(s, measure.momentum_1, k)
