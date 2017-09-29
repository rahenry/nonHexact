import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, sys, struct, os, argparse, scipy.optimize
import matplotlib.pyplot as plt
import solver, writer

parser = argparse.ArgumentParser()
parser.add_argument('inputs', nargs = '+')
args = parser.parse_args(sys.argv[1:])
POLYFIT_DEGREE = 5
POLYFIT_INF = 1000

def diff_function(x, y):
    return abs(x-y)/(abs(x)+abs(y))

def expfit_func(x, a, b, c):
    return a * numpy.exp(-b * x) + c

def expfit_series(sols, key):
    data = []
    for s in sols[-10:]:
        data.append(s[key])

    popt, pcov = scipy.optimize.curve_fit(expfit_func, range(len(data)), data)
    return popt

def polyfit_series(sols, key):
    data = []
    for s in sols:
        data.append(s[key])

    fit = numpy.polyfit(range(len(data)), data, POLYFIT_DEGREE)
    return fit

def extrapolate(sols, key):
    if key in writer.INDEX_KEYS:
        return 'inf'
    fit = expfit_series(sols, key)
    return expfit_func(POLYFIT_INF, *fit)

def is_numeric(s):
    numtypes = [int, float]
    for ntype in numtypes:
        try:
            ntype(s)
            return True
        except ValueError:
            pass
    return False

def is_line_numeric(l):
    for x in l.split():
        if not is_numeric(x): return False
    return True

def get_headings(l):
    l = l.strip().split()
    return l

def compare_files(file_list):
    line_data = {}
    headings = ''
    for fname in file_list:
        f = open(fname)
        split = f.readlines()
        if not is_line_numeric(split[0]):
            if headings == '':
                headings = get_headings(split[0])
        if not is_line_numeric(split[0]):
            split = split[1:]
        line_data[fname] = split
        f.close()

    if headings == '':
        headings = ['index']
        heading_index = 1
        z = line_data[fname][0]
        z = z.strip().split()
        while len(headings) < len(z):
            headings.append('q' + str(heading_index))
            heading_index += 1

    sols = {}
    for fname in file_list:
        sols[fname] = {}
        for l in line_data[fname]:
            new_sol = {}
            l = l.strip().split()
            for name, x in zip(headings, l):
                if '.' not in x:
                    new_sol[name] = int(x)
                else: new_sol[name] = float(x)
            sols[fname][writer.make_ident(new_sol)] = new_sol

    sols_composite = {}
    output_scheme = []
    for h in headings:
        if h in writer.INDEX_KEYS: output_scheme.append(h)
        else:
            file_index = 1
            for fname in file_list:
                output_scheme.append(h+ '_' + str(file_index))
                if file_index > 1:
                    output_scheme.append(h + '_diff' + str(file_index))
                file_index += 1

    composite_sols = {}
    sol_list = []
    for ident, sol in sols[file_list[0]].iteritems():
        new_sol = dict(sol)
        for key in sol:
            if key in writer.INDEX_KEYS: continue

            file_index = 1
            for fname in file_list:
                new_sol[key + '_' + str(file_index)] = sols[fname][ident][key]
                if file_index > 1:
                    v1 = sols[fname][ident][key]
                    v2 = sols[file_list[0]][ident][key]
                    new_sol[key + '_diff' + str(file_index)] = diff_function(v1, v2)
                file_index += 1

        new_sol['output'] = writer.generate_output(new_sol, output_scheme)
        composite_sols[ident] = new_sol
        sol_list.append(new_sol)

    #print composite_sols
    res = writer.generate_header(random.choice(composite_sols.values()), output_scheme) + '\n'
    sol_list = writer.sort_sols(sol_list, writer.SORTING_SCHEME_DEFAULT)
    for sol in sol_list:
        res += sol['output'] + '\n'

    sol_extrapolated = {}
    for s in output_scheme:
        sol_extrapolated[s] = extrapolate(sol_list, s)
    print res
    #print sol_extrapolated
    #print writer.generate_output(sol_extrapolated, output_scheme)

if len(args.inputs) > 1:
    compare_files(args.inputs)
else:
    compare_files(args.inputs + args.inputs)

