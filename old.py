import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special
import matplotlib.pyplot as plt

TOLERANCE = 1E-5
SIGMA_OFFSET = 0 

Ns = [3]
Ls = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
Ls = [2, 3, 4, 5]
Ls = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
Ls = [8]
Ls = [2, 3, 4]
gammas = [0.1, 0.5, 1., 2., 10.]
gammas = [0.1, 0.5, 1., 2.]
bcs = [0]
bcs = [0, 1]

def sigma(i, j):
    return (i / (N ** j)) % N

def normalised(s):
    norm = numpy.linalg.norm(s)
    if norm < 1E-12: return s
    return s / norm

def index(initial, j_modified, sigma_j_f):
    sigma_j_i = sigma(initial, j_modified)
    return initial + ((N ** j_modified)  * (sigma_j_f - sigma_j_i))
    
def cZ_expected(sol):
    gamma = sol['gamma']
    p = -2.
    if sol['bc'] == 1: p = -4.
    if gamma <= 1: return 0.

    else: return (1. - (gamma ** -p))

def cX_expected(sol):
    gamma = sol['gamma']
    if gamma > 1.: return 0.
    return pow(1. - gamma*gamma, 0.25)

def generate_mappings(N, L, gamma, bc):
    omega = cmath.exp(2.0 * math.pi * 1.0J / N)
    M = N**L # size of hilbert space
    Q = range(M) # indices of all states
    s0 = []
    X_mapping = []
    X_conj_mapping = []
    Z_mapping = []
    for j in range(L):
        X_mapping.append([])
        X_conj_mapping.append([])
    for i in range(N ** L):
        #s0.append(1.0 + 0.51 * random.random())
        s0.append(1.0)
        for j in range(L):
            sigma_j = sigma(i, j)
            if (sigma_j < 1):
                X_mapping[j].append(index(i, j, sigma_j+1))
                X_conj_mapping[j].append(index(i, j, 0))
            elif (sigma_j < N - 1):
                X_mapping[j].append(index(i, j, sigma_j+1))
                X_conj_mapping[j].append(index(i, j, sigma_j-1))
            else:
                X_mapping[j].append(index(i, j, 0))
                X_conj_mapping[j].append(index(i, j, sigma_j-1))

        Z_coef = 0.
        for j in range(L-1):
            Z_coef += (omega ** (sigma(i, j) -SIGMA_OFFSET )) * ((omega ** (sigma(i, j+1) - SIGMA_OFFSET)) ** (N-1))

        if (bc):
            j = L-1
            Z_coef += (omega ** (sigma(i, j) - SIGMA_OFFSET)) * ((omega ** (sigma(i, 0) -SIGMA_OFFSET)) ** (N-1))
        Z_mapping.append(Z_coef)
    return (X_mapping, X_conj_mapping, Z_mapping)

def solve(N, L, gamma, bc):
    omega = cmath.exp(2.0 * math.pi * 1.0J / N)
    M = N**L # size of hilbert space
    Q = range(M) # indices of all states

    ident = (N, L, gamma, bc)

    X_mapping, X_conj_mapping, Z_mapping = generate_mappings(N, L, gamma, bc)

    s0 = []
    rows = []
    cols = []
    data = []
    count = 0

    for i in Q:
        rows.append(i)
        cols.append(i)
        data.append(-gamma * Z_mapping[i])

        for x in range(L):
            j = X_mapping[x][i]
            rows.append(i)
            cols.append(j)
            data.append(-1)
    #print "Data created, building matrix..."
    m = scipy.sparse.coo_matrix((data, (rows, cols)), (M, M))
    #print "Matrix built, starting diag..."

    e = scipy.sparse.linalg.eigs(m, k=2, maxiter = 10000000, tol=TOLERANCE, which='SR')

    evec_ind = 0
    if (e[0][1] < e[0][0]): evec_ind = 1
    #print e[0][evec_ind]
    evec = []
    for z in e[1]:
        evec.append(z[evec_ind])
    cZ = 0
    cX = 0
    evec = normalised(evec)

    bra = numpy.conj(evec)
    for i in Q:
        #cZ += bra[i] * pow(omega, sigma(i,0)-SIGMA_OFFSET) * pow(omega, (N-1) * (sigma(i, L-1)-SIGMA_OFFSET)) * evec[i]
        cZ += bra[i] * pow(omega, sigma(i,0)-SIGMA_OFFSET) * pow(omega, (1) * (sigma(i, L-1)-SIGMA_OFFSET)) * evec[i]

    ket = numpy.zeros_like(evec)
    x_inds = [0, L-1]
    for i in Q:
        #ket[X_conj_mapping[0][i]] += evec[i]
        ket[X_mapping[0][i]] += evec[i]
    ket_old =  numpy.array(ket)

    if True:
        ket = numpy.zeros_like(ket_old)
        for i in Q:
            ket[X_mapping[L-1][i]] += ket_old[i]

    for q in range(N-1):
        continue
        ket = numpy.zeros_like(ket_old)
        for i in Q:
            #ket[X_conj_mapping[L-1][i]] += ket[i]
            ket[X_mapping[L-1][i]] += ket_old[i]
        ket_old = numpy.array(ket)
    for i in Q:
        cX += bra[i] * ket[i]


    sol = {
            'energy' : min(e[0]),
            'evec' : evec,
            'cZ' : cZ,
            'cX': cX,
            'N' : N,
            'L' : L,
            'gamma' : gamma,
            'bc' : bc,
            'ident' : ident,
            'e' : min(e[0]) / L,
            }

    sol.update({
        'cZ_expected': cZ_expected(sol),
        'cX_expected': cX_expected(sol),})
    
    return sol
results = {}
counter = 0
print_list = ['bc', 'gamma', 'L', 'cZ', 'cX', 'cZ_expected', 'cX_expected', 'e', 'e_exact']
print_mappings = {
        'bc' : 'BCs',
        'gamma' : 'lambda',
        }

for N in Ns:
    for L in Ls:
        for gamma in gammas:
            for bc in bcs:
                if L == Ls[-1]:


                    comp = float(counter) / float(len(gammas) * len(bcs))
                    print "%.5f complete" % comp
                    counter = counter + 1
                sol = solve(N, L, gamma, bc)
                results[sol['ident']] = sol

def eps(j, L, N, gamma):
    kj = 2.*j*math.pi/(2.*L+1.)
    return pow(1. + gamma**N + 2.*(gamma**(.5*N)) * math.cos(kj), 1./N)


def exact_eigenvalue(L, N, gamma):
    if gamma == 1.:
        res = 0.0
        for j in range(L):
            res -= eps(j+1, L, N, gamma)
        return res / L
    if (gamma > 1.):
        return -gamma * scipy.special.hyp2f1(-1./N, -1./N, 1., (1./gamma)**N)
    else:
        return -scipy.special.hyp2f1(-1./N, -1./N, 1., gamma**N)

#print '%-3s %10s %10s %10s %10s %10s' % ('L', 'lambda', 'cZ', 'cX', 'e', 'e_exact')

data_file = 'newest_data.txt'
f = open(data_file, 'w+')
print '%-10s %5s %-3s %8s %8s %8s %8s %8s %8s' % ('BCs', 'lambda', 'L', 'cZ_real', 'cZ_imag', 'cX_real', 'cX_imag', 'e', 'e_exact')
f.write('%-10s %5s %-3s %8s %8s %8s %8s %8s %8s\n' % ('BCs', 'lambda', 'L', 'cZ_real', 'cZ_imag', 'cX_real', 'cX_imag', 'e', 'e_exact'))
for gamma in gammas:
    for N in Ns:
        for bc in bcs:
            for L in Ls:
                idn = (N, L, gamma, bc)
                for ident, r in results.iteritems():
                    if r['ident'] != idn: continue
                    #print '%-3i %10f %10f %10f %10f %10f' % (r['L'], r['gamma'], r['cZ'], cZ_expected(r['gamma']), r['e'], exact_eigenvalue(r['L'], r['N'], r['gamma']))
                    cZ = r['cZ']
                    cX = r['cX']
                    bc_string = 'closed'
                    if bc == 1: bc_string = 'periodic'

                    print '%-10s %5.2f %-3i %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f' % (bc_string,r['gamma'], r['L'], cZ.real, cZ.imag, cX.real, cX.imag, r['e'].real, exact_eigenvalue(r['L'], r['N'], r['gamma']))
                    f.write('%-10s %5.2f %-3i %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n' % (bc_string,r['gamma'], r['L'], cZ.real, cZ.imag, cX.real, cX.imag, r['e'].real, exact_eigenvalue(r['L'], r['N'], r['gamma'])))




