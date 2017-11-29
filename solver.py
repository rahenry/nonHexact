import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special
import matplotlib.pyplot as plt
import writer

SIGMA_OFFSET = 0 

def sigma(i, j, N):
    return (i / (N ** j)) % N

def normalised(s):
    norm = numpy.linalg.norm(s)
    if norm < 1E-12: return s
    return s / norm

def index(initial, j_modified, sigma_j_f, N):
    sigma_j_i = sigma(initial, j_modified, N)
    return initial + ((N ** j_modified)  * (sigma_j_f - sigma_j_i))
    
def cZ_expected(sol):
    gamma = sol['lambda']
    p = -2.
    if sol['bc'] == 1: p = -4.
    if gamma <= 1: return 0.

    else: return (1. - (gamma ** -p))

def cX_expected(sol):
    gamma = sol['lambda']
    if gamma > 1.: return 0.
    return pow(1. - gamma*gamma, 0.25)

def calc_charge_mapping(sol, i):
    res = 0
    L = sol['L']
    N = sol['N']
    for j in range(L):
        sigma_j = sigma(i, j, sol['N'])
        z = (sigma_j - 1) % N
        res += (N**j) * z
    return res

def generate_mappings(sol):
    N = sol['N']
    L = sol['L']
    bc = sol['bc']
    gamma = sol['lambda']

    omega = cmath.exp(2.0 * math.pi * 1.0J / N)
    M = N**L # size of hilbert space
    Q = range(M) # indices of all states
    s0 = []
    X_mapping = []
    X_conj_mapping = []
    Z_mapping = []
    charge_mapping = []

    for j in range(L):
        X_mapping.append([])
        X_conj_mapping.append([])
    for i in range(N ** L):
        if (bc == 3):
            charge_mapping.append(calc_charge_mapping(sol, i))
            #print (i, calc_charge_mapping(sol, i))
        #s0.append(1.0 + 0.51 * random.random())
        s0.append(1.0)
        for j in range(L):
            sigma_j = sigma(i, j, N)
            if (sigma_j < 1):
                X_mapping[j].append(index(i, j, sigma_j+1, N))
                X_conj_mapping[j].append(index(i, j, 0, N))
            elif (sigma_j < N - 1):
                X_mapping[j].append(index(i, j, sigma_j+1, N))
                X_conj_mapping[j].append(index(i, j, sigma_j-1, N))
            else:
                X_mapping[j].append(index(i, j, 0, N))
                X_conj_mapping[j].append(index(i, j, sigma_j-1, N))


        Z_coef = 0.
        R = [1]
        if sol['model'] == 'SICP':
            R = range(1, N)
        for r in R:
            u = 1.
            if sol['model'] == 'SICP':
                u = SICP_factor(r, N)

            for j in range(L-1):
                Z_coef += u * pow((omega ** (sigma(i, j, N) -SIGMA_OFFSET )) * ((omega ** (sigma(i, j+1, N) - SIGMA_OFFSET)) ** (N-1)), r)

            if (bc == 1 or bc == 2):
                j = L-1
                g = 1.0
                if bc == 2: g = omega
                Z_coef += g * u * pow((omega ** (sigma(i, j, N) - SIGMA_OFFSET)) * ((omega ** (sigma(i, 0, N) -SIGMA_OFFSET)) ** (N-1)), r)

        Z_mapping.append(Z_coef)
    return (X_mapping, X_conj_mapping, Z_mapping, charge_mapping)

def SICP_factor(r, N):
    u = numpy.exp(math.pi * 1.J * (2.*r-N)/2./N)/numpy.sin(math.pi*r/N)
    return u

def solve(sol):
    N = sol['N']
    L = sol['L']
    gamma = sol['lambda']
    bc = sol['bc']
    TOLERANCE = sol['TOLERANCE']
    omega = cmath.exp(2.0 * math.pi * 1.0J / N)
    M = N**L # size of hilbert space
    Q = range(M) # indices of all states

    ident = writer.make_ident(sol)

    X_mapping, X_conj_mapping, Z_mapping, charge_mapping = generate_mappings(sol)

    s0 = []
    rows = []
    cols = []
    data = []
    count = 0

    R = range(1, N)
    if (sol['model'] == 'SICP'):
        for i in Q:
            rows.append(i)
            cols.append(i)
            x = 0.
            x -= gamma * (Z_mapping[i])
            data.append(x)

            for j in range(L):
                for r in R:
                    k = i
                    z = 0
                    while (z < r):
                        k = X_mapping[j][k]
                        z += 1

    
                    Z = 0
                    K = i
                    while (Z < r+5):
                        #print K
                        K = X_mapping[j][K]
                        Z += 1
                    #print '...'
                    rows.append(i)
                    cols.append(k)
                    data.append(-SICP_factor(r, N))


    else:

        for i in Q:
            rows.append(i)
            cols.append(i)
            data.append(-gamma * Z_mapping[i])

            for x in range(L):
                j = X_mapping[x][i]
                rows.append(i)
                cols.append(j)
                data.append(-1)
            if sol['bc'] == 3:
                rows.append(i)
                j = charge_mapping[i]
                cols.append(j)
                g = (omega ** (sigma(i, L-1, N) - SIGMA_OFFSET)) * ((omega ** (sigma(i, 0, N) -SIGMA_OFFSET)) ** (N-1))
                g = (omega ** (sigma(i, 0, N) - SIGMA_OFFSET)) * ((omega ** (sigma(i, L-1, N) -SIGMA_OFFSET)) ** (N-1))
                a = 1.
                data.append(-gamma * g)

    e = []
    if sol['eig_method'] == 'sparse':
        m = scipy.sparse.coo_matrix((data, (rows, cols)), (M, M))
        e = scipy.sparse.linalg.eigs(m, k=2, maxiter = 10000000, tol=TOLERANCE, which='SR')
    elif sol['eig_method'] == 'sparse_eigs':
        m = scipy.sparse.coo_matrix((data, (rows, cols)), (M, M))
        neigs = 1200
        e = scipy.sparse.linalg.eigs(m, k=neigs, maxiter = M*10, tol=TOLERANCE, which='LM', return_eigenvectors=False)
        e = (e, [])
    else:
        m = numpy.zeros((M,M), dtype=numpy.complex_)
        for (i,j,d) in zip(rows, cols, data):
            m[i, j] = d
        e = scipy.linalg.eig(m)
        for q in e[0]:
            1
            #print q
    #e = scipy.sparse.linalg.eigs(m, k=5, maxiter = 100000, tol=1E-2, which='SM')
    #for q in e[0]:
        #print q

    cZ = 0
    cX = 0
    evec = []

    if not sol['eig_method'] == 'sparse_eigs':
        evec_ind = 0
        if (e[0][1] < e[0][0]): evec_ind = 1
        #print e[0][evec_ind]
        for z in e[1]:
            evec.append(z[evec_ind])
        evec = normalised(evec)

        bra = numpy.conj(evec)
        for i in Q:
            #cZ += bra[i] * pow(omega, sigma(i,0)-SIGMA_OFFSET) * pow(omega, (N-1) * (sigma(i, L-1)-SIGMA_OFFSET)) * evec[i]
            cZ += bra[i] * pow(omega, sigma(i,0, N)-SIGMA_OFFSET) * pow(omega, (1) * (sigma(i, L-1, N)-SIGMA_OFFSET)) * evec[i]

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


    ind = 0
    for x in e[0]:
        break
        if abs(x) < 1E-5:
        #if ind == 0:
            xdata = []
            ydata = []
            for y in e[1][ind]:
                xdata.append(y.real)
                ydata.append(y.imag)
            plt.plot(xdata, ydata, linestyle='', marker='o', ms=3)
            plt.axes().set_aspect('equal', 'datalim')
            #print x
            #print e[1][ind]
            #plt.show()

            break
        ind += 1

    plt.clf()
    ind = 0
    xdata = []
    ydata = []
    count1 = 0
    for x in e[0]:
        #print x
        if abs(x) < 1E-3:
            count1 += 1
        #if ind == 0:
            #print x
            for y in e[1][ind]:
                xdata.append(y.real)
                ydata.append(y.imag)
        ind += 1
    print 'count1 = ', count1
    plt.plot(xdata, ydata, linestyle='', marker='o', ms=3)
    plt.axes().set_aspect('equal', 'datalim')
    #plt.show()
    sol.update({
            'energy' : min(e[0]),
            'evec' : evec,
            'cZ' : cZ,
            'cX': cX,
            'N' : N,
            'L' : L,
            'gamma' : gamma,
            'bc' : bc,
            #'ident' : ident,
            'e' : min(e[0]) / L,
            'energy_expected' : exact_eigenvalue(sol),
            'eigenvalues' : e[0],
            })

    sol.update({
        'cZ_expected': cZ_expected(sol),
        'cX_expected': cX_expected(sol),
        'e_infinity' : sol['energy_expected']
        })
    
    return sol

def eps(j, L, N, gamma):
    kj = 2.*j*math.pi/(2.*L+1.)
    return pow(1. + gamma**N + 2.*(gamma**(.5*N)) * math.cos(kj), 1./N)


def exact_eig2(L, N, gamma):
    res = 0.0
    for j in range(L):
        res -= eps(j+1, L, N, gamma)
    return res / L
def exact_eigenvalue(sol):
    L = sol['L']
    N = sol['N']
    gamma = sol['lambda']
    if False:
    #if gamma == 1.:
        res = 0.0
        for j in range(L):
            res -= eps(j+1, L, N, gamma)
        return res / L

    if sol['model'] == 'SICP':
        res = 0.
        for l in range(1, N):
            res += (1.+gamma) * scipy.special.hyp2f1(-0.5, float(l)/N, 1., 4.*gamma/((1.+gamma)**2.))
        return -res
    else:
        if (gamma >= 1.):
            return -(1+gamma**N)**(1./N) * scipy.special.hyp2f1(-1./2./N, (-1./N+1.)/2., 1., 4.*gamma**N/((1.+gamma**N)**2))
            return -(1.-gamma**N)**(1./N) * scipy.special.hyp2f1(-1./N,1.-1./N,1.,1./(1.-gamma**(-N)))
            return -gamma * scipy.special.hyp2f1(-1./N, -1./N, 1., (1./gamma)**N)

        return -scipy.special.hyp2f1(-1./N, -1./N, 1., gamma**N)

