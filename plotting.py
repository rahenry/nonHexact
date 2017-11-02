import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, sys, struct, os, operator
import matplotlib.pyplot as plt
import solver, data_processing, writer, data_processing
from matplotlib.backends.backend_pdf import PdfPages
SORTING_SCHEME_DEFAULT = writer.SORTING_SCHEME_DEFAULT

labels = ['N', 'bc', 'L', 'lambda']

def plot_eigenvalues(data):
    for input_file, d in data.iteritems():
        print input_file
        data_dir = os.path.abspath(input_file + '_dir')
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
        plot_dir = os.path.abspath(data_dir + '/plots')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        pp = PdfPages(os.path.abspath(plot_dir + '/eigenvalues.pdf'))

        sol_list = []
        for ident, sol in d.iteritems():
            sol_list.append(sol)
        sol_list = writer.sort_sols(sol_list, SORTING_SCHEME_DEFAULT)
        for sol in sol_list:
            print ident
            if sol['eig_method'] == 'full' or sol['eig_method'] == 'sparse_eigs':
                xdata = []
                ydata = []
                for e in sol['eigenvalues']:
                    xdata.append(e.real)
                    ydata.append(e.imag)

                plt.plot(xdata, ydata, linestyle='', marker='o', ms=3)
                plt.axes().set_aspect('equal', 'datalim')

                labeltext = ''
                for l in labels:
                    if l in sol:
                        #labeltext += writer.write_name(l, l, p1=len(l), justify_left=True)
                        labeltext += l
                        labeltext += ' = ' + writer.write_any(sol[l], l, justify_left=True) + '\n'
                if labeltext != '':
                    plt.text(0.1, 0.9, labeltext, ha='left', va='top', transform=plt.axes().transAxes)
                pp.savefig()
                #plt.savefig(pp, format='pdf')
                if len(sol_list) == 1: plt.show()
                plt.clf()
        pp.close()
