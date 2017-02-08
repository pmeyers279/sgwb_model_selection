#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals, print_function
__doc__ = """
Script that does default visualizations (marginal plots, 1-d and 2-d).

Author: Johannes Buchner (C) 2013
edited by Thomas Callister and Pat Meyers
"""
import numpy
from numpy import exp, log
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import sys, os
import json
import pymultinest

def marginals(prefix):
    print('model "%s"' % prefix)
    parameters = json.load(open(prefix + 'params.json'))
    n_params = len(parameters)

    # Get posterior data
    a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename = prefix)
    s = a.get_stats()
    json.dump(s, open(prefix + 'stats.json', 'w'), indent=4)

    # Print posterior data
    print('  marginal likelihood:')
    print('    ln Z = %.1f +- %.1f' % (s['global evidence'], s['global evidence error']))
    print('  parameters:')

    # Loop across parameters, obtain median and std.dev
    for p, m in zip(parameters, s['marginals']):
            lo, hi = m['1sigma']
            med = m['median']
            sigma = (hi - lo) / 2
            i = max(0, int(-numpy.floor(numpy.log10(sigma))) + 1)
            fmt = '%%.%df' % i
            fmts = '\t'.join(['    %-15s' + fmt + " +- " + fmt])
            print(fmts % (p, med, sigma))

    # Initialize PlotMarginal object
    print('creating marginal plot ...')
    p = pymultinest.PlotMarginal(a)

    # Get posterior and mode info
    values = a.get_equal_weighted_posterior()
    assert n_params == len(s['marginals'])
    modes = s['modes']

    # Set $D = 1 if more than 20 parameters, 2 otherwise
    dim2 = os.environ.get('D', '1' if n_params > 20 else '2') == '2'
    nbins = 100 if n_params < 3 else 20

    # If $D==2...
    if dim2:

            # Initialize plot
            fig = plt.figure(figsize=(5*n_params, 5*n_params))

            # Loop across parameters, create 1D plot
            for i in range(n_params):

                    ax1 = fig.add_subplot(n_params, n_params, n_params*i+(i+1))
                    ax1.set_xlabel(parameters[i].replace('_','\_'))
            
                    # Get parameter statistics, set bounds to +/- 5sigma
                    m = s['marginals'][i]
            
                    # Histogram posterior
                    x,w,patches = ax1.hist(values[:,i], bins=nbins, edgecolor='black', color='mediumslateblue',\
                            histtype='stepfilled', alpha=0.2,normed=True)
                    ax1.set_ylim(0, x.max())
            
                    # y-value at which to place error-bar (1/20th of vertical extent)
                    ylim = ax1.get_ylim()
                    y = ylim[0] + 0.05*(ylim[1] - ylim[0])

                    # Place error-bar
                    center = m['median']
                    low1, high1 = m['1sigma']
                    ax1.errorbar(x=center, y=y,
                            xerr=numpy.transpose([[center - low1, high1 - center]]), 
                            color='blue', linewidth=2, marker='s')

                    # Axis formatting
                    ax1.set_ylabel("PDF")
                    #ax1.set_xlim(m['5sigma'])
            
                    # Loop across other parameters
                    for j in range(i+1,n_params):

                            # Set up 2D posterior plot
                            ax2 = fig.add_subplot(n_params, n_params, n_params*j+(i+1),sharex=ax1)
                            p.plot_conditional(i, j, cmap = 'Blues',
                                    bins='log',)

                            # Plot error bars for each mode
                            for m in modes:
                                    ax2.errorbar(x=m['mean'][i], y=m['mean'][j], xerr=m['sigma'][i],\
                                            yerr=m['sigma'][j])

                            # Axes labels
                            ax2.set_xlabel(parameters[i].replace('_','\_'))
                            ax2.set_ylabel(parameters[j].replace('_','\_'))

            # Save figure
            plt.savefig(prefix + 'marg.pdf')
            plt.savefig(prefix + 'marg.png')
            plt.close()

    # If $D==1...
    else:

            # Print stuff
            from matplotlib.backends.backend_pdf import PdfPages
            sys.stderr.write('1dimensional only. Set the D environment variable \n')
            sys.stderr.write('to D=2 to force 2d marginal plots.\n')
            pp = PdfPages(prefix + 'marg1d.pdf')

            # Proceed as above, but only plot 1D posterior plots
            for i in range(n_params):

                    plt.figure(figsize=(3, 3))
                    plt.xlabel(parameters[i].replace('_','\_'))
                    plt.locator_params(nbins=5)

                    m = s['marginals'][i]
                    iqr = m['q99%'] - m['q01%']
                    xlim = m['q01%'] - 0.3 * iqr, m['q99%'] + 0.3 * iqr
                    #xlim = m['5sigma']
                    plt.xlim(xlim)

                    oldax = plt.gca()
                    x,w,patches = oldax.hist(values[:,i],
                            bins=numpy.linspace(xlim[0], xlim[1], 20),
                            edgecolor='black', color='grey', histtype='stepfilled', alpha=0.2)
                    oldax.set_ylim(0, x.max())

                    newax = plt.gcf().add_axes(oldax.get_position(), sharex=oldax, frameon=False)
                    p.plot_marginal(i, ls='-', color='blue', linewidth=3)
                    newax.set_ylim(0, 1)

                    ylim = newax.get_ylim()
                    y = ylim[0] + 0.05*(ylim[1] - ylim[0])
                    center = m['median']
                    low1, high1 = m['1sigma']
                    #print center, low1, high1
                    newax.errorbar(x=center, y=y,
                            xerr=numpy.transpose([[center - low1, high1 - center]]),
                            color='blue', linewidth=2, marker='s')
                    oldax.set_yticks([])
                    newax.set_ylabel("Probability")
                    ylim = oldax.get_ylim()
                    newax.set_xlim(xlim)
                    oldax.set_xlim(xlim)
                    plt.savefig(pp, format='pdf', bbox_inches='tight')
                    plt.close()
            pp.close()

if __name__=="__main__":
    if len(sys.argv) != 2:
            sys.stderr.write("""SYNOPSIS: %s <output-root>

            output-root :Where the output of a MultiNest run has been written to.
                            Example: chains/1-
    %s""" % (sys.argv[0], __doc__))
            sys.exit(1)

    # Get parameters
    prefix = sys.argv[1]
    marginals(prefix)
