from __future__ import division
import matplotlib
matplotlib.use('agg')
import numpy as np
import pymultinest
import json
from sbms.models import omgwf
from sbms.models import power_law
from sbms.io import read_ini

# get input distribution: two power laws
# each with 1e-9 for omega_alpha
# one with alpha=3, one with alpha=0
params = read_ini('../sbms/tests/data/test_ini.ini')
omega_gw, f = omgwf(params)

def prior(cube, ndim, nparams):
    """
    set up prior probability distributions

    Parameters
    ----------
    cube : array
        variables we want to set prior probabilities on
    ndim : array
        number of dimensions
    nparams : int
        number of parameters

    Returns
    -------
    none

    """
    cube[0] = cube[0]*1e-8 # omega_alpha for model 1
    cube[1] = cube[1]*10 - 5 # alpha for model 1


def loglike(cube, ndim, nparams):
    """
    log likelihood

    Parameters
    ----------
    cube : array
        array of parameters
    ndim : int
        number of dimensions
    nparams : int
        number of parameters

    Returns
    -------
    likelihood : int
        value of the likelihood given a realization
        of values chosen from the priors for each
        variable
    """
    omg_alpha1, alpha1 = cube[0], cube[1]
    omgw_model,f = power_law.omega_gw_spectrum(omg_alpha1, alpha1, fref=100)
    sigma = omega_gw / 10
    return -(np.sum((omgw_model - omega_gw)**2 / (2*sigma**2)))


parameters = ["omg_alpha1", "alpha1"]
n_params = len(parameters)
pymultinest.run(loglike, prior, n_params, outputfiles_basename='out/test_1_',
        resume=False, verbose=True, n_live_points=2000)
json.dump(parameters, open('out/test_1_params.json','w'))
a = pymultinest.Analyzer(outputfiles_basename='out/test_1_', n_params = n_params)
a_lnZ = a.get_stats()['global evidence']
print
print '************************'
print 'MAIN RESULT: Evidence Z '
print '************************'
print '  log Z for model with 1 line = %.1f' % (a_lnZ / np.log(10))
print

# TODO: implement a model with 2 lines?



