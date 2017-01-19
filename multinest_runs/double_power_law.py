from __future__ import division
import matplotlib
matplotlib.use('agg')
import numpy as np
import pymultinest
import json
from sbms.models import omgwf, ptEst_sigma
from sbms.models import power_law
from sbms.io import read_ini
from sbms.orfs import ORF

# get input distribution: two power laws
# each with 1e-9 for omega_alpha
# one with alpha=3, one with alpha=0
params = read_ini('../sbms/tests/data/test_ini.ini')
Y, sig2, f = ptEst_sigma(params)

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
    cube[2] = cube[2]*1e-8 # omega_alpha for model 2
    cube[3] = cube[3]*(5-cube[1]) + cube[1]


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
    omg_alpha1, alpha1, omg_alpha2, alpha2 = cube[0], cube[1], cube[2], cube[3]
    # model
    Y_model, f = power_law.get_data(omg_alpha1, alpha1, fref=100)
    Y_model2, f = power_law.get_data(omg_alpha2, alpha2, fref=100)
    Y_model += Y_model2
    # Y and sigma (global) are our variables
    return -(np.sum((Y_model - Y)**2 / (2*sig2)))


parameters = ["omg_alpha1", "alpha1", "omg_alpha2","alpha2"]
n_params = len(parameters)
pymultinest.run(loglike, prior, n_params, outputfiles_basename='out/test_2_',
        resume=False, verbose=True, n_live_points=2000)
json.dump(parameters, open('out/test_2_params.json','w'))
a = pymultinest.Analyzer(outputfiles_basename='out/test_2_', n_params = n_params)
a_lnZ = a.get_stats()['global evidence']
print
print '************************'
print 'MAIN RESULT: Evidence Z '
print '************************'
print '  log Z for model with 1 line = %.1f' % (a_lnZ / np.log(10))
print

# TODO: implement a model with 2 lines?



