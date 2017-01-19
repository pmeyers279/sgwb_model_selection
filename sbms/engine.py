#! /usr/bin/python
from __future__ import division
import matplotlib
import time
matplotlib.use('agg')
import os
import matplotlib.pyplot as plt
import numpy as np
import pymultinest
import json
from sbms.models import omgwf, omgw_f_registry
from sbms.models import power_law
from sbms.io_sbms import read_ini
from sbms.orfs import ORF
from sbms.priors.priors import GeneralPrior
import optparse

def parse_command_line():
    """
    parse command line
    """
    parser = optparse.OptionParser()
    parser.add_option("--injection-file", "-i",
        help="injection param file", default=None,
        dest="inj_file", type=str)
    parser.add_option("--recovery-file", "-r",
        help="recovery param file", default=None,
        dest="recovery_file", type=str)
    parser.add_option("--output-prefix", "-o",
        help="output prefix. path needs to exist right now.",
        default='./multinest_files_',
        dest="output_pref", type=str)
    params, args = parser.parse_args()
    return params

def engine(injection_file, recovery_file, output_prefix='./multinest_files_'):
    # get input distribution: two power laws
    # each with 1e-9 for omega_alpha
    # one with alpha=3, one with alpha=0
    params = read_ini(injection_file)
    params2 = read_ini(recovery_file)
    Omgw, sig2, f = omgwf(params)
    Omgw += np.random.randn(Omgw.size) * np.sqrt(sig2)

    plt.figure()
    plt.plot(f, Omgw, label=r'Input $Y(f)$')
    plt.plot(f,np.sqrt(sig2), label=r'Input $\sigma_Y(f)$')
    plt.yscale('log')
    plt.legend()
    plt.savefig(output_prefix + 'input_data')

    parameters = []
    nparams = 0
    for key in params2.keys():
        nparams = params2[key]['nparams']
        for ii in range(int(nparams)):
            parameters.append(params2[key]['param'+str(ii+1)].split(',')[0])
    print parameters

    print 'OK...we are recovering with...'
    for key in params2.keys():
        print '=================='
        print 'Model ' + key
        print 'Model type: ' + params2[key]['model_type']
        for ii in range(int(params2[key]['nparams'])):
            sp = params2[key]['param'+str(ii+1)].split(',')
            print 'Parameter %d: %s' % (ii+1, sp[0])
            print '\t%s Prior' % (sp[1])
            print '\t%s Prior Params: %f %f' % (sp[1],float(sp[2]),float(sp[3]))
    print '=================='
    print 'If you see anything wrong with this, stop it now!'
    print 'Starting in...'
    for ii in range(5):
        print '%d...' % (5-ii)
        time.sleep(1)


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
        # construct prior from recovery file
        counter = 0
        for key in params2.keys():
            nparams_tmp = int(params2[key]['nparams'])
            for ii in range(nparams_tmp):
                # sp = [name, prior type, x1, x2]
                sp =\
                    params2[key]['param'+str(ii+1)].split(',')
                cube[counter] = GeneralPrior(cube[counter], sp[1], float(sp[2]),
                        float(sp[3]))
                counter += 1

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
        # construct our model
        Y_model = np.zeros(Omgw.size)
        counter = 0
        for key in params2.keys():
            model = params2[key]['model_type']
            nparams = int(params2[key]['nparams'])
            args = [cube[i] for i in range(counter,counter+nparams)]
            Y_temp,f = omgw_f_registry[model.replace('_',' ')](*args)
            Y_model += Y_temp
            counter += nparams
        return -(np.sum((Y_model - Omgw)**2 / (4*sig2)))

    #parameters = ["omg_alpha1", "alpha1", "omg_alpha2","alpha2"]
    n_params = len(parameters)
    pymultinest.run(loglike, prior, n_params,
            outputfiles_basename=output_prefix,
            resume=False, verbose=True, n_live_points=2000)
    json.dump(parameters, open(output_prefix + 'params.json','w'))
    a = pymultinest.Analyzer(outputfiles_basename=output_prefix, n_params = n_params)
    a_lnZ = a.get_stats()['global evidence']
    print
    print '************************'
    print 'MAIN RESULT: Evidence Z '
    print '************************'
    print '  log Z for model with 1 line = %.1f' % (a_lnZ / np.log(10))
    print

    ## TODO: implement a model with 2 lines?

if __name__=="__main__":
    from sbms.marginals import marginals
    params = parse_command_line()
    # run multinest first
    if not os.path.exists(params.inj_file):
        raise IOError('Your injection file doesnt exist!')
    if not os.path.exists(params.recovery_file):
        raise IOError('Your recovery file doesnt exist!')
    if not os.path.exists(os.path.dirname(params.output_pref)):
        raise IOError("""The directory where you want to save things doesnt
                exist yet!""")
    engine(params.inj_file, params.recovery_file, params.output_pref)
    print '=============================='
    print 'DONE RUNNING MULTINEST        '
    print ''
    print 'RUNNING MULTINEST MARGINALS   '
    print '=============================='
    # run multinest marginals as well
    marginals(params.output_pref)
