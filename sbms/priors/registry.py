from ..io_sbms import read_ini
from .priors import *

def get_model_params(ini_file):
    params = read_ini(ini_file)
    new_params = OrderedDict()
    for model in params.keys():
        for key in params[model].keys():
            sp = params[model].split(',')
            new_params[key]['type'] = GeneralPrior



def prior_maker(cube, ndim, nparams):
    counter = 0
    registry = {}
    for model in params.keys():
        for key in params[model].keys():
            sp = params[model].split(',')
            if sp[0]=='lt':
                # if we want to make this thing less than
                # some other value, then 'lt' is first param,
                # and we assume it's uniform otherwise
                cube[counter] =\
                UniformPrior(cube[counter],cube[sp[1]],float(s[2]))
            else:
                cube[counter] = GeneralPrior(float(sp[1]), sp[0], float(sp[2]),
                        float(sp[3]))
            registry[key] = cube[counter]
            counter +=1
