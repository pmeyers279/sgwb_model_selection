from .power_law import omega_gw_spectrum as pl_omgwf

# dict of omgw_f models
omgw_f_registry = {'power law' : pl_omgwf}

# get omega_gw(f) for some set of models
# as specified in a parameter file
def omgwf(params):
    First = True
    for model in params.keys():
        if model[-2]=='_':
            modname=model[:-2]
        else:
            modname=model
        if First:
            omgwf,f = omgw_f_registry[modname](params[model])
            First = False
        else:
            print 'HI'
            omgwf_temp,f = omgw_f_registry[modname](params[model])
            omgwf += omgwf_temp
    return omgwf, f
