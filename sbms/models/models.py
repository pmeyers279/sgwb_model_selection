from .power_law import omega_gw_spectrum as pl_omgwf
from ..noise import get_sigma_from_noise
from .broken_power_law import omega_gw_spectrum as bpl_omgwf

# dict of omgw_f models
omgw_f_registry = {'power law' : pl_omgwf,
                   'broken power law' : bpl_omgwf}

# get omega_gw(f) for some set of models
# as specified in a parameter file
def omgwf(params):
    First = True
    for model in params.keys():
        if model=='noise':
            continue
        if model[-2]=='_':
            modname=model[:-2]
        else:
            modname=model
        if First:
            omgwf,f = omgw_f_registry[modname](params[model])
            try:
                str_format=params['noise']['simulation_style']
            except KeyError:
                str_format='file'
            try:
                spectrum_type=params['noise']['spectrum_type']
            except KeyError:
                spectrum_type='psd'
            sig2 =\
                get_sigma_from_noise(int(params['noise']['Tobs']),params['noise']['noise_str'],
                outfs=f, format=str_format, type=spectrum_type)
            First = False
        else:
            omgwf_temp,f = omgw_f_registry[modname](params[model])
            omgwf += omgwf_temp
    return omgwf, sig2, f
