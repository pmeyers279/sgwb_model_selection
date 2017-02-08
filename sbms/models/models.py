from .power_law import omega_gw_spectrum as pl_omgwf
from .power_law import unpack_dict as pl_params
from .broken_power_law import omega_gw_spectrum as bpl_omgwf
from .broken_power_law import unpack_dict as bpl_params
from .simple_line import omega_gw_spectrum as sl_omgwf
from .simple_line import unpack_dict as sl_params
from ..noise import get_sigma_from_noise


# dict of omgw_f models
omgw_f_registry = {'power law' : pl_omgwf,
                   'broken power law' : bpl_omgwf,
                   'simple line' : sl_omgwf}
param_registry = {'power law' : pl_params,
                  'broken power law' : bpl_params,
                  'simple line' : sl_params}

# get omega_gw(f) for some set of models
# as specified in a parameter file
def omgwf(params):
    First = True
    for model in params.keys():
        if model=='noise':
            noise_dict = params.pop('noise')
            continue
        if model[-2]=='_':
            modname=model[:-2]
        else:
            modname=model
        if First:
            a,k = param_registry[modname](params[model])
            omgwf,f = omgw_f_registry[modname](*a, **k)
            try:
                str_format=noise_dict['simulation_style']
            except KeyError:
                str_format='file'
            try:
                spectrum_type=noise_dict['spectrum_type']
            except KeyError:
                spectrum_type='psd'
            sig2 =\
                get_sigma_from_noise(int(noise_dict['Tobs']),noise_dict['noise_str'],
                outfs=f, format=str_format, type=spectrum_type)
            First = False
        else:
            a,k = param_registry[modname](params[model])
            omgwf_temp,f = omgw_f_registry[modname](*a, **k)
            omgwf += omgwf_temp
    return omgwf, sig2, f
