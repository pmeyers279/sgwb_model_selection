import numpy as np

def omega_gw_spectrum(log_omg_1, fref=20, flow=20, fhigh=100, df=0.25):
    f = np.arange(flow, fhigh+df, df)
    return (10**log_omg_1)*(f/fref), f

def unpack_dict(param_dict):
    kwargs = {}
    args = []
    try:
        args.append(float(param_dict.pop('log_omg_1')))
    except:
        raise ValueError('Must specify Omega_alpha for power law model')
    for key in param_dict.keys():
        if key != 'frequencies':
            param_dict[key] = float(param_dict[key])
    return args, param_dict
