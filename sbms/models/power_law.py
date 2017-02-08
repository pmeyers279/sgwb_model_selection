import numpy as np
from ..noise import get_noise_psd
from ..orfs import ORF

def omega_gw_spectrum(omg_alpha, alpha=0, fref=20, flow=20, fhigh=100, df=0.25,
        frequencies=None):
    """
    Parameters
    ----------
    omg_alpha : `float`
        :math:`\Omega_{\alpha}`, amplitude for power spectrum
    alpha : `float`
        spectral index
    fref : `float`
        reference frequency
    flow : `float`
        low end of frequency spectrum
    fhigh : `float`
        high end of frequency spectrum
    df : `flaot`
        frequency bin width

    Returns
    -------
    omega_gw_f : `numpy.ndarray`
        :math:`\Omega_{gw}(f)`
    f : `numpy.ndarray`
        frequency array
    """
    f = np.arange(float(flow), float(fhigh)+float(df), float(df))
    if frequencies is not None:
        f = frequencies
    return 10**omg_alpha * (f/float(fref))**alpha, f

def unpack_dict(param_dict):
    kwargs = {}
    args = []
    try:
        args.append(float(param_dict.pop('omg_alpha')))
    except:
        raise ValueError('Must specify Omega_alpha for power law model')
    for key in param_dict.keys():
        if key != 'frequencies':
            param_dict[key] = float(param_dict[key])
    return args, param_dict
