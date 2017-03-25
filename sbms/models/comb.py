from __future__ import division
import numpy as np

def omega_gw_spectrum(comb_height=0, offset=0.5, separation=0.5, flow=20, fhigh=100,
        df=0.25, frequencies=None, decimals=5):
    """
    Parameters
    ----------
    omg_ref : `float`
        :math:`\Omega_{\alpha}`, amplitude for power spectrum at the cross-over
        from one power law to the other
    f_break : `float`
        frequency at which cross-over occurs
    alpha1 : `float`
        power law index for first section
    alpha2 : `float`
        power law index for second section
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
    # fix these
    offset_old = offset
    separation_old = separation
    if frequencies is not None:
        f = frequencies
        df = np.min(f[1:] - f[:-1])
    else:
        f = np.arange(float(flow), float(fhigh)+float(df), float(df))
    omgwf = np.zeros(f.size)
    offset = np.round(np.round(offset/df) * df, decimals)
    separation = np.round(np.round(offset/df) * df, decimals)
    try:
        mask_idxs = np.arange(offset/df, omgwf.size-separation/df, separation/df,
             dtype=int)
    except:
        raise ValueError('this is failing...offset %4.5f, df=%4.5f,\
        separation=%4.5f, old_offset=%4.5f, old_separation=%4.5f' % (offset,
            df, separation, offset_old, separation_old))
    omgwf[mask_idxs] = 10**comb_height
    return omgwf, f

def unpack_dict(param_dict):
    kwargs = {}
    args = []
    try:
        args.append(float(param_dict.pop('comb_height')))
    except:
        raise ValueError('Must specify Omega_alpha for power law model')
    for key in param_dict.keys():
        if key != 'frequencies':
            param_dict[key] = float(param_dict[key])
    return args, param_dict

