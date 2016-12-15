import numpy as np

def omega_gw_spectrum(omg_alpha, alpha=0, fref=20, flow=20, fhigh=100, df=0.25):
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
    if isinstance(omg_alpha, dict):
        # unpack params
        try:
            flow = omg_alpha['flow']
        except:
            pass
        try:
            fhigh = omg_alpha['fhigh']
        except:
            pass
        try:
            df = omg_alpha['df']
        except:
            pass
        try:
            alpha=omg_alpha['alpha']
        except:
            pass
        try:
            fref = omg_alpha['fref']
        except:
            pass
        try:
            omg_alpha2=float(omg_alpha['omg_alpha'])
            omg_alpha = omg_alpha2
        except:
            raise ValueError('Must specify Omega_alpha for power law model')

    f = np.arange(flow, fhigh+df, df)
    return omg_alpha * (f/fref)**alpha, f
