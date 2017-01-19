import numpy as np

def omega_gw_spectrum(omg_ref,f_break=50,alpha1=3,alpha2=-3, flow=20, fhigh=100, df=0.25):
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
    if isinstance(omg_ref, dict):
        # unpack params
        try:
            flow = float(omg_ref['flow'])
        except:
            pass
        try:
            fhigh = float(omg_ref['fhigh'])
        except:
            pass
        try:
            df = float(omg_ref['df'])
        except:
            pass
        try:
            alpha1=float(omg_ref['alpha1'])
        except:
            print 'WARNING: alpha1 is not set...setting it to 3'
            pass
        try:
            alpha2=float(omg_ref['alpha2'])
        except:
            print 'WARNING: alpha2 is not set...setting it to -3'
            pass
        try:
            f_break = float(omg_ref['f_break'])
        except:
            print """WARNING: break frequency "f_break" is not set...setting it
            to 50 Hz"""
            pass
        try:
            omg_ref2=float(omg_ref['omg_ref'])
            omg_ref = omg_ref2
        except:
            raise ValueError('Must specify Omega_alpha for power law model')

    f = np.arange(flow, fhigh+df, df)
    omgwf = np.zeros(f.size)
    slope1 = np.where(f<=f_break)
    slope2 = np.where(f>f_break)
    omgwf[slope1] = omg_ref * (f[slope1]/f_break)**alpha1
    omgwf[slope2] = omg_ref * (f[slope2]/f_break)**alpha2
    return omgwf, f

