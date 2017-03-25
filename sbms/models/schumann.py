import numpy as np
from astropy import cosmology, units
from ..orfs import ORF

def omega_gw_spectrum(kappa_amp=2, beta_amp=2.67, kappa_ang=0.25, beta_ang=1.74,
                      M_f=None, flow=20, fhigh=100, df=0.25, frequencies=None):
    """
    Parameters
    ----------
    kappa_amp : `float`
        length coupling transfer function amplitude
        parameter.
    beta_amp : `float`
        length coupling transfer function power law
        index.
    kappa_ang : `float`
        angular coupling transfer function amplitude
    beta_ang : `float`
        angular coupling transfer function power law
        index.
    M_f : `numpy.ndarray`
        Correlated magnetic field power spectrum
        in units of :math:`pT^2`
    flow : `float`
        low frequency
    fhigh : `float`
        high frequency
    df : `float`
        bin width
    frequencies : `numpy.ndarray`
        list of frequencies to use for analysis.
        if specified this supersedes flow, fhigh, df
        specifications

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
    # if nothing is supplied, then use 1pT^2
    if M_f is None:
        M_f = np.ones(f.size)
    Omgw = np.zeros(f.size)
    H_amp = kappa_amp**2 * M_f * (1e-46) * (f/10)**(-2*beta_amp)
    H_ang = kappa_ang**2 * M_f * (1e-46) * (f/10)**(-2*beta_ang)
    #T,V,S = ORF(f)
    H0_s = cosmology.WMAP9.H0.to(units.s**-1).value
    Omgw += H_amp * (f)**(3) * (2 * np.pi**2) / (3 * H0_s**2)
    Omgw += H_ang * (f)**(3) * (2 * np.pi**2) / (3 * H0_s**2)
    return Omgw, f

def unpack_dict(param_dict):
    return param_dict