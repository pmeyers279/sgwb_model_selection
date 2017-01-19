import os
import numpy as np
from scipy.interpolate import interp1d

H0 = 2.2685e-18

def get_sigma_from_noise(Tobs, noise_curve_string, outfs=None):
    Tobs *= 24 * 365 * 3600
    norm = (10.*np.pi**2.)/(3.*H0**2.)
    f, psd = get_noise_psd(noise_curve_string, outfs=outfs)
    df = f[1] - f[0]
    sigma2s = np.power(norm,2.)*(1./(2.*Tobs*df))*np.power(f,6.)*np.power(psd,2.)
    return sigma2s

def get_noise_psd(noise_curve_string, outfs=None, flow=20, fhigh=200):
    """
    get noise curve using lalsim-detector-noise.
    This function makes a temporary file to store output and then
    removes it.

    Parameters
    ----------
    noise_curve : `str`
        flag for lalsim-detector-noise
    outputtype : `str`
        output type. 'psd' or 'ts'

    Returns
    -------
    freqs : `numpy.ndarray` (if 'psd')
        frequencies
    times : `numpy.ndarray` (if 'ts')
        times
    noise : `numpy.ndarray`
        noise values
    """
    # set up arguments
    args = ''
    args += ' -P '
    args += noise_curve_string
    args += ' -f %f' % flow
    # run lalsim-detector-noise
    os.system('lalsim-detector-noise %s > tmp.txt' % args)
    # load in psd
    data = np.loadtxt('tmp.txt', comments='#')
    freqs = data[:,0]
    noise = data[:,1]
    os.remove('tmp.txt')
    if outfs is None:
        return freqs, noise
    else:
        #interpolate
        newvals = interp1d(freqs, noise)
        # clean up
        return outfs, newvals(outfs)**2
