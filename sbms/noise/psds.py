import os
import numpy as np
from scipy.interpolate import interp1d
from ..orfs import ORF

H0 = 2.2685e-18

def get_sigma_from_noise(Tobs, noise_curve_string, outfs=None,
        format='lal',type='asd'):


    Tobs *= 24 * 365 * 3600
    norm = (10.*np.pi**2.)/(3.*H0**2.)
    f, psd = get_noise_psd(noise_curve_string, outfs=outfs, format=format, type=type)

    torf, vorf, sorf = ORF(f)
    df = f[1] - f[0]
    sigma2s = np.power(norm,2.)*(1./(2.*Tobs*df*np.abs(torf)**2))*np.power(f,6.)*np.power(psd,2.)
    return sigma2s

def get_noise_psd(noise_curve_string, outfs=None,
        format='lal', type='asd'):
    """
    get noise curve using lalsim-detector-noise.
    This function makes a temporary file to store output and then
    removes it.

    Parameters
    ----------
    noise_curve : `str`
        flag for lalsim-detector-noise or link to file
    out_fs : `numpy.ndarray`
        list of frequencies to output noise curve with
    format : `str`, options=['lal','file']
        what is the format of the noise_curve_str?
        'lal' means its an argument that can be passed to
        'lalsim-detector-noise'. 'file' means it's a file
        we want to load
    type : `str`, options=['asd','psd']
        What are we getting from the file or from lalsim-detector?
        Is it an asd or a psd? If it's an asd we know to square the result
        otherwise we don't square it.

    Returns
    -------
    freqs : `numpy.ndarray` (if 'psd')
        frequencies
    noise : `numpy.ndarray`
        noise values
    """
    if format=='lal':
        # set up arguments
        args = ''
        args += ' -P '
        args += noise_curve_string
        args += ' -f %f' % flow
        # run lalsim-detector-noise
        os.system('lalsim-detector-noise %s > tmp.txt' % args)
        fname='tmp.txt'
    elif format=='file':
        fname = noise_curve_string
        power = 1
    if type=='asd':
        power = 2
    elif type=='psd':
        power = 1
    # load in psd
    data = np.loadtxt(fname, comments='#')
    freqs = data[:,0]
    noise = data[:,1]
    if format=='lal':
        # remove temp noise curve
        os.remove('tmp.txt')
    if outfs is None:
        return freqs, noise**power
    else:
        #interpolate
        newvals = interp1d(freqs, noise)
        # clean up
        return outfs, newvals(outfs)**power
