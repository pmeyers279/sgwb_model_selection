import numpy as np
import scipy as sp
from constants import *
from astrofunct import *
import math
from operator import *
import time

def omega_gw_spectrum(M=10.0, lamb= 1e-6,flow = 20, fhigh=100, df=1.0, sfr ='h', frequencies = None):
    """
    Parameters
    ---------
    M: 'float'
       chirp mass, in solar masses [1, 2.5]
    lamb: 'float'
       mass fraction, in solar masses [0, .1]
    flow: 'float'
       low end of frequency spectrum
    fhigh: 'float'
       high end of frequency spectrum
    df: 'float'
       frequency bin width
    sfr: 'string'
       star formation rate with options Hopkins and Beacom 2006('h'), Fardal 2007 ('f'), Wilken 2008 ('w'), Nagamine 2006 ('n'), Springel and Hernquist 2003 ('s'), Behroozi et. al 2013 ('b'), Kistler et al 2013 ('k'), Low Metallicity ('z')
    
     Returns                                                            
    -------                                                             
    omega_gw_f : `numpy.ndarray`                                        
        :math:`\Omega_{gw}(f)`                                          
    f : `numpy.ndarray`                                                 
        frequency array                                                
    """

    start = time.time()

    mm = M*Msolar # convert unit to gram 
    Kb = (G*np.pi)**(2.0/3.)*mm**(5.0/3.)/3.  # cm**2*g/s**(4/3)

    # Cutoffs of recieved frequence range 
    # per second
    # here 2**(6/5)*m is the total mass of the two NS
    
    f = np.arange(float(flow),float(fhigh)+float(df),float(df))
    fgmin = 10
    fgmax = c**3/(6.**1.5*np.pi*G*(2.**(6.0/5.)*mm))

    
    for i in range(0, len(f)- 1):
        if f[i]<fgmin and f[i+1]>fgmin:
            f[i]= fgmin
        elif f[i]>fgmin and f[i]<fgmax and f[i+1]>fgmax:
            continue;
        elif f[i]<fgmax and f[i+1]>fgmax:
            f[i] = fgmax
        elif f[i]>fgmax:
            f[i]=fgmax
            break;
        if len(f) ==2:
            f[0] = fgmin
            f[1] = fgmax
    
    # dimensionless 
    # max of the redshift
    zmax = 6.  

    # normalizer of delay time distribution
    tmin = 0.02

    Const = (8.*np.pi*G)/(3.*c**2*H0**3)*lamb*Kb/yr/Mpc**3
    #omega = np.array(len(f)*[0], float)
    omgwf = np.zeros((len(f),1))

    if sfr != "z":
        
        integrand = np.vectorize(lambda z:1./Ez(OmegaM, OmegaV, z)/(1+z)**(1.0/3.0))
        def zrange2(fval):
            if fval >= fgmin and fval <= fgmax: 
                if fval != 0 and fgmax/fval - 1 < 6 : 
                    zsup = fgmax/fval - 1
                else:
                    zsup = 6.
            zrange = np.arange(0, zsup, 0.01)
            return zrange
        zrange = np.arange(0, 6., 0.1)
        rates = map(lambda zval: BC_rate(zval,tmin,sfr) ,zrange)
        vfunc = np.vectorize(lambda zval: np.interp(zval, zrange, rates)*integrand(zval))
        vfunc2 = np.vectorize(lambda x: (vfunc(zrange2(x)[0])+vfunc(zrange2(x)[-1])+sum(2.0*vfunc(zrange2(x)[1:-1])))/2.0*0.01*Const*x**(2.0/3.0)) 
        omgwf = vfunc2(f);
        
    else:
         #Low Mettalicity sfr- This needs to be revised, Following code is too slow...  
        lines = []
        zvalues = []
        sfrs = []
        fil = open('/home/fitzaxen/Workspace/spectrum/EfficiencyZ_2.txt')
        for line in fil.readlines():
            lines.append(line.strip().split('\t'))
        for element in lines:
            zvalues.append(float(element[0]))
            sfrs.append(star_form_rate('k', float(element[0]))* float(element[1]))
        fil.close()
        rate = lambda z:BC_rate(z, tmin,sfr, zvalues, sfrs)
        zrange = np.arange(0, 6+0.01, 0.01)
        rates = []
        for zval in zrange:
            rates.append(rate(zval))
        integrand_part = lambda z:1/Ez(OmegaM,OmegaV,z)/(1+z)**(1.0/3.)
        omega = np.array(len(f)*[0], float)
        for x in range(len(f)):
            if  ((f[x] >=  fgmin) and (f[x]<= fgmax)):         
                if f[x] != 0 and fgmax/f[x] - 1 < 6 :
                    zmax = fgmax/f[x] - 1
                else:
                    zmax = 6
                vals = []
                zrange = np.arange(0, zmax+0.01, 0.01)
                for index,zval in enumerate(zrange):
                    vals.append(integrand_part(zval)*rates[index])
                val = np.sum(vals) * 0.01
                omgwf[x] = Const*f[x]**(2.0/3.)*val 
        
        
    if sfr != "b":
        #print(time.time()-start)
        return omgwf, f
    else:

        rates = map(lambda zval: BC_rate(zval,tmin,'3') ,zrange)
        vfunc = np.vectorize(lambda zval: np.interp(zval, zrange, rates)*integrand(zval))
        vfunc2 = np.vectorize(lambda x: (vfunc(zrange2(x)[0])+vfunc(zrange2(x)[-1])+sum(2.0*vfunc(zrange2(x)[1:-1])))/2.0*0.01*Const*x**(2.0/3.0))
        omega1 = vfunc2(f);
        omgwf = map(add, omgwf, omega1)
        #print(time.time()-start)
        return omgwf,f

def unpack_dict(param_dict):
    kwargs = {}
    args = []
    for key in param_dict.keys():
        if key != 'frequencies':
            param_dict[key] = float(param_dict[key])
    return args, param_dict
