import numpy as np
import scipy as sp
from constants import *
from astrofunct import *
import math
from operator import *

def omega_gw_spectrum(M=10.0, lamb= 1e-6,eta=0.1, zi=0.0, flow = 20, fhigh=100, df=0.25, sfr ='h', frequencies = None):
    """
    Parameters
    ---------
    M: 'float'
       chirp mass, in solar masses [2.5, 20]
    lamb: 'float'
       mass fraction, in solar masses [0, .1]
    eta: 'float'
       mass ratio (0, 0.25]
    zi: 'float'
       spin parameter [-0.85, 0.85]
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

    #this follows the analysis in https://arxiv.org/abs/1104.3565

    f = np.arange(float(flow),float(fhigh)+float(df),float(df))
    if frequencies is not None:
        f = frequencies
    omgwf = np.zeros(f.size)
    
    Kb = (G*np.pi)**(2.0/3.0)*(M*Msolar)**(5.0/3)/3

    Const = (8*np.pi*G)/(3*c**2*H0**3)*lamb*Kb/yr/Mpc**3

    zmin = 0.0
    zmax = 10.0 #already accounted for in BC_Rate function
    tmin = 0.05
    tmax = 13.5
    fgmin = 0.0
    #fgmax = v3 #maximum frequency

    v_ = lambda f: (np.pi*Mtot*f*LAL_MTSUN_SI)**(1.0/3.0)

    #calulating total mass of black hole
    M = float(M)
    root = (0.25-eta)**(1.0/2.0)#eta cannot be more than 0.25
    fraction = (0.5 + root) / (0.5 - root)
    invfraction = 1.0 / fraction
    m2 = (M*(1+fraction)**0.2)/fraction**0.6
    m1 = (M*(1+invfraction)**0.2)/invfraction**0.6
    Mtot = m1 + m2
    
    v_ = lambda f: (np.pi*Mtot*f*LAL_MTSUN_SI)**(1.0/3.0)

    v1 = lambda z: (404*(0.66389*eta**2 - 0.10321*eta + 0.10979)/0.125481)*(20.0/Mtot)/(1+z);
    v2 = lambda z: (807*(1.3278*eta**2 - 0.20642*eta + 0.21957)/0.250953)*(20.0/Mtot)/(1+z);
    v3 = lambda z:(1153*(1.7086*eta**2 - 0.26592*eta + 0.28236)/0.322668)*(20.0/Mtot)/(1+z);
    sigma = lambda z:(237 * (1.1383*eta**2 - 0.177*eta + 0.046834)/0.0737278)*(20.0/Mtot)/(1+z);

    mod1 = lambda v_:(1 + (-323.0/224.0+ 451.0*eta/168.0)*v_**2 + zi*(27.0/8.0 - 11.0*eta/6.0)*v_**3)**2
    mod2 = lambda v_: (1+ (1.4545*zi - 1.8897)*v_ + (-1.8153*zi + 1.6557)*v_**2)**2

    w1 =lambda v1, mod1, mod2:(1 / v1) * (mod1 / mod2)
    w2 = lambda v1, v2, mod1:(1 / v1)*(1/v2**(4.0/3.0))*mod1

    def getv(z, f):
        v = 0
        v__ = v_(f)
        v1_ = v1(z)
        v2_ = v2(z)
        v3_ = v3(z)
        sigma_ = sigma(z)
        mod1_ = mod1(v__)
        mod2_ = mod2(v__)
        w1_ = w1(v1_, mod1_, mod2_)
        w2_ = w2(v1_, v2_, mod1_)
        if f >=  fgmin and f <= v3_:
            if 0 < f and f < v1_:
                v = f**(-1.0/3.0)*mod1_
            elif v1_ <= f and f <= v2_:
                v = f**(2.0/3.0)*w1_*mod2_
            else: # v2 < f[x]  and f[x] < v3:
                v = f**2.0*w2_/(1+ 4*(f-v2_)**2/sigma_**2)**2
        return v

    if sfr != "z":
         
        integrand = lambda z,f: getv(z, f)/Ez(OmegaM, OmegaV, z)/(1+z)**(1.0/3.0)
        zrange = np.arange(zmin, zmax+0.1, 0.1)
        rates = map(lambda zval: BC_rate(zval,tmin,sfr) ,zrange)
        zrange2 = np.arange(zmin, zmax+0.01, 0.01)
        vfunc = np.vectorize(lambda zval,x: np.interp(zval, zrange, rates)*integrand(zval, x))
        vfunc2 = np.vectorize(lambda x: sum(vfunc(zrange2, x))*0.01*Const*x)
        omgwf = vfunc2(f);

    else:
       
        lines = []
        zvalues = []
        sfrs = []
        fil = open('EfficiencyZ_2.txt')
        for line in fil.readlines():
            lines.append(line.strip().split('\t'))
        for element in lines:
            zvalues.append(float(element[0]))
            sfrs.append(star_form_rate('k', float(element[0]))* float(element[1]))
        fil.close()
        rate = lambda z:BC_rate(z, tmin, sfr, zvalues, sfrs)
        integrand_part = lambda z,f: getv(z,f)/Ez(OmegaM, OmegaV, z)/(1+z)**(1.0/3.0)
        zrange = np.arange(zmin, zmax+0.01, 0.01)
        rates = []
        for zval in zrange:
            rates.append(rate(zval)) 
        for x in range(len(f)):
            vals = []
            for index,zval in enumerate(zrange):
                vals.append(integrand_part(zval, f[x])*rates[index])
            val = np.sum(vals) * 0.01
            omgwf[x] = Const*f[x]* val
    
    if sfr != "b":
        return omgwf, f
    else:
        zrange = np.arange(zmin, zmax+0.1, 0.1)
        rates = map(lambda zval: BC_rate(zval,tmin,'3') ,zrange)
        zrange2 = np.arange(zmin, zmax+0.01, 0.01)
        vfunc = np.vectorize(lambda zval,x: np.interp(zval, zrange, rates)*integrand(zval, x))
        vfunc2 = np.vectorize(lambda x: sum(vfunc(zrange2, x))*0.01*Const*x)
        omega1 = vfunc2(f);
        omgwf = map(add, omgwf, omega1)

        return omgwf,f

def unpack_dict(param_dict):
    kwargs = {}
    args = []
    for key in param_dict.keys():
        if key != 'frequencies':
            param_dict[key] = float(param_dict[key])
    return args, param_dict
