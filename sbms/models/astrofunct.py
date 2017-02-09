# -*- coding: utf-8 -*-
import numpy as np
import scipy.integrate as integrate
import scipy.optimize as sop
import time

def star_form_rate(sfr, z, zvalues = None, sfrs = None):
    if   sfr == 'h': 
        return 0.7*(0.017+0.13*z)/(1+((z/3.3)+0j)**5.3) # I added "+ 0j" to deal with error.  
    elif sfr == 'f': 
        return (0.0103+0.088*z)/(1+((z/2.4)+0j)**2.8)
    elif sfr == 'w': 
        return (0.014+0.11*z)/(1+((z/1.4)+0j)**2.2)
    elif sfr == 'n': 
        # disk
        rhod = lambda t: 0.056*(t/4.5)*np.exp(-t/4.5)
        # Bulk
        rhob = lambda t: 0.198*(t/1.5)*np.exp(-t/1.5)
        
        t = cosmic_time(z)
        return rhod(t) + rhob(t)
        
    elif sfr == 's': 
        r_num = 0.15*(14.0/15)*np.exp(0.6*(z - 5.4))
        r_den = (14.0/15) - 0.6 + 0.6*np.exp((14.0/15)*(z - 5.4))
        return r_num/r_den

    elif sfr == 'b':
        r_num = 0.178*2.37*np.exp(1.80*(z-2.00))
        r_den = 2.37-1.80+1.80*np.exp(2.80*(z-2.00))
        return r_num/r_den

    elif sfr == 'k':
        r_num = 0.146*2.80* np.exp(2.46*(z-1.72))
        r_den = 2.80- 2.46+2.46*np.exp(2.80*(z-1.72))
        return r_num/r_den

    elif sfr == '3':
        r_num = 0.00218*13.81*np.exp(13.36*(z-11.87))
        r_den = 13.81 -13.36+13.36*np.exp(13.81*(z-11.87))
        return r_num/r_den

    elif sfr == 'z':
        return np.interp(z, zvalues, sfrs)

def lookback_time(z, M=0.3, V=0.7):
    # Cosmic (look back) time at the redshift of coalescence z 
    # in unit of Gyr
    return cosmic_time(0, M=0.3, V=0.7) - cosmic_time(z, M=0.3, V=0.7)


def cosmic_time(z, M=0.3, V=0.7):
    # Cosmic (look back) time at the redshift of coalescence z 
    # in unit of Gyr

    Gyr= 3.1536e7*10**9;     # second per gillion year
    H0Mpc = 0.7*10**7;      # cm/(s*Mpc)    -- million parsec 
    Mpc = 3.085e24;         # cm            --
    H0 = H0Mpc/Mpc;         # per second -- hubble constant
    # temp = lambda x: 1.0/(Ez(0.3,0.7,x)*(1+x))
    # time, err = integrate.quad(temp, 0, z)
    # time = 2*np.arcsinh(1/np.sqrt(M*(1+z)**3/V))/3/np.sqrt(V)/(H0*Gyr)
    return 2*np.arcsinh(1/np.sqrt(M*(1+z)**3/V))/3/np.sqrt(V)/(H0*Gyr)
    
def cosmic_time_to_z(t, M=0.3, V=0.7):
    Gyr= 3.1536e7*10**9;     # second per gillion year
    H0Mpc = 0.7*10**7;      # cm/(s*Mpc)    -- million parsec 
    Mpc = 3.085e24;         # cm            --
    H0 = H0Mpc/Mpc;         # per second -- hubble constant
    
    return (1/np.sinh(t*(H0*Gyr)*np.sqrt(V)*3/2)**2*V/M)**(1.0/3) -1

def lookback_time_to_z(t, M=0.3, V=0.7):
    ct = cosmic_time(0, M=0.3, V=0.7) - t
    return cosmic_time_to_z(ct, M=0.3, V=0.7)
    
def hubble_const(z):
    # HUBBLE_CONST Summary of this function goes here
    # Detailed explanation goes here
    return H0*np.sqrt(0.7+0.3*(1+z)**3)

def Ez(OmegaM, OmegaV, z):
    return (OmegaM*(1+z)**3+OmegaV)**0.5  # dimensionless

def BC_rate(z,tmin,sfr, zvalues=None, sfrs=None):
    #   BINARY_COAL_RATE Summary of this function goes here
    #   Detailed explanation goes here

    #Rzf = @(z) star_form_rate(sfr,z)./(1+z)

    ## unnormalized probability distribution of the delay time ______
    # tmin = 0.02
    # tmax = cosmic_time(0)
    zmax = 10.0
    tmax = cosmic_time(zmax)  # tmax is not max, it is smaller than tmin
    # Prob_distri = @(t) 1./t


    ## ________the coalescence rate per comoving volumn__________
    # tsup = (cosmic_time(6)-cosmic_time(z)) #upperbound of integral
   
    tz = cosmic_time(z)
    tsup = (tz - tmax) #upperbound of integral
    P_norm = np.log(13.5/tmin);

    if tsup > tmin:
        #tem = lambda t: star_form_rate(sfr, cosmic_time_to_z(tz-t))/(1+cosmic_time_to_z(tz-t))/t
        if sfr != 'z':
            tem = lambda t: star_form_rate(sfr, cosmic_time_to_z(tz-t))/(1+z)/t
            value, err = integrate.quad(tem, tmin , tsup)
            rate = value/P_norm
        else:
            trange = np.arange(tmin, tsup+0.01, 0.01)
            z_used = []
            for tval in trange:
                z_used.append(cosmic_time_to_z(tz-tval))
            sfr_used = np.interp(z_used, zvalues, np.asarray(sfrs))
            rates = []
            for val in sfr_used:
                rates.append(val/(1+z)/trange[np.where(sfr_used == val)]/P_norm)
            rate =sum(rates) *0.01
    else:
        rate = 0

    
    return rate
