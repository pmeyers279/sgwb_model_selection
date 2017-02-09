# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp
from numpy import interp
from scipy import integrate
from scipy.integrate import odeint
from scipy.integrate import ode
from scipy.integrate import quad
from scipy.integrate import dblquad
from scipy.integrate import nquad
from scipy.optimize import fminbound
from scipy.optimize import fsolve
from scipy.optimize import minimize_scalar
from scipy import interpolate
from scipy import *
from scipy.special import *
from constants import *
from astrofunct import *
from ompclib_demo import *
import os
import math
import subprocess
import time
##from subprocess import call
##from subprocess import check_output
from subprocess import *

def hanford(f):
    fgmin = 460
    fgmax = 1000
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    return [ (x >=  fgmin and x <= fgmax)*7.7e-4*(x/900)**3.0 for x in g]

def pulsar_2015(f):
    fgmin = 2.0e-9
    fgmax =  5.0e-9
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    return [ (x >=  fgmin and x <= fgmax)*(4.2e-10)/(0.678)**2.0 for x in g]
    

def pulsar_2013(f):
    fgmin = 2.0e-9
    fgmax =  5.0e-9
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    return [ (x >=  fgmin and x <= fgmax)*5.3e-9 for x in g] 

def LIGO_S6_1(f):
    fgmin = 41.5
    fgmax = 169.25
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    return [ (x >=  fgmin and x  <= fgmax)*5.6e-6 for x in g]

def LIGO_S6_2(f):
    fgmin = 170
    fgmax = 600
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    return [ (x >=  fgmin and x < fgmax)*1.8e-4 for x in g]

def LIGO_S6_3(f):
    fgmin = 600
    fgmax = 1000
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    return [ (x >=  fgmin and x < fgmax)*0.14*(x/900)**3.0 for x in g]

def LIGO_S6_4(f):
    fgmin =1000
    fgmax =1726
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    return [ (x >=  fgmin and x < fgmax)*1.0*(x/1300)**3.0 for x in g]

def LIGO_S5(f):
    fgmin = 41.5
    fgmax = 169.25
    #if f[0]>fgmax or f[len(f)-1]<fgmin:
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax

    #print(' final f is ')
    #print f
    #omega = np.zeros((len(f),1))
   # for x in range(fgmin,fgmax):
	#    omega[x,0]= 6.9e-6
    #return omega[:,0]
    return [ (x >=  fgmin and x <= fgmax)*6.9e-6 for x in g] 

def LIGO_S4(f):
    fgmin = 51
    fgmax = 150
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	   # print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	  # print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	   # print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    #print(' final f is ')
    #print f
    return [ (x >=  fgmin and x <= fgmax)*6.5e-5 for x in g] 

def LIGO_S3(f):
    fgmin = 69
    fgmax = 156
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    #print(' final f is ')
    #print f
    return [ (x >=  fgmin and x <= fgmax)*8.4e-4 for x in g]

def Advan_LIGO_Mid(f):
    fgmin = 10
    fgmax = 200
    for i in range(0, len(f)- 1):
	if f[i]<fgmin and f[i+1]>fgmin:
	    #print('MinHit')
	    f[i]= fgmin
	elif f[i]>fgmin and f[i]<fgmax and f[i+1]>fgmax:
	    continue;
	elif f[i]<fgmax and f[i+1]>fgmax:
	    #print('MaxHit')
	    f[i] = fgmax
	elif f[i]>fgmax:
	    #print('BreakHit')
	    f[i]=fgmax
	    break;
	if len(f) ==2:
	    f[0] = fgmin
	    f[1] = fgmax
    # print(' final f is ')
    # print f
    return [ (x >=  fgmin and x <= fgmax)*9e-9 for x in f]

def Advan_LIGO_Final(f):
    fgmin = 10
    fgmax = 200
    for i in range(0, len(f)- 1):
	if f[i]<fgmin and f[i+1]>fgmin:
	    #print('MinHit')
	    f[i]= fgmin
	elif f[i]>fgmin and f[i]<fgmax and f[i+1]>fgmax:
	    continue;
	elif f[i]<fgmax and f[i+1]>fgmax:
	    #print('MaxHit')
	    f[i] = fgmax
	elif f[i]>fgmax:
	    #print('BreakHit')
	    f[i]=fgmax
	    break;
	if len(f) ==2:
	    f[0] = fgmin
	    f[1] = fgmax
    # print(' final f is ')
    # print f
    return [ (x >=  fgmin and x <= fgmax)*1.1e-9 for x in f]

def Advan_LIGOp(f):
    fgmin = 10
    fgmax = 200
    for i in range(0, len(f)- 1):
	if f[i]<fgmin and f[i+1]>fgmin:
	    #print('MinHit')
	    f[i]= fgmin
	elif f[i]>fgmin and f[i]<fgmax and f[i+1]>fgmax:
	    continue;
	elif f[i]<fgmax and f[i+1]>fgmax:
	    #print('MaxHit')
	    f[i] = fgmax
	elif f[i]>fgmax:
	    #print('BreakHit')
	    f[i]=fgmax
	    break;
	if len(f) ==2:
	    f[0] = fgmin
	    f[1] = fgmax
    # print(' final f is ')
    # print f
    return [ (x >=  fgmin and x <= fgmax)*3.7e-10 for x in f]

def Advan_LIGOpp(f):
    fgmin = 10
    fgmax = 200
    for i in range(0, len(f)- 1):
	if f[i]<fgmin and f[i+1]>fgmin:
	    #print('MinHit')
	    f[i]= fgmin
	elif f[i]>fgmin and f[i]<fgmax and f[i+1]>fgmax:
	    continue;
	elif f[i]<fgmax and f[i+1]>fgmax:
	    #print('MaxHit')
	    f[i] = fgmax
	elif f[i]>fgmax:
	    #print('BreakHit')
	    f[i]=fgmax
	    break;
	if len(f) ==2:
	    f[0] = fgmin
	    f[1] = fgmax
    # print(' final f is ')
    # print f
    return [ (x >=  fgmin and x <= fgmax)*1.5e-10 for x in f]

def Voyager(f):
    fgmin = 10
    fgmax = 200
    for i in range(0, len(f)- 1):
	if f[i]<fgmin and f[i+1]>fgmin:
	    #print('MinHit')
	    f[i]= fgmin
	elif f[i]>fgmin and f[i]<fgmax and f[i+1]>fgmax:
	    continue;
	elif f[i]<fgmax and f[i+1]>fgmax:
	    #print('MaxHit')
	    f[i] = fgmax
	elif f[i]>fgmax:
	    #print('BreakHit')
	    f[i]=fgmax
	    break;
	if len(f) ==2:
	    f[0] = fgmin
	    f[1] = fgmax
    # print(' final f is ')
    # print f
    return [ (x >=  fgmin and x <= fgmax)*1.7e-11 for x in f]


def Cosmic_Explorer(f):
    fgmin = 1
    fgmax = 200
    for i in range(0, len(f)- 1):
	if f[i]<fgmin and f[i+1]>fgmin:
	    #print('MinHit')
	    f[i]= fgmin
	elif f[i]>fgmin and f[i]<fgmax and f[i+1]>fgmax:
	    continue;
	elif f[i]<fgmax and f[i+1]>fgmax:
	    #print('MaxHit')
	    f[i] = fgmax
	elif f[i]>fgmax:
	    #print('BreakHit')
	    f[i]=fgmax
	    break;
	if len(f) ==2:
	    f[0] = fgmin
	    f[1] = fgmax
    # print(' final f is ')
    # print f
    return [ (x >=  fgmin and x <= fgmax)*7.9e-14 for x in f]

def Einstein_Telescope(f):
    fgmin = 1
    fgmax = 200
    for i in range(0, len(f)- 1):
	if f[i]<fgmin and f[i+1]>fgmin:
	    #print('MinHit')
	    f[i]= fgmin
	elif f[i]>fgmin and f[i]<fgmax and f[i+1]>fgmax:
	    continue;
	elif f[i]<fgmax and f[i+1]>fgmax:
	    #print('MaxHit')
	    f[i] = fgmax
	elif f[i]>fgmax:
	    #print('BreakHit')
	    f[i]=fgmax
	    break;
	if len(f) ==2:
	    f[0] = fgmin
	    f[1] = fgmax
    # print(' final f is ')
    # print f
    return [ (x >=  fgmin and x <= fgmax)*2.5e-13 for x in f]

def CMB_LSS(f):
    fgmin = 1e-15
    fgmax = 1e10
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    #print(' final f is ')
    #print f
    #return [ (x >=  fgmin and x <= fgmax)*1.3e-5 for x in g]
    return [ (x >=  fgmin and x <= fgmax)*(1.0e-6/0.678**2) for x in g]
    
def BBN(f):
    fgmin = 1e-10
    fgmax = 1e10
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    #print(' final f is ')
    #print f
    return [ (x >=  fgmin and x <= fgmax)*1.5e-5 for x in g]
    
def Pulsar_Timing(f):
    fgmin = 1.0/20.0/(365*24*3600)
    fgmax = 1.0/(365*24*3600)
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    #print(' final f is ')
    #print f
    return [ (x >=  fgmin and x <= fgmax)*(2.0e-8/(0.678**2)) for x in g]

def Planck_Proj(f):
    fgmin = 1e-15
    fgmax = 1e10
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    #print(' final f is ')
    #print f
    return [ (x >=  fgmin and x <= fgmax)*2.7e-6 for x in g]
    
def COBE(f):
    fgmin = 3.0e-18
    fgmax = 1.0e-16
    omin = 1.4e-13*((1.0e-16/3.0e-18)**2)
    omax = 1.4e-13
    slop = math.log(omin/omax)/math.log(fgmin/fgmax)
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    #print(' final f is ')
    #print f
    return [((x>= fgmin and x <= fgmax))* (omin)*((x/fgmin)**slop) for x in g] 

def Inflaction(f, r):
    Mpl = 1.22e19 # planck mass in GeV/c^2
    Mplred = 2.435e18 # reduced planck mass in GeV/c^2
    As = math.exp(3.089) * 1e-10 # scalar spectrum amplitude, Planck XXII paper for Lambda CDM model
    OmegaM = 0.315  #using planck paper XVI results
    #OmegaM = 1; #Turner's value
    cc = 3e8  # speed of light in m/s
    hh = 0.673  #reduced Hubble parameter, using planck paper XVI result
    #hh = 1; #Turner's value
    H0 = 100.0 * hh  # Hubble parameter km/s/Mpc
    Mpc = 3.08567758 * 1e16 * 1e6  # in meters
    kH0 = 2* pi * H0 / Mpc * 1000.0 / cc  # wavenumber corresponding to H0. in m^{-1}

    gstar = 3.36; #effective number of relativistic DoF for CMB and neutrinos

    #calculate keq - there are several sources, 
    OmegaL = 0.685  #from Planck paper XVI
    zeq = 3402.0 #from Planck paper XVI
    #keq = H0 * (OmegaL + 2*OmegaM*(1+zeq)^3)^0.5 / (1 + zeq) * 1000 / Mpc / cc; %from Marco Peloso
    #keq = 0.03 * hh / Mpc; %in m^{-1}, from Annu. Rev. Astron. Astrophys. 1994. 32: 319-70
    #keq = 9.6e-17 * 2 * pi / cc; %Turner's value
    keq = 6.22e-2 * (OmegaM * hh**2 / np.sqrt(gstar / 3.36)) / Mpc #in m^{-1}, from Turner paper

    Vstar = 3 * pi**2 * As * r * Mplred**4 / 2.0 # from Planck XXII paper
    #const = 3 * pi^2 * As * r / 2;
    const = (1.94e16)**4 * r / 0.12 / Mpl**4
    nt = -r/8.0


    #ff = logspace(-18,6,100); #frequency in Hz 'logspace' is matlab function, find py equivalent if this line is needed
    kk = 2 * pi * f / cc # corresponding wavenumbers in m^{-1}

    #Turner paper
    omega = OmegaM**2 * const / (2*pi*kk/kH0)**(2-nt) * (1 + (4.0/3.0)*(kk/keq) + (5.0/2.0)*(kk/keq)**2)
    return omega

def Mag_Model(f, B=1E+15, eps=0.0002, P0=0.001, lamb=0.001, Iz=1E+45):

    # M = (1.4)*1.989e33     # gram
    R = 1.0e6                   # cm
    # Iz = 0.4*M*R**2         # g*cm**2
    # moment of inertia along z axis assumming a uniform perfect sphere.

    # Parameters
    fgmin = 0
    fgmax = 2.0/P0
    Kb = (192*np.pi**4*G)/(5*c**2*R**6)*(Iz)**3*eps**2/B**2
    Rb = Kb/(np.pi**2*Iz)
    Const = (8*np.pi*G)/(3*c**2*H0**3)*lamb*Kb/yr/Mpc**3 #Non integral part
    
    def integrand(z,v):
        numer = (1+z)**2*star_form_rate(sfr,z)
        rest_denom = 1 + Rb*v**2*(1+z)**2		# dimensionless
        denom = Ez(OmegaM,OmegaV,z)*rest_denom		# dimensionless
        return numer/denom

    if sfr != 'z':    
        omega = np.array(len(f)*[0], float)
        for x in range(len(f)):
            if  ((f[x] >=  fgmin) and (f[x]<= fgmax)):         
                if f[x] != 0 and fgmax/f[x] - 1 < 6 :
                    zmax = fgmax/f[x] - 1
                else:
			zmax = 6
        # z = 0:zmax/1000:zmax
                val, err = integrate.quad(integrand, 0, zmax, args=(f[x]))
                omega[x] = Const*f[x]**4*val
    else:
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
        rate = lambda z:star_form_rate(sfr, z, zvalues, sfrs)
        zrange = np.arange(0, 6+0.01, 0.01)
        rates = []
        for zval in zrange:
            rates.append(rate(zval)) 
        def integrand_part(z,v):
            numer = (1+z)**2
            rest_denom = 1 + Rb*v**2*(1+z)**2		# dimensionless
            denom = Ez(OmegaM,OmegaV,z)*rest_denom		# dimensionless
            return numer/denom
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
                    vals.append(integrand_part(zval, f[x])*rates[index])
                val = np.sum(vals) * 0.01
                omega[x] = Const*f[x]**4*val
    
    if sfr != 'b':
        return omega
    else:
        def integrand1(z,v):
            numer = (1+z)**2*star_form_rate('3',z)
            rest_denom = 1 + Rb*v**2*(1+z)**2		# dimensionless
            denom = Ez(OmegaM,OmegaV,z)*rest_denom		# dimensionless
            return numer/denom
        omega1 = np.array(len(f)*[0], float)
        for x in range(len(f)):
            if  ((f[x] >=  fgmin) and (f[x]<= fgmax)): 
                val, err = integrate.quad(integrand1, 0, zmax, args=(f[x]))
                omega1[x]=Const*f[x]**4*val
                omega[x] = omega[x] + omega1[x]
        return omega

    
def r_mod(f,P0,lamb):
    # R_MOD Summary of this function goes here
    # gravitational wave background from r mode
    m = 1.4
    M = m*Msolar
    R = 1.0E6
    P0 = P0*1.0e-3 # ms into s
    
    # keplerian velocity  (Ferrari 99) 
    # http://arxiv.org/pdf/astro-ph/9806357v2
    CC = 7.8*10**3 # (Friedmann, Ipser & Parker 1989)

    Kepler_v = CC*np.sqrt(m*10**6/R)
    Pk = 2*np.pi/Kepler_v*10**3  # of unit ms

    Kepler_f = 0.076*Kepler_v
    Pf = 2*np.pi/Kepler_f*10**3
    
    # gravitational frequency as a function of the period (kHz)
    fgmax = 4.0/3.0/Pk*10**3
    fgmin = 0
    Er = lambda x: (0.261 - 1.5*1.635*10**(-2))*M*(R**2)*(2*np.pi)**2*(1/x**2 - 1/Pf**2)
    Kr = lambda x: 2*Er(x)/((4.0/3.0/x)**2)
    dEr =lambda f,fsup: f*(f < fsup)
    Const = (8*np.pi*G)/(3*c**2*H0**3)*lamb/yr*Kr(P0)/Mpc**3
    
    def integrand(z,v):
        denom = Ez(OmegaM,OmegaV,z)
        integrand = dEr(v*(1+z),fgmax)*star_form_rate(sfr,z)/denom
        return integrand
        
    ## Finally get omega    
    omega = np.array(len(f)*[0], float)
    if sfr != 'z':
        for x in range(len(f)):
            if f[x] >=  fgmin and f[x] <= fgmax:
                if f[x] != 0 and fgmax/f[x] - 1 < 6 :
        		zmax = fgmax/f[x] - 1
                else:
    			zmax = 6
                val, err = integrate.quad(integrand, 0, zmax, args=(f[x]))
                omega[x] = Const*f[x]*val
    else:
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
        rate = lambda z:star_form_rate(sfr, z, zvalues, sfrs)
        zrange = np.arange(0, 6+0.01, 0.01)
        rates = []
        for zval in zrange:
            rates.append(rate(zval)) 
        def integrand_part(z,v):
            denom = Ez(OmegaM,OmegaV,z)
            integrand = dEr(v*(1+z),fgmax)/denom
            return integrand
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
                    vals.append(integrand_part(zval, f[x])*rates[index])
                val = np.sum(vals) * 0.01
                omega[x] = Const*f[x]*val
    if sfr != 'b':
        return omega
    elif sfr == 'b':
    
        def integrand(z,v):
            denom = Ez(OmegaM,OmegaV,z)
            integrand = dEr(v*(1+z),fgmax)*star_form_rate('3',z)/denom
            return integrand
        
    ## Finally get omega    
        omega3 = np.array(len(f)*[0], float)
   
        for x in range(len(f)):
            if f[x] >=  fgmin and f[x] <= fgmax:
                    val, err = integrate.quad(integrand, 0, zmax, args=(f[x]))
                    omega3[x] = Const*f[x]*val
                    omega[x] = omega[x] + omega3[x]
    return omega

    
def bar_Model(f,P0,lamb):
    # GW background from secular bar modes
    #   Detailed explanation goes here
    # M = 1.4*1.989*10**33
    # R = 10*10**5
    # betamin = 0.1375
    # betamax = 0.2738

    ## Fits from the paper of Lai and Shapiro
    # period = @(b) (1.127 - 4.621*b + 17.108*b**2)
    nusup = lambda P: (10**4*(-5.3 + 15.76*P - 15.5*P**2 + 5.1*P**3))
    
    Pmin = 0.8151
    Pmax = 1.1443
    fgmin = 0
    fgmax = nusup(Pmax)
    
    dEr = lambda x,y: x*(x < y)
    E0 = lambda P: 10**52*(-6.12 + 10.64*P - 3.84*P**2)
    Kr = lambda x: 2*E0(x)/(nusup(x)**2)*(x > Pmin)*(x < Pmax)
    
    Const = (8*np.pi*G)/(3*c**2*H0Mpc**3)*lamb/yr*Kr(P0)
    
    def integrand(z,v):
        denom = Ez(OmegaM,OmegaV,z)
        integrand = dEr(v*(1+z),fgmax)*star_form_rate(sfr,z)/denom
        return integrand
        
    ## Finally get omega    
    omega = np.array(len(f)*[0], float)
    
    if sfr != 'z':
        for x in range(len(f)):
            if f[x] >=  fgmin and f[x] <= fgmax:
                if f[x] != 0 and fgmax/f[x] - 1 < 6 :
        		zmax = fgmax/f[x] - 1
                else:
    			zmax = 6
	 
                val, err = integrate.quad(integrand, 0, zmax, args=(f[x]))
                omega[x] = Const*f[x]*val
    else:
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
        rate = lambda z:star_form_rate(sfr, z, zvalues, sfrs)
        zrange = np.arange(0, 6+0.01, 0.01)
        rates = []
        for zval in zrange:
            rates.append(rate(zval)) 
        def integrand_part(z,v):
            denom = Ez(OmegaM,OmegaV,z)
            integrand = dEr(v*(1+z),fgmax)/denom
            return integrand
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
                    vals.append(integrand_part(zval, f[x])*rates[index])
                val = np.sum(vals) * 0.01
                omega[x] = Const*f[x]*val
    if sfr != 'b':
           return omega
    else:

        def integrand(z,v):
            denom = Ez(OmegaM,OmegaV,z)
            integrand = dEr(v*(1+z),fgmax)*star_form_rate('3',z)/denom
            return integrand
        
    ## Finally get omega    
        omega1 = np.array(len(f)*[0], float)
    
        for x in range(len(f)):
            if f[x] >=  fgmin and f[x] <= fgmax:
                        val, err = integrate.quad(integrand, 0, zmax, args=(f[x]))
                        omega1[x] = Const*f[x]*val
            omega[x] = omega[x] + omega1[x]
        return omega
        
def PBB_Model(f,Mu,F1,Fs):
    H = 3.24078e-18 #Hubble
    M = 1.851256e43 #Planck Mass
    
    omega = np.array(len(f)*[0], float)
    for i in range(len(f)):
        a = 1/(1+np.sqrt(3))
        b = (a/48)*2**(2*Mu)*(2*Mu-1)**2*(sp.special.gamma(Mu))**2
        c = ((2*np.pi*Fs)**2/(H*M))**2*(F1/Fs)**(2*Mu+1)*(f[i]/Fs)**(5-2*Mu)
        HJ1 = hankel2(0,a*f[i]/Fs)*(jn(Mu-1,f[i]/Fs)-(Mu/f[i])*jn(Mu,f[i]/Fs))
        HJ2 = hankel2(1,a*f[i]/Fs)*jn(Mu,f[i]/Fs)
        HJ3 = ((1-a)/(2*a))*(Fs/f[i])*hankel2(0,a*f[i]/Fs)*jn(Mu,f[i]/Fs)
        omega[i] = b*c*(HJ1+HJ2-HJ3)*(HJ1+HJ2-HJ3).conj()*.5184
    
    return omega
    
def CS_Model(f, Prob, Size,Gmu):
    import subprocess
    EXCU_PATH = os.path.join(THIS_PATH,  "cs_lambda_stochastic")
    args = (EXCU_PATH, '-a', str(np.log10(f[1])), '-b', str(np.log10(f[-1])))
    args += ('-c', str(len(f)), '-d', str(np.log10(Gmu)))
    args += ('-e', str(np.log10(Gmu)), '-f', '1', '-g', str(np.log10(Size)))
    args += ('-i', str(np.log10(Size)), '-j', '1', '--log-pstart', str(np.log10(Prob)))
    args += ('--log-pend', str(np.log10(Prob)), '--np', '1', '--ln-zstart', '-10', '--ln-zend', '64', '--dlnz', '0.01', '--index', '1')
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    popen.wait()
    cc = np.loadtxt('stochastic_OmegaGW.dat',  comments='%')
    xp = cc[:,3]
    fp = cc[:,5]
    return np.interp(f,xp,fp, 0 , 0)

def CS_Model2(f, Prob, Size,Gmu):
    #doesnt work! wont update dat file from here (works from terminal)
    import subprocess
    EXCU_PATH = os.path.join(THIS_PATH,  "/home/fitzaxen/public_html/Kink_and_Cusp_Files/KINK_SL/myprogram")
    args = (EXCU_PATH, '-a', str(np.log10(f[1])), '-b', str(np.log10(f[-1])))
    args += ('-c', str(len(f)), '-d', str(np.log10(Gmu)))
    args += ('-e', str(np.log10(Gmu)), '-f', '1', '-g', str(np.log10(Size)))
    args += ('-i', str(np.log10(Size)), '-j', '1', '--log-pstart', str(np.log10(Prob)))
    args += ('--log-pend', str(np.log10(Prob)), '--np', '1', '--ln-zstart', '-10', '--ln-zend', '64', '--dlnz', '0.01', '--index', '1')
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    popen.wait()
    cc = np.loadtxt('/home/fitzaxen/public_html/Kink_and_Cusp_Files/KINK_SL/stochastic_OmegaGW2.dat',  comments='%')
    xp = cc[:,3]
    fp = cc[:,5]
    return np.interp(f,xp,fp, 0 , 0)

def CS_Model3(f, Prob,Gmu):
    import subprocess
    EXCU_PATH = os.path.join(THIS_PATH, "cs_lambda_stochasticLL")
    args = (EXCU_PATH, '-a', str(np.log10(f[1])), '-b', str(np.log10(f[-1])))
    args += ('-c', str(len(f)), '-d', str(np.log10(Gmu)))
    args += ('-e', str(np.log10(Gmu)), '-f', '1')
    args += ( '--log-pstart', str(np.log10(Prob)))
    args += ('--log-pend', str(np.log10(Prob)), '--np', '1', '--ln-zstart', '-10', '--ln-zend', '64', '--dlnz', '0.01')#, '--index', '1')
    #popen = subprocess.Popen(args, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
    #popen.wait()
    cc = np.loadtxt('/home/fitzaxen/public_html/Kink_and_Cusp_Files/CUSP_LL/stochastic_OmegaGW.dat',  comments='%')
    xp = cc[:,1]
    fp = cc[:,3]
    return np.interp(f,xp,fp, 0 , 0)

def CS_Model4(f, Prob,Gmu):
    import subprocess
    EXCU_PATH = os.path.join(THIS_PATH,  "cs_lambda_stochasticLLKINK2")
    args = (EXCU_PATH, '-a', str(np.log10(f[1])), '-b', str(np.log10(f[-1])))
    args += ('-c', str(len(f)), '-d', str(np.log10(Gmu)))
    args += ('-e', str(np.log10(Gmu)), '-f', '1')
    args += ( '--log-pstart', str(np.log10(Prob)))
    args += ('--log-pend', str(np.log10(Prob)), '--np', '1', '--ln-zstart', '-10', '--ln-zend', '64', '--dlnz', '0.01')
    #popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    #popen.wait()
    cc = np.loadtxt('/home/fitzaxen/public_html/Kink_and_Cusp_Files/KINK_LL/stochastic_OmegaGWKINK.dat',  comments='%')
    xp = cc[:,1]
    fp = cc[:,3]
    return np.interp(f,xp,fp, 0 , 0)

def axion(f, p ,cmb, Ncmb):
    import subprocess
    #print(THIS_PATH);  #output is /home/user1/gwplotter/Workspace/spectrum
#    os.system(' "module load mathematica/10.2" ')
#    os.system(' "pwd" ')
 #   os.system(' "./home/user1/gwplotter/public_html/module.sh" ')
    EXCU_PATH = os.path.join(THIS_PATH,  "axion_math")
    
##    print(EXCU_PATH)
  #  Mathematica = '/usr/local/bin/math'
    Mathematica = '/local/site/pkg/Wolfram/Mathematica/10.2/bin/math'  
    args = (Mathematica, '-script', EXCU_PATH, str(p), str(cmb), str(Ncmb), str(len(f)))
  ##  args = ('mathematica', '-script', EXCU_PATH, str(p), str(cmb), str(Ncmb), str(len(f)))
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    popen.wait()
    cc = np.loadtxt('axion_output.dat',  comments=' %')
    xp = cc[:,0]
    fp = cc[:,1]
    return np.interp(f,xp,fp, 0 , 0)
    
    # if p == 1:
    #     ll = np.loadtxt(os.path.join(THIS_PATH,'axtion_dat', 'll1.dat'), comments=' %')
    # else:
    #     ll = np.loadtxt(os.path.join(THIS_PATH,'axtion_dat', 'll2.dat'), comments=' %')
    # Xival = 4
    # Nval = 3
    # Xitab = [ ll[5*Nval*i + 1] for i in range(Xival)]  
    # Ntab = [ ll[5*i + 2,0] for i in range(Xival)]  
    # 
    # ll = os.path.join(THIS_PATH, 'axion_dat' "ll.dat")
    # f2 = os.path.join(THIS_PATH, 'axion_dat' "f2.dat")
    # f3 = os.path.join(THIS_PATH, 'axion_dat' "f3.dat")
    # fhl = os.path.join(THIS_PATH, 'axion_dat' "fhl.dat")
    # cc = np.loadtxt('axion_output.dat',  comments=' %')
    
def Dual_NS(f,M,lamb):
    #print(f)
    mm = M*Msolar # convert unit to gram 
    Kb = (G*np.pi)**(2.0/3.)*mm**(5.0/3.)/3.  # cm**2*g/s**(4/3)

    # Cutoffs of recieved frequence range 
    # per second
    # here 2**(6/5)*m is the total mass of the two NS
    

    fgmin = 10
    fgmax = c**3/(6.**1.5*np.pi*G*(2.**(6.0/5.)*mm)) 
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    #print(' final f is ')
    #print f
    #if f[1]==fgmax:
#	print('f[1]==fgmax')
    




    # dimensionless 
    # max of the redshift
    zmax = 6.  

    # normalizer of delay time distribution
    tmin = 0.02

    Const = (8.*np.pi*G)/(3.*c**2*H0**3)*lamb*Kb/yr/Mpc**3
    #omega = np.array(len(f)*[0], float)
    omega = np.zeros((len(g),1))
    
    if sfr != 'z':
        integrand = lambda z: BC_rate(z,tmin,sfr)/Ez(OmegaM,OmegaV,z)/(1+z)**(1.0/3.)
        for x in range(len(g)):
            if g[x] >=  fgmin and g[x] <= fgmax:
                if g[x] != 0 and fgmax/g[x] - 1 < 6 :
		#	print('inside hit ', x)
        		zsup = fgmax/g[x] - 1
                else:
    			zsup = 6.
                val, err = integrate.quad(integrand, 0, zsup)
                omega[x] = Const*g[x]**(2.0/3.)*val
    else:
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
        omega = np.array(len(g)*[0], float)
        for x in range(len(g)):
            if  ((g[x] >=  fgmin) and (g[x]<= fgmax)):         
                if g[x] != 0 and fgmax/g[x] - 1 < 6 :
                    zmax = fgmax/g[x] - 1
                else:
			zmax = 6
                vals = []
                zrange = np.arange(0, zmax+0.01, 0.01)
                for index,zval in enumerate(zrange):
                    vals.append(integrand_part(zval)*rates[index])
                val = np.sum(vals) * 0.01
                omega[x] = Const*g[x]**(2.0/3.)*val 
    if sfr != 'b':
        return omega
    else:
        omega1 = np.zeros((len(g),1))
    
        integrand = lambda z: BC_rate(z,tmin,'3')/Ez(OmegaM,OmegaV,z)/(1+z)**(1.0/3.)
   
        for x in range(len(g)):
                if g[x] >=  fgmin and g[x] <= fgmax:
                    val, err = integrate.quad(integrand, 0, zsup)
                    omega1[x] = Const*g[x]**(2.0/3.)*val
                    omega[x] = omega[x] +omega1[x]
        return omega

def Dual_BH(f, M, lamb, eta, zi):

    Kb = (G*np.pi)**(2.0/3.0)*(M*Msolar)**(5.0/3)/3

    Const = (8*np.pi*G)/(3*c**2*H0**3)*lamb*Kb/yr/Mpc**3

    zmin = 0.0
    zmax = 10.0 #already accounted for in BC_Rate function
    tmin = 0.05
    tmax = 13.5
    fgmin = 0.0
    #fgmax = v3 #maximum frequency

    #calulating total mass of black hole
    M = float(M)
    root = (0.25-eta)**(1.0/2.0)#eta cannot be more than 0.25
    fraction = (0.5 + root) / (0.5 - root)
    invfraction = 1.0 / fraction
    m2 = (M*(1+fraction)**0.2)/fraction**0.6
    m1 = (M*(1+invfraction)**0.2)/invfraction**0.6
    Mtot = m1 + m2

    v1 = lambda z: (404*(0.66389*eta**2 - 0.10321*eta + 0.10979)/0.125481)*(20.0/Mtot)/(1+z);
    v2 = lambda z: (807*(1.3278*eta**2 - 0.20642*eta + 0.21957)/0.250953)*(20.0/Mtot)/(1+z);
    v3 = lambda z:(1153*(1.7086*eta**2 - 0.26592*eta + 0.28236)/0.322668)*(20.0/Mtot)/(1+z);
    sigma = lambda z:(237 * (1.1383*eta**2 - 0.177*eta + 0.046834)/0.0737278)*(20.0/Mtot)/(1+z);

    omega = np.array(len(f)*[0], float)

    v_ = lambda f: (np.pi*Mtot*f*LAL_MTSUN_SI)**(1.0/3.0)

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

    if sfr != 'z':
        #integrand = lambda z,f: getv(z, f)*BC_rate(z,tmin, sfr)/Ez(OmegaM, OmegaV, z)/(1+z)**(1.0/3.0)
        integrand = lambda z,f: getv(z, f)/Ez(OmegaM, OmegaV, z)/(1+z)**(1.0/3.0)
        zrange = np.arange(zmin, zmax+0.1, 0.1)
        rates = []
        for zval in zrange:
            rates.append(BC_rate(zval, tmin, sfr))
        for x in range(len(f)):
            zrange2 = np.arange(zmin, zmax+0.01, 0.01)
            #val, err = integrate.quad(integrand, zmin, zmax, args=(f[x]))
            #omega[x]=Const*f[x]*val
            vals = []
            for index,zval in enumerate(zrange2):
                val1 = np.interp(zval, zrange, rates)
                vals.append(integrand(zval, f[x])*val1)
            val = np.sum(vals) * 0.01
            omega[x] = Const*f[x]* val
    else:
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
            omega[x] = Const*f[x]* val

    if sfr != 'b':
        return omega
    else:
        omega1 = np.array(len(f)*[0], float)
        #integrand1 =lambda z,f: getv(z, f)*BC_rate(z,tmin, '3')/Ez(OmegaM, OmegaV, z)/(1+z)**(1.0/3.0)
        integrand1 = lambda z, f: getv(z,f)/Ez(OmegaM, OmegaV, z)/(1+z)**(1.0/3.0)
        zrange = np.arange(zmin, zmax+0.1, 0.1)
        rates = []
        for zval in zrange:
            rates.append(BC_rate(zval, tmin, '3'))
        for x in range(len(f)):
            #val, err = integrate.quad(integrand1, zmin, zmax, args=(f[x]))
            zrange2 = np.arange(zmin, zmax+0.01, 0.01)
            vals = []
            for index,zval in enumerate(zrange2):
                val1 = np.interp(zval, zrange, rates)
                vals.append(integrand(zval, f[x])*val1)
            val = np.sum(vals) * 0.01
            omega1[x]=Const*f[x]*val
            omega[x] = omega[x] + omega1[x]
        return omega

def NSBH(f, M,lamb):
    # Stochastic Energy Spectrum of NS-BH Model
    # f - recieved frequency
    # m - average chirp mass in solar mass 
    # lamb - mass fraction rate
    # sfr - star formation history 
    
    mm = M*Msolar # convert unit to gram 
    Kb = (G*np.pi)**(2.0/3.0)*mm**(5.0/3.0)/3.0  # cm**2*g/s**(4/3)

    # Cutoffs of recieved frequence range 
    # per second
    # here 2**(6/5)*m is the total mass of the two NS
    fgmin = 10
    fgmax = c**3/(6**1.5*np.pi*G*(2**(6.0/5)*mm)) 

    # dimensionless 
    # max of the redshift
    zmax = 6   

    # normalizer of delay time distribution
    tmin = 0.1

    Const = (8*np.pi*G)/(3*c**2*H0**3)*lamb*Kb/yr/Mpc**3
    omega = np.array(len(f)*[0], float)
    if sfr != 'z':
        integrand = lambda z: BC_rate(z,tmin,sfr)/Ez(OmegaM,OmegaV,z)/(1+z)**(1.0/3)
        for x in range(len(f)):
            if f[x] >=  fgmin and f[x] <= fgmax:
                if f[x] != 0 and fgmax/f[x] - 1 < 6 :
        		zsup = fgmax/f[x] - 1
                else:
    			zsup = 6
                val, err = integrate.quad(integrand, 0, zsup)
                omega[x] = Const*f[x]**(2.0/3.0)*val
    else:
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
                omega[x] = Const*f[x]**(2.0/3.)*val   
    if sfr != 'b':
        return omega
    else:
        omega1 = np.array(len(f)*[0], float)
        integrand = lambda z: BC_rate(z,tmin,'3')/Ez(OmegaM,OmegaV,z)/(1+z)**(1.0/3)
        for x in range(len(f)):
            if f[x] >=  fgmin and f[x] <= fgmax:
                val, err = integrate.quad(integrand, 0, zsup)
                omega1[x] = Const*f[x]**(2.0/3.0)*val
                omega[x] = omega[x] + omega1[x]
        return omega
    
def NS_turb(f, d0m=100, N=1.0, v=1.0, M_star=1.4, R_star=10.0):###find c to c_
    # Stochastic Energy Spectrum of the NS Turbulence Model
    #  d0m   - Crust core differential shear
    #   N    - Neutron Star Fraction
    #   v    - kinematic viscosity
    # M_star - NS mass
    # R_star - NS radius
    #d0m=200
    #N=2.0
    #v=2.0
    #M_star=2.8
    #R_star=20.0
    ztest = arange(0,105,5)
    M_sun = 1.9891e30			#mass of the sun in kg
    R_star = R_star * 1.0e3 #convert to meters
    # Calculate some constants...
    ks = 2.0 * np.pi / R_star		# Stirring wavenumber - Defined after eqn (7)
    eps = (R_star ** 2) * (d0m ** 3)	# energy dissipation rate
    eta = (1.0 / np.sqrt(2 * np.pi)) * (eps ** (1.0 / 3.0)) * (ks ** (2.0 / 3.0))# defined in eqn (7)
    kd = (8 ** (1.0 / 4.0)) * (eps ** (1.0 / 4.0)) * (27 ** (-1.0 / 4.0)) * (v ** (-3.0 / 4.0))# viscous dissipation wavenumber

    ##Calculate N(z)... (comments found in astrofunct.py file)
    H = 73000.0 # Hubble Constant in m/(s*Mpc) - note the weird units. converted to float
    rhoc = (3 * (H ** 2)) / (8 * np.pi * G)
    c_ = 299792458.0 #m/s
    Om_m = 0.26
    Om_L = 0.74
  
    fEE = lambda z_: np.sqrt(Om_m * ((1 + z_)**3) + Om_L)

    ##Distances.
    dH = c_ / H #Hubble distance (Mpc).
    ddCdz = lambda z_, dC_: dH / fEE(z_)
    zspan = arange(0,100.2, .1)
    solver = ode(ddCdz).set_integrator('vode', atol = 1, rtol = 1e-6, max_step = .1) 
    solver.set_initial_value(0,0)
    dC = np.zeros((1001,2))
    dt = .1
    iterator = 1
    while solver.successful() and solver.t<99.9: 
	solver.integrate(solver.t+dt)
	dC[iterator, 0]=solver.t
	dC[iterator, 1] = solver.y[0]
	iterator = iterator+1
    end
    fdC = lambda z_: np.interp(z_, dC[:,0], dC[:,1])    
    dM = lambda z_: fdC(z_)
    dL = lambda z_: (1 + z_) * dM(z_)
    dA = lambda z_: dM(z_) / (1 + z_)
    dVdz = lambda z_: 4.0 * np.pi * ((fdC(z_)) ** 2) * dH / fEE(z_)
    Vc = lambda z_: (4.0 * np.pi / 3.0) * (dM(z_) ** 3)
    
    ##IMF
    IMFIntegrand = lambda M: M ** -2.35
    IMFscale = 1.0 / double(quad(IMFIntegrand, 0.1, 100)[0]) 
    fIMF = lambda M_: IMFscale * (M_ ** -2.35) 
    Mmin = 8.0 #unit of M_sun
    Mmax = 40.0 #unit of M_sun
    IMFint = double(quad(fIMF, 8.0, 40.0)[0])
    
    ##starformationrate
    aHB6 = 0.0170
    bHB6 = 0.13
    cHB6 = 3.3
    dHB6 = 5.3
    h100 = H / (1000.0 * 100.0) #converted to float
    fSFR = lambda z_: (aHB6 + bHB6 * z_) * h100 / (1 + (z_ / cHB6) ** dHB6)

    z1 = arange(1,100.000000000001,.1)
    SFR = map(fSFR, z1)
   
    ##diffsourceformrate
    fdNdz = lambda z_: IMFint * fSFR(z_) * dVdz(z_) / (1 + z_)
    z = arange(0,100.1, .1)
    dNdz = map(fdNdz, z)

    ##numneutronstar
    Hyr = (73.0 / (3.0857e19)) * 60.0 * 60.0 * 24.0 * 365.25# Hubble constant in units 1/yr
    fdNumNSdz = lambda z_, NumNS_: (1 / (Hyr * M_star)) * fdNdz(z_) / ((1 + z_) * fEE(z_))
    ztest2= arange(0,20,1)
    solver2 = ode(fdNumNSdz).set_integrator('vode', atol = 100, rtol = 1e-6)
    solver2.set_initial_value(0,0)
    NumberNS2 = np.zeros((1001,2))
    iterator2 = 0
    while solver2.successful() and solver2.t<100: 
	  solver2.integrate(solver2.t+dt)
	  NumberNS2[iterator2, 0] = solver2.t
	  NumberNS2[iterator2, 1] = solver2.y[0]
	  iterator2 = iterator2 + 1
    end
    zadjusted = arange(0,100.2, .1)
    fNumberNS = lambda z_: np.interp(100, NumberNS2[7:,0], NumberNS2[7:,1]) - np.interp(z_, NumberNS2[7:,0], NumberNS2[7:,1])
    fNumNSdVdz = lambda z_: fNumberNS(z_) / dVdz(z_)
    z__ = np.linspace(0,15,1000)
    NumNSdVdz = map(fNumNSdVdz, z__)
    fdNumDenNSdz = lambda z_, NumDenNS_: (IMFint / (Hyr * M_star)) * (fSFR(z_) / (1 + z_))
    solver3 = ode(fdNumDenNSdz).set_integrator('vode', atol = 1, rtol = 1e-6)
    solver3.set_initial_value(0,0)
    EENumDenNS = np.zeros((1001, 2))
    iterator3 = 0
    while solver3.successful() and solver3.t<100:#never goes through loop
	solver3.integrate(solver3.t+dt)	
	EENumDenNS[iterator3,0]=solver3.t
	EENumDenNS[iterator3,1] = solver3.y[0]
	iterator3 = iterator3 + 1
    end
    fNumDenNS = lambda z_: (1. / (fEE(z_) * (1.0 + z_))) * (np.interp(20, EENumDenNS[:,0], EENumDenNS[:,1]) - np.interp(z_, EENumDenNS[:,0], EENumDenNS[:,1]))
    z__ = np.linspace(0,100,1000) 
    NumDenNS = map(fNumDenNS, z__)
    ztest = arange(0,1050,50)
    #Calculate Omega_gw... (comments in astrofunct.py define calculatedEdSdv)
    Mpc2m = 3.08e22
    fdVdzmetres = lambda z_: dVdz(z_) * ((Mpc2m) ** 3)
    fVcmetres = lambda z_: Vc(z_) * ((Mpc2m) ** 3)
    fVLmetres = lambda z_: (4.0 * np.pi / 3.0) * (dL(z_) ** 3) * ((Mpc2m) ** 3)
    fdCmetres = lambda z_: fdC(z_) * (Mpc2m)
    fdLmetres = lambda z_: dL(z_) * (Mpc2m) 
    Hsec = H / Mpc2m
    pc = 3.0 * (Hsec ** 2) / (8.0 * np.pi * G)
    fh_rms_redshift = lambda z_: 5.0e-28 * (M_star / 1.4) * ((R_star / 1.0e4) ** 3) * ((d0m / 10.0) ** 3) * (1. / (dL(z_) * 1000))
    fvv = lambda z_: Om_m * ((1 + z_) ** 3)
    fLL = lambda z_: log((np.sqrt(fvv(z_) + Om_L) + np.sqrt(Om_L)) / (np.sqrt(fvv(z_) + Om_L) - np.sqrt(Om_L)))
    fTLookBack = lambda z_: (fLL(0) - fLL(z_)) / (3.0 * Hsec * np.sqrt(Om_L))
    fdTdz = lambda z_: 1.0 / (Hsec * (1.0 + z_) * np.sqrt(fvv(z_) + Om_L))
    fdEdSdv1 = lambda ve_, z_: (1.0 / ((4.0 * sqrt(2.0)) * (np.pi ** (7.0 / 2.0)) * (d0m ** 3))) * exp(-(ve_ ** 2) / ((2 ** (1.0 / 3.0)) * (np.pi ** (4.0 / 3.0)) * (d0m ** 2))) * (fh_rms_redshift(z_) ** 2) * (ve_ ** 2) * (-3.0 * (2 ** (1.0 / 3.0)) * (np.pi ** (4.0 / 3.0)) * (d0m ** 2) + (7.0 * (ve_ ** 2)))
    B1 = 0.5 * np.sqrt(np.pi) * eta * (kd ** (2.0 / 3.0)) * (ks ** (-2.0 / 3.0))
    fdEdSdv2 = lambda ve_, z_: exp(-(ve_ ** 2) / (4.0 * (B1 ** 2))) * (fh_rms_redshift(z_) ** 2) * np.sqrt(np.pi / 2.0) * (d0m ** 3) * (48.0 * (B1 ** 6) + 12.0 * (B1 ** 4) * (ve_ ** 2) + (ve_ ** 6)) / (16.0 * (B1 ** 6) * (ve_ ** 2))
    B2 = 0.5 * np.sqrt(np.pi) * eta
    fdEdSdv3 = lambda ve_, z_: -exp(-(ve_ ** 2) / (4.0 * (B2 ** 2))) * (fh_rms_redshift(z_) ** 2) * np.sqrt(np.pi / 2.0) * (d0m ** 3) * (48.0 * (B2 ** 6) + 12.0 * (B2 ** 4) * (ve_ ** 2) + (ve_ ** 6)) / (16.0 * (B2 ** 6) * (ve_ ** 2))
    fdEdSdv = lambda ve_, z_: (fdTdz(z_) * (c_ ** 3) / G) * (fdEdSdv1(ve_, z_) + fdEdSdv2(ve_, z_) + fdEdSdv3(ve_, z_))
  
    # Use exponential steps to product the constant steps on log-scale
    # The min and max values may need to be input from the website.
    step = 0.1	# step size as appear on log-scale
    min = 1.0e-5	# minimum value of v_ (must be positive)
    max = 1.0e5		# maximum value of v_
    v_peak = d0m / (10.0 * (26.0e-3))
    if min < 1.0e-2 * v_peak:    
        min = 1.0e-2 * v_peak    
    end
    
    i = arange(log10(min), log10(max)+ 0.0001, step)
    v_ = 10. ** i# generate exponential step for v

    

    Omegagw2 = np.zeros((len(v_), 1))# reserve the memory for computation result
 

    MPE = 1e-10# Maximum Percentage Error
    tolerance = 1e-25# initial integration error tolerance
    
    # Main loop for calculating Omega_gw as function of frequency:    
    for i_ in range(0,len(v_)):
        ffdEdSdv = lambda z_: fdEdSdv(v_[i_] * (1 + z_), z_)
	integrand = lambda z_: (5.0 / (4. * np.pi)) * (4. * np.pi * (fdLmetres(z_) ** 2) * ffdEdSdv(z_) * v_[i_] * (1.0 + z_)) * (fNumDenNS(z_) / ((1.0 + z_) * ((Mpc2m) ** 3.0))) #error shows up here, test other functions
	a = (1 / (pc * (c_ ** 2))) * double(quad(integrand,0,20, epsabs = tolerance))
	Omegagw2[i_] = a[0]
        tolerance = MPE * Omegagw2[i_]
    end
    omega = lambda v__: N * np.interp(v__, v_, Omegagw2[:,0])
    CalledOmega = omega(f)
    counter = 0
    for i in range(0,len(CalledOmega)-2):			#for loop removes unwanted tails.
	if CalledOmega[i] == CalledOmega[i+1] and counter == i:
	    CalledOmega[i] = 0
	    counter = counter +1
	elif CalledOmega[i] ==CalledOmega[i+1] and counter !=i:
	    for x in range(i+1, len(CalledOmega)-1):
		  CalledOmega[x] = 0
	    end
	    break
    end
    #CalledOmega[len(CalledOmega)-1] = 0 #this cuts off data point we want in some cases (depending on range). fix by new if statement that tests last two, deletes last if they are equal.
    #print(CalledOmega)

    return CalledOmega

def eqstate_model(f, r, w, n):
    #This is a function to make plot as in astro-ph/0708.2279v2
    #Written by shivaraj 05/2008

    # C2 and C3 are taken to be order 1 and other values are from that paper
    C2=1.0
    C3=1.0
    DelR=sqrt(2.24*10**(-9))
    h=0.72;
    #r=0.1;
    # g1zc is taken as g1zbbn in the limiting case; Here 1 corresponds to *
    g1zc=10.75;
    g1szeq=3.9091;
    g1zeq=3.3626;
    # g1szc is taken as g1szbbn in the limiting case
    g1szc=10.75
    H0=3.24*h*10**(-18)
    # zbbn1=zbbn+1
    zbbn1=5.9*10**9

    # f=fend; The following expression is an approximate one
    #f=4.5*10^8*r^(1/4);

    # Original expression is "gamma=(Omegamat/(1+zeq))*g1zc*g1szeq^(4/3)/
    #                                               (g1zeq*g1szc^(4/3))"
    # with 1+zeq is defined as 2.3*10**4*Omegamat*h**2
    gamma=(1.0/(2.3*10**4*h**2))*g1zc*g1szeq**(4.0/3.0)/(g1zeq*g1szc**(4.0/3.0))
    A1=C2*C3*DelR**2*gamma/24.0
    # Original expression is "A2=(2*pi*f/H0)*1/((1+zc)*gamma^0.5)" with
    # 1+zc=1+zbbn=zbbn1
    A2=(2.0*np.pi*f/H0)*1.0/(zbbn1*gamma**0.5)
    # Original expression is "A3=(2*pi*f/H0)*(H0*a/k)" where k/(H0*a)=150.0/h 
    # with k=0.05 Mpc^(-1);
    A3=(2.0*np.pi*f/H0)*(h/150.0)

    alpha=2.0 * (3.0*w - 1.0) / (3.0*w + 1.0)

    omega = A1 * A2**alpha * A3**n * r
    return omega

def eLisa(f):
    #sensitivity to Omega_GW for eLISA
    fgmin = 1.0e-5
    fgmax = 1.0
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax 


    #f = 1e-5: 1e-5: 1;

    s1 = 2.13e-29 * (1.0 + 1.0e-4/f)
    s2 = 1.37e-32 * (1.0 + 1.0e-4/f) / f**4
    s3 = 5.25e-23
    s4 = 6.28e-23

    ss = sqrt( (20.0/3.0) * (4.0*s1/(2.0*np.pi*f)**4 + s3 + s4) / (1.0e9)**2.0 * (1.0 + (f / (0.41 * 3.0e8 / 2.0/ 1.0e9))**2 ) )
    #ss = sqrt( (20/3) * (4*s1./(2*pi*ff).^4 + s3 + s4) / (1e9)^2 .* (1 + (ff ./ (0.41 * 3e8 / 2/ 1e9)).^2 ) ) / 7.57e3;

    sens = (ss / 6.3e-22 * (g / 100)**(3.0/2.0))**2
    omega = sens
    for i in range(len(omega)):
	if g[i]<fgmin or g[i]>fgmax:
	    omega[i] = 0
    return omega
    
def White_Dwarf_Binary(f):
    fgmin = 1.0e-4
    fgmax = 1.0e-2
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    freq = [1.0e-4, 3.0e-4, 1.0e-3, 2.0e-3, 3.0e-3, 1.0e-2]
    strain = [3.0e-18, 2.0e-18, 5.0e-19, 4.0e-19, 1.5e-20, 1.0e-21] #   Hills-Bender
    freq_array = np.zeros((6,1))
    for i in range(len(freq_array)):
	freq_array[i] = freq[i]
    #print freq_array[0,0]
    strain_array = np.zeros((6,1))
    for i in range(len(strain_array)):
	strain_array[i] = strain[i]
    #print strain_array[0,0]
    omeg = (strain_array / 6.3e-22 * (freq_array / 100.0)**(3.0/2.0))**2
    #print omeg[0,0]
    omega = np.interp(g, freq_array[:,0], omeg[:,0])
    for i in range(len(omega)):
	if g[i]<fgmin or g[i]>fgmax:
	    omega[i] = 0
    return omega

def Supermassive_BBH(f):
    fgmin = 3.0e-9
    fgmax = 2.0e-7
    g = f
    for i in range(0, len(g)- 1):
	if g[i]<fgmin and g[i+1]>fgmin:
	    #print('MinHit')
	    g[i]= fgmin
	elif g[i]>fgmin and g[i]<fgmax and g[i+1]>fgmax:
	    continue;
	elif g[i]<fgmax and g[i+1]>fgmax:
	    #print('MaxHit')
	    g[i] = fgmax
	elif g[i]>fgmax:
	    #print('BreakHit')
	    g[i]=fgmax
	    break;
	if len(g) ==2:
	    g[0] = fgmin
	    g[1] = fgmax
    freq = [3.0e-9, 2.0e-7]
    hc = [3.0e-15, 2.0e-16]
    H0 = 70000.0/3.1e22
    freq_array = np.zeros((2,1))
    for i in range(len(freq_array)):
	freq_array[i] = freq[i]
    #print freq_array[0,0]
    hc_array = np.zeros((2,1))
    for i in range(len(hc_array)):
	hc_array[i] = hc[i]
    #print strain_array[0,0]
    omeg = (2.0*np.pi**2 / 3.0 / H0**2) * freq_array**2 * hc_array**2
    #print omeg[0,0]
    omega = np.interp(g, freq_array[:,0], omeg[:,0])
    for i in range(len(omega)):
	if g[i]<fgmin or g[i]>fgmax:
	    omega[i] = 0
    return omega


def test_func(f, x, d0m):
    omega = np.zeros((len(f),1))
    for i in range(0, len(f)):
	omega[i]=f[i]*f[i]-x*8*d0m
    print omega
    #print f[5]
    return omega

def CCBH(f, epsilon, alpha, a, m_low, m_high):
  m_min = 0.1*Msolar #dummy variable for m_min

  const_ = 8.0*pi*G*epsilon*alpha / (3.0*H0**3) ##constant term in front of integral

  #initial mass function definition--Salpeter
  N = 0.35*m_min**0.35; #normalization constant
  IMF = lambda m: N*m**(-2.35)#steleter mass function
  v_ = lambda m: (c**3)*(1-0.63*(1-a)**0.3)/(2.0*pi*G*alpha*m)#dominant emission frequency
  Omega = np.zeros ((len(f), 1));
  m_low = m_low*Msolar*1.0 #take m_low in terms of solar mass.
  #convert star formation rate to g/s/cm^3
  conv_fact = Msolar/3.15569e7/((3.08567758e24)**3)
  #display(v_(m_low));
  if sfr != 'z':
      for x in range(len(f)): 
	z_ = lambda m: v_(m)/f[x] - 1	  
	if z_(m_low) <= 0:   #make sure z prime is positive
	    continue;
	if m_low < 4.0*Msolar/alpha:       #make sure there is enough mass to create black hole
	    m_low = 4.0*Msolar/alpha
	#low bound for mass integral_result

	m_h = (c**3)*(1-0.63*(1-a)**0.3)/(2.0*pi*G*alpha*f[x])
	if m_h > m_high*Msolar:
	    m_h = m_high*Msolar
	if m_low>m_h:
	    #display('hit');
	    continue;
	integrand = lambda m: conv_fact*star_form_rate(sfr, z_(m))*IMF(m)*m/(1+z_(m))/Ez(OmegaM, OmegaV, z_(m))
	integral_result = quad(integrand,m_low, m_h)[0]
	Omega[x] = const_*integral_result
  else:
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
      for x in range(len(f)):
	z_ = lambda m: v_(m)/f[x] - 1	  
	if z_(m_low) <= 0:   #make sure z prime is positive
	    continue;
	if m_low < 4.0*Msolar/alpha:       #make sure there is enough mass to create black hole
	    m_low = 4.0*Msolar/alpha
	#low bound for mass integral_result

	m_h = (c**3)*(1-0.63*(1-a)**0.3)/(2.0*pi*G*alpha*f[x])
	if m_h > m_high*Msolar:
	    m_h = m_high*Msolar
	if m_low>m_h:
	    #display('hit');
	    continue;
        #integrand = lambda m: conv_fact*star_form_rate(sfr, z_(m), zvalues, sfrs)*IMF(m)*m/(1+z_(m))/Ez(OmegaM, OmegaV, z_(m))
	def integrand(m):
            if z_(m) <= 10:
                inte = conv_fact*star_form_rate(sfr, z_(m), zvalues, sfrs)*IMF(m)*m/(1+z_(m))/Ez(OmegaM, OmegaV, z_(m));
            else:
                inte = conv_fact*star_form_rate('k', z_(m))*IMF(m)*m/(1+z_(m))/Ez(OmegaM, OmegaV, z_(m));
            return inte
        integral_result = quad(integrand,m_low, m_h)[0]
	Omega[x] = const_*integral_result
  if sfr != 'b':
            return Omega
  else:
      m_min = 36.0*Msolar #dummy variable for m_min

      const_ = 8.0*pi*G*epsilon*alpha / (3.0*H0**3) ##constant term in front of integral

  #initial mass function definition--Salpeter
      N = 0.35*m_min**0.35; #normalization constant
      IMF = lambda m: N*m**(-2.35)
      v_ = lambda m: (c**3)*(1-0.63*(1-a)**0.3)/(2.0*pi*G*alpha*m)
      Omega1 = np.zeros ((len(f), 1));
      m_low = m_low*Msolar*1.0 #take m_low in terms of solar mass.
  #convert star formation rate to g/s/cm^3
      conv_fact = Msolar/3.15569e7/((3.08567758e24)**3)
  #display(v_(m_low));
      for x in range(len(f)): 
          z_ = lambda m: v_(m)/f[x] - 1	  
          if z_(m_low) <= 0:   #make sure z prime is positive
              continue;
          if m_low < 4.0*Msolar/alpha:       #make sure there is enough mass to create black hole
              m_low = 4.0*Msolar/alpha
	#low bound for mass integral_result

	#z_ must be greater than 0
          m_h = (c**3)*(1-0.63*(1-a)**0.3)/(2.0*pi*G*alpha*f[x])
          if m_h > m_high*Msolar:
              m_h = m_high*Msolar
          if m_low>m_h:
	    #display('hit');
              continue;
	#while z_(m_high+1.0e33)>0
	 #   m_high = m_high + 1.0e35;
	 #   display(z_(m_high));
	#end
	#3define z prime, the z for which the fermi-dirac is nonzero
	
	#print(z_(m_low));
	#print(star_form_rate(sfr,z_(m_low)));
	#print(IMF(m_low));
	#display(m_low./(1+z_(m_low)));
	#display(Ez(OmegaM, OmegaV, z_(m_low)));
          integrand = lambda m: conv_fact*star_form_rate('3', z_(m))*IMF(m)*m/(1+z_(m))/Ez(OmegaM, OmegaV, z_(m))
	# print integrand(m_low)
          integral_result = quad(integrand,m_low, m_h)[0]
	#print(integral_result)
          Omega1[x] = const_*integral_result
          Omega[x] = Omega1[x] + Omega[x]
      return Omega

def CCBH_full(f, aa, bb, lamb, Eneut, q):
    #simulations suggest that should fall aa:5-150, b:40-400- if user enters negative value, put into reasonable range
    if aa <= 0:
        aa = 5
    if bb <= 0:
        bb = 40
    zi = G*lamb*Eneut**2*q**2*10**4/c**5 #includes converting zi to correct units
    #zi = 1
    const = 8.0*pi*G*zi/(3*H0**3*c**2)
    Omega = np.zeros((len(f), 1))
    z_low = 0.0
    z_high = 20.0
    conv_fact = Msolar/3.15569e7/((3.08567758e24)**3);
    for x in range(len(f)):
        f_p = lambda z: f[x]*(1+z)
    if sfr != 'z':
        for x in range(len(f)):
           integrand = lambda z:star_form_rate(sfr, z)*conv_fact*(1+f_p(z)/aa)**6*exp(-2*f_p(z)/bb)/(1+z)/Ez(OmegaM, OmegaV, z)
        #print(integrand(z_low))
           integral_result = quad(integrand, z_low, z_high)[0]
           Omega[x] = const*f[x]*integral_result
    else:
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
        rate = lambda z:star_form_rate(sfr, z, zvalues, sfrs)
        integrand_part = lambda z:conv_fact*(1+f_p(z)/aa)**6*exp(-2*f_p(z)/bb)/(1+z)/Ez(OmegaM, OmegaV, z)
        zrange = np.arange(z_low, z_high+0.01, 0.01)
        rates = []
        for zval in zrange:
            rates.append(rate(zval)) 
        for x in range(len(f)):
            vals = []
            for index,zval in enumerate(zrange):
                vals.append(integrand_part(zval)*rates[index])
            val = np.sum(vals) * 0.01
            Omega[x] = const*f[x]* val
   
    if sfr != 'b':
        return Omega
    else:
        Omega1 = np.zeros((len(f), 1))
        for x in range(len(f)):
            f_p = lambda z: f[x]*(1+z)
            integrand1 = lambda z:star_form_rate('3', z)*conv_fact*(1+f_p(z)/aa)**6*exp(-2*f_p(z)/bb)/(1+z)/Ez(OmegaM, OmegaV, z)
            integral_result1 = quad(integrand1, z_low, z_high)[0]
            Omega1[x] = const* f[x]* integral_result1
            Omega[x] = Omega1[x]+ Omega[x]
        return Omega

def CCBH_metal(f, a, epsilon):
    LENGTH = len(f);
    Omega = np.zeros(LENGTH);
    fil = open("/home/user1/gwplotter/Workspace/spectrum/input/parameters.dat","w")
    #fil = open("/home/fitzaxen/Workspace/spectrum/input/parameters.dat","w")
    fil.write(str(LENGTH))
    fil.write('\n')
    for i in range(0, LENGTH):
        fil.write(str(f[i]))
        fil.write('\n')
    fil.write(str(a))
    fil.write('\n')
    fil.write(str(epsilon))
    fil.write('\n')
    fil.write(sfr)
    fil.write('\n')
    fil.close();

    if fil.closed:

        os.system(' "./unique_executable" ') #in gwplotter
        #os.system(' "./unique_executable4" ') # in fitzaxen

        #Note: to create file from personal directory, duplicate/modify file CCBH_metal_low_mem.cpp in public_html. Then create executable file with the following command: g++ -m32 CCBH_metal_low_mem2.cpp -o unique_executable -L/usr/lib -lstdc++
        #Remember if you are copying this file over to gwplotter directory to have appropriate gwplotter directories uncommented or the spectrum wont work

    fil2=open("/home/user1/gwplotter/Workspace/spectrum/results/Omega_from_C.dat","r")
    #fil2=open("/home/fitzaxen/Workspace/spectrum/results/Omega_from_C.dat","r")
    if fil2.closed:
        print "Error: unable to open data file."
    if sfr == 'k' or sfr == 'b':
        for i in range (0, LENGTH):
            Omega[i] = float(fil2.readline())
    fil2.close();
    return Omega;

def phase_transition(f, T, g, v, betaH):

    K = 1.0
    m_Pl = 1.22e19 #Gev

    if betaH < 10 or betaH > 1000:
        betaH = np.log(m_Pl/T)

    fbeta = 0.62/(1.8 - 0.1*v + v**2)

    f0 = (1.65e-7)*fbeta*betaH*T*(g/100.0)**(1.0/6.0)
    omega0h2 = (1.67e-5)*K**2 *1.0*((0.11*v**3)/(0.42+v**2))*(1.0/betaH)**2*(100.0/g)**(1.0/3.0)

    p = 2.8
    q = 1.0

    Omega = np.zeros(len(f))
    for x in range(len(f)):
        Omega[x] = omega0h2*(p+q)*(f[x]/f0)**p/(q+ p*(f[x]/f0)**(p+q))/h0**2
 
    return Omega

def cosmic_strings_small_cusps(f,Gu,epsilon, p): 
    
    gamma = 50.0
    t0 = 10**17.5
    alpha = epsilon * gamma * Gu 
    c = 1
    constant = 2*c*Gu*(np.pi)**2.0/3.0/p*H0**(1.0/3.0)/gamma/alpha**(1.0/3.0)
    z_eq = 5440.0 

    Ez = lambda OmegaM, OmegaV,z:(OmegaM*(1.0+z)**3.0+OmegaV)**0.5

    c_z = lambda z:1.0+ (9.0*z/(z+z_eq))

    def cut(z,f):
        if 1.0-(f*(1.0+z)*alpha*phi_t(z)/H0)**(-1.0/3.0) < 0.0:
           return 0.0
        else:
           return 1.0

    def cut2(z,f):
        if z < 1:
            cons = 54.0*np.pi/4.0
        elif z < z_eq:
            cons = 54.0 * np.pi
        else:
            cons = 72.0*np.pi
        if cons/t0/alpha**(8.0/3.0)/(f*t0)**(2.0/3.0)*phi_n(z)/f < 1.0:
            return 0.0
        else:
            return 1.0

    cc = np.loadtxt('/home/fitzaxen/Workspace/spectrum/cosmo_approx.txt')
    zvals = cc[:,0]
    cc2 = np.loadtxt('/home/fitzaxen/Workspace/spectrum/cosmo_approx2.txt')
    zvals2 = cc2[:,0]
    

    def phi_t(z):
        if z < z_eq:
            #return cc[:,0][z]
            return np.interp(z, zvals, cc[:,1])
        else:
            #return integrate.quad(lambda z_:1.0/(1.0+z_)/Ez(OmegaM, OmegaV, z_),z,np.inf)[0]
            return (z_eq)**(1.0/2.0)* z**-2.0
           #return np.interp(z, zvals2, cc2[:,1])

    def phi_r(z):
        if z < z_eq:
         #return integrate.quad(lambda z_:1.0/Ez(OmegaM, OmegaV, z_), 0.0,z)[0]
            #return cc[:,1][z]
            return np.interp(z, zvals, cc[:,2])
        else:
            #return integrate.quad(lambda z_:1.0/Ez(OmegaM, OmegaV, z_), 0.0,z)[0]
            return 3.6
            #return np.interp(z, zvals2, cc2[:,2])

    def phi_v(z):
        if z < z_eq:
        #return 4.0* np.pi*(phi_r(z))**2.0/(1+z)**3.0/Ez(OmegaM, OmegaV, z)
            #return cc[:,2][z]
            return np.interp(z, zvals, cc[:,3])
        else:
            #return 4.0* np.pi*(phi_r(z))**2.0/(1+z)**3.0/Ez(OmegaM, OmegaV, z)
            return 325*(z_eq)**(1.0/2.0)/ z**5.0
            #return np.interp(z, zvals2, cc2[:,3])

    def phi_n(z):
        return z**3.0 / (1.0+z)**(7.0/6.0)*(1+ (z/z_eq))**(11.0/6.0)
            

    integrand = lambda z,f: c_z(z)*phi_v(z)*cut(z,f)*cut2(z,f)/(phi_r(z))**2.0/(1.0+z)**(7.0/3.0)/(phi_t(z))**(10.0/3.0)

    zsup = 10e28
    zmaxfunc = lambda z,f:abs((f*(1.0+z)*alpha*phi_t(z)/H0)**(-1.0/3.0)-1.0)
    zmax = lambda f: fminbound(zmaxfunc,0.0,zsup, args = (f,))
    #zmax = lambda f: minimize_scalar(zmaxfunc, bounds = (0.0, zsup), method = 'bounded', args = (f,))

    Omega = np.zeros(len(f))
    for x in range(len(f)):
        if zmax(f[x]) < z_eq:
            result = integrate.quad(integrand, 0.0, zmax(f[x]), args = f[x])[0]
        else:
            result = integrate.quad(integrand, 0.0, z_eq, args = f[x])[0] + integrate.quad(integrand, z_eq, zmax(f[x]), args = f[x])[0]
        Omega[x] = constant*result *f[x]**(-1.0/3.0)
    return Omega
    
def cosmic_string_small_kinks(f,Gu,epsilon, p): 
    
    gamma = 50.0
    t0 = 10**17.5
    alpha = epsilon * gamma * Gu 
    c = 1
    constant = 4*c*Gu*(np.pi)**2.0/3.0/p*H0**(2.0/3.0)/gamma/alpha**(2.0/3.0)
    z_eq = 5440.0 
    t_eq=8.122570474611143e+11

    Ez = lambda OmegaM, OmegaV,z:(OmegaM*(1.0+z)**3.0+OmegaV)**0.5

    c_z = lambda z:1.0+ (9.0*z/(z+z_eq))

    def cut(z,f):
        if 1.0-(f*(1.0+z)*alpha*phi_t(z)/H0)**(-1.0/3.0) < 0.0:
           return 0.0
        else:
           return 1.0

    def cut2(z,f):
        if 100.0/t0/alpha**(8.0/3.0)/(f*t0)**(2.0/3.0)*phi_n(z)/f < 1.0:
            return 0.0
        else:
            return 1.0

    cc = np.loadtxt('/home/fitzaxen/Workspace/spectrum/cosmo_approx.txt')
    zvals = cc[:,0]

    def phi_t(z):
        if z < z_eq:
            #return integrate.quad(lambda z_:1.0/(1.0+z_)/Ez(OmegaM, OmegaV, z_),z,np.inf)[0]
            #return cc[:,0][z]
            return np.interp(z, zvals, cc[:,1])
        else:
            return (z_eq)**(1.0/2.0)* z**-2.0


    def phi_r(z):
        if z < z_eq:
            #return integrate.quad(lambda z_:1.0/Ez(OmegaM, OmegaV, z_), 0.0,z)[0]
            #return cc[:,1][z]
            return np.interp(z, zvals, cc[:,2])
        else:
            return 3.6

    def phi_v(z):
        if z < z_eq:
            #return 4.0* np.pi*(phi_r(z))**2.0/(1+z)**3.0/Ez(OmegaM, OmegaV, z)
            #return cc[:,2][z]
            return np.interp(z, zvals, cc[:,3])
        else:
            return 325*(z_eq)**(1.0/2.0)/ z**5.0

    def phi_n(z):
        return z**3.0 / (1.0+z)**(7.0/6.0)*(1+ (z/z_eq))**(11.0/6.0)
            

    integrand = lambda z,f: c_z(z)*phi_v(z)*cut(z,f)*cut2(z,f)/(phi_r(z))**2.0/(1.0+z)**(8.0/3.0)/(phi_t(z))**(11.0/3.0)

    zsup = 10e28
    zmaxfunc = lambda z,f:abs((f*(1.0+z)*alpha*phi_t(z)/H0)**(-1.0/3.0)-1.0)
    zmax = lambda f: fminbound(zmaxfunc,0.0,zsup, args = (f,))
    #zmax = lambda f: minimize_scalar(zmaxfunc, bounds = (0.0, zsup), method = 'bounded', args = (f,))

    Omega = np.zeros(len(f))
    for x in range(len(f)):
        if zmax(f[x]) < z_eq:
            result = integrate.quad(integrand, 0.0, zmax(f[x]), args = f[x])[0]
        else:
            result = integrate.quad(integrand, 0.0, z_eq, args = f[x])[0] + integrate.quad(integrand, z_eq, zmax(f[x]), args = f[x])[0]
        Omega[x] = constant*result *f[x]**(-2.0/3.0)
    return Omega

def cosmic_string_large_cusps(f,Gu,p): 

    gamma = 50.0
    alpha = 0.1 
    c = 1
    constant = 4.0*np.pi**2.0/3.0/H0**2.0
    #z_eq = 5440.0 
    #t_eq=8.122570474611143e+11
    z_eq = 3369
    t_eq = 1.63919e+12

    Ez = lambda OmegaM, OmegaV,z:(OmegaM*(1.0+z)**3.0+OmegaV)**0.5

    t = lambda z: phi_t(z) /H0

    def cut(z,f,l):
        if 1.0-(f*(1.0+z)*l)**(-1.0/3.0) < 0.0:
           return 0.0
        else:
           return 1.0

    def cut2(z,l,phit):
        if l > alpha*phit/H0:
           return 0.0
        else:
           return 1.0

    cc0 =  np.loadtxt('/home/fitzaxen/Workspace/spectrum/cosmo_approx0.txt')
    zvals0 = cc0[:,0]
    cc = np.loadtxt('/home/fitzaxen/Workspace/spectrum/cosmo_approx.txt')
    zvals = cc[:,0]
    cc2 = np.loadtxt('/home/fitzaxen/Workspace/spectrum/cosmo_approx2.txt')
    zvals2 = cc2[:,0]
    cc3 = np.loadtxt('/home/fitzaxen/Workspace/spectrum/cosmo_approx3.txt')
    zvals3 = cc3[:,0]

    def phi_t(zrange):
        #if z < z_eq:
         #   return integrate.quad(lambda z_:1.0/(1.0+z_)/Ez(OmegaM, OmegaV, z_),z,np.inf)[0]
            #return cc[:,1][int(z)]
          #  return np.interp(z, zvals, cc[:,1])
        #else:
         #   return (z_eq)**(1.0/2.0)* z**-2.0
        #uppervals = [integrate.quad(lambda z_:1.0/(1.0+z_)/Ez(OmegaM, OmegaV, z_),z,np.inf)[0] for z in upper]
        #uppervals = np.interp(upper, zvals2, cc2[:,1])

        #lowest = [z for z in zrange if z <= 1]
        #lower = [z for z in zrange if z <= z_eq and z > 1]
        #upper = [z for z in zrange if z > z_eq]
        #lowestvals = np.interp(lowest, zvals0, cc0[:,1])
        #lowervals = np.interp(lower, zvals, cc[:,1])
        #uppervals = np.interp(upper, zvals3, cc3[:,1])

        #uppervals = []
        #prevt = lowervals[len(lowervals)-1]
        #zvalprev = lower[len(lower)-1]
        #for n in range(0, len(upper)):
         #   uppervals.append(prevt -integrate.quad(lambda z_:1.0/(1.0+z_)/Ez(OmegaM, OmegaV, z_),zvalprev,upper[n])[0] )
          #  prevt = uppervals[n]
           # zvalprev = upper[n]
        #uppervals = [(z_eq)**(1.0/2.0)* z**-2.0 for z in upper]
        #uppervals = np.interp(upper, zvals3, cc3[:,1])
        #uppervals = [integrate.quad(lambda z_:1.0/(1.0+z_)/Ez(OmegaM, OmegaV, z_),z,np.inf)[0] for z in upper]
        #uppervals = [integrate.quad(lambda z_:1.0/(1.0+10**z_)/Ez(OmegaM, OmegaV,10** z_)*10**z_*np.log(10),np.log10(z),50.0)[0] for z in upper]
        #return cc3[:,1]
        
        return np.interp(zrange,zvals3, cc3[:,1]) 
        #return np.concatenate([lowestvals,lowervals,uppervals])

    def phi_r(zrange):
        #if z < z_eq:
            #return integrate.quad(lambda z_:1.0/Ez(OmegaM, OmegaV, z_), 0.0,z)[0]
            #return cc[:,2][int(z)]
            #return np.interp(z, zvals, cc[:,2])
        #else:
         #   return 3.6

        #lowest = [z for z in zrange if z <= 1]
        #lower = [z for z in zrange if z <= z_eq and z > 1]
        #upper = [z for z in zrange if z > z_eq]
        #lowestvals = np.interp(lowest, zvals0, cc0[:,2])
        #lowervals = np.interp(lower, zvals, cc[:,2])
        #uppervals = np.interp(upper, zvals3, cc3[:,2])

        #uppervals = []
        #prevr = lowervals[len(lowervals)-1]
        #zvalprev = lower[len(lower)-1]
        #for n in range(0, len(upper)):
         #   uppervals.append(prevr -integrate.quad(lambda z_:1.0/(1.0+z_)/Ez(OmegaM, OmegaV, z_),zvalprev,upper[n])[0] )
          #  prevr = uppervals[n]
           # zvalprev = upper[n]

        #uppervals = [3.6 for z in upper]
        #uppervals = np.interp(upper, zvals3, cc3[:,2])
        #uppervals = np.interp(upper, zvals2, cc2[:,2])
        #uppervals = [integrate.quad(lambda z_:1.0/Ez(OmegaM, OmegaV, z_), 0.0,z)[0] for z in upper]
        #uppervals = [integrate.quad(lambda z_:1.0/Ez(OmegaM, OmegaV, 10**z_)*10**z*np.log(10),-50,np.log10(z))[0] for z in upper]
        
        return np.interp(zrange,zvals3, cc3[:,2]) 
        #return np.concatenate([lowestvals,lowervals,uppervals])
        #return cc3[:,2]

    def phi_v(zrange):
        #if z < z_eq:
            #return 4.0* np.pi*(phi_r(z))**2.0/(1+z)**3.0/Ez(OmegaM, OmegaV, z)
            #return cc[:,3][int(z)]
         #   return np.interp(z, zvals, cc[:,3])
        #else:
         #   return 325*(z_eq)**(1.0/2.0)/ z**5.0

        #lowest = [z for z in zrange if z <= 1]
        #lower = [z for z in zrange if z <= z_eq and z > 1]
        #upper = [z for z in zrange if z > z_eq]
        #lowestvals = np.interp(lowest, zvals0, cc0[:,3])
        #lowervals = np.interp(lower, zvals, cc[:,3])
        #uppervals = np.interp(upper, zvals3, cc3[:,3])

        #phirs = phi_r(zrange)
        #uppervals = [325*(z_eq)**(1.0/2.0)* z**-5.0 for z in upper]
        #uppervals = np.interp(upper, zvals3, cc3[:,3])
        #uppervals = np.interp(upper, zvals2, cc2[:,3])
        #uppervals = []
        #for index, z in enumerate(upper):
         #   uppervals.append(4.0* np.pi*(phirvals[index])**2.0/(1+z)**3.0/Ez(OmegaM, OmegaV, z))
        
        return np.interp(zrange,zvals3, cc3[:,3]); 
        #return np.concatenate([lowestvals,lowervals,uppervals])
        #return cc3[:,3]

    Zi_r = 0.4*15.0*alpha**(1.0/2.0)
    Zi_m = 0.12 * 4.0

    def n(l, z,phit):
        if (phit/H0) <= t_eq:
            return Zi_r*(phit/H0)**(-3.0/2.0)*(l+(gamma*Gu/H0*phit))**(-5.0/2.0)/p
        else:
            if alpha*t_eq-gamma*Gu*((phit/H0)-t_eq) < l:
                return Zi_m *(phit/H0)**(-2.0)*(l+(gamma*Gu/H0*phit))**(-2.0)/p
            else:
                return Zi_r * (t_eq)**(1.0/2.0)*(phit/H0)**-2.0*(l+(gamma*Gu/H0*phit))**(-5.0/2.0)/p

    zsup = 10e28

    theta = lambda z, f, l: 1.0/(f*l*(1.0+z))**(1.0/3.0)
    delta = lambda z, f, l: (theta(z,f,l))**2.0/4.0
    #h = lambda l,z,f: Gu*H0*l**(2.0/3.0)/(1.0+z)**(1.0/3.0)/phi_r(z)/f**(4.0/3.0)
    h = lambda l, z, f, phir:Gu*H0*l**(2.0/3.0)/(1.0+z)**(1.0/3.0)/phir/f**(4.0/3.0)
    #dRdzdl = lambda z,f,l:phi_v(z)/H0**3.0/(1.0+z)*2.0*n(l,z,f)/l*delta(z,f,l)
    dRdzdl = lambda z,f,l, phiv,phit: phiv/H0**3.0/(1.0+z)*2.0*n(l,z,phit)/l*delta(z,f,l)
    integrand = lambda l,z,f, phit, phir, phiv:dRdzdl(z,f,l,phiv,phit)*(h(l,z,f,phir))**2.0*cut(z,f,l)*cut2(z,l,phit)

    theta2 = lambda z, l: 1.0/(l*(1.0+z))**(1.0/3.0)
    delta2 = lambda z, l: (theta2(z,l))**2.0/4.0
    #h = lambda l,z,f: Gu*H0*l**(2.0/3.0)/(1.0+z)**(1.0/3.0)/phi_r(z)/f**(4.0/3.0)
    h2 = lambda l, z, phir:Gu*H0*l**(2.0/3.0)/(1.0+z)**(1.0/3.0)/phir
    #dRdzdl = lambda z,f,l:phi_v(z)/H0**3.0/(1.0+z)*2.0*n(l,z,f)/l*delta(z,f,l)
    dRdzdl2 = lambda z,l, phiv,phit: phiv/H0**3.0/(1.0+z)*2.0*n(l,z,phit)/l*delta2(z,l)
    integrand2 = lambda l,z,phit, phir, phiv:dRdzdl2(z,l,phiv,phit)*(h2(l,z,phir))**2.0*cut2(z,l,phit)

    lmin = lambda z, f: 1.0/f/(1.0+z)
        
    lmax = lambda z,phit: alpha * phit/H0

    equalfunc = lambda z, f, phit: (1.0/f/(1.0+z))-(alpha * phit/H0)
    equal = lambda f, phit: fsolve(equalfunc, 2000.0, args = (f,phit,))

    zrange = np.logspace(-10.0,28.0, num=870)
    #print(zrange)
    #np.savetxt('zvals.txt', zrange)
    phitvals = phi_t(zrange)
    phirvals = phi_r(zrange)
    phivvals = phi_v(zrange)


    #lrange = np.logspace(-48.0, 17.0, num=500)

    #m = np.zeros((len(zrange), len(lrange)))
    #for index1, zval in enumerate( zrange):
     #   for index2,lval in enumerate(lrange):
      #      m[index1][index2] = integrand2(lval, zval, phitvals[index1], phirvals[index1], phivvals[index1])
    #np.savetxt('integrand_vals.txt', m)

    #start = time.time()
    #fil2 = np.loadtxt('integrand_vals.txt')
    #for index1, zval in enumerate(zrange):
     #   for index2, lval in enumerate(lrange):
      #      e = fil2[index1][index2]
    #end = time.time()-start
    #print(end)

    #test below f[x] = 10e-10

    #print(Gu)
    Omega = np.zeros(len(f))
    for x in range(len(f)):
        index = 0
        result = 0.0
        prevzresult = 0.0
        zvalprev = 0.0
        for index1,zval in enumerate(zrange):
            if lmin(zval, f[x]) > lmax(zval, phitvals[index]):
                result += 0.0
            else:
                lrange = np.logspace(np.log10(lmin(zval,f[x])), np.log10(lmax(zval, phitvals[index])),num=100)
                #lrange = np.logspace(np.log10(lmin(zsup, f[x])), np.log10(lmax(zsup,phitvals[199])),num=100)
            #lvalprev= lmin(zval, f[x])
            #print(lrange[0])
            #print(lrange[99])
                sum1 = 0.0
                #lvalprev= lrange[0]
                lvalprev = 0.0
                prevlresult = 0.0
                for index2,lval in enumerate(lrange):
                    #if zval < 2000:
                        #print(integrand(lval, zval, f[x], phitvals[index], phirvals[index], phivvals[index]))
                        #print(integrand2(lval, zval,phitvals[index], phirvals[index], phivvals[index])/ f[x]**(10.0/3.0)* cut(zval, f[x], lval))
                    dl = lval-lvalprev
                    val = (integrand(lval, zval, f[x], phitvals[index], phirvals[index], phivvals[index]))#*dl#+ integrand(lvalprev, zval, f[x], phitvals[index-1], phirvals[index-1], phivvals[index-1]))/2.0*dl
                    #val = (integrand2(lval, zval,phitvals[index], phirvals[index], phivvals[index])/ f[x]**(10.0/3.0)* cut(zval, f[x], lval))
                    #val = fil2[index1][index2]/ f[x]**(10.0/3.0)* cut(zval, f[x], lval)
                    if lvalprev == 0.0:
                        prevlresult = val
                    else:
                        sum1 += (prevlresult + val)/2.0 * dl
                        prevlresult = val
                    lvalprev = lval
                dz = zval - zvalprev
                if zvalprev == 0.0:
                    prevzresult = sum1
                else:
                    result += (sum1+prevzresult)/2.0 * dz
                    prevzresult = sum1
                #result += sum1*dz
            zvalprev = zval
            index += 1

        Omega[x] = constant*result*f[x]**3.0
        print('f', f[x])
        print('Omega',Omega[x])

    return Omega


   
