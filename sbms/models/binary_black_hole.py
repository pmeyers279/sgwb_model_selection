import numpy as np

def omega_gw_spectrum(M=8.7, lamb= 0.01,eta=0.1, zi=0.1, flow = 20, fhigh=100, df=0.25):
    """
    Parameters
    ---------
    M: 'float'
       chirp mass, in solar masses [1, 2.5]
    lamb: 'float'
       mass fraction, in solar masses [0, .1]
    eta: 'float'
       mass ratio [0, 0.25]
    zi: 'float'
       spin parameter [-0.85, 0.85]
    flow: 'float'
       low end of frequency spectrum
    fhigh: 'float'
       high end of frequency spectrum
    df: 'float'
       frequency bin width
    
     Returns                                                            
    -------                                                             
    omega_gw_f : `numpy.ndarray`                                        
        :math:`\Omega_{gw}(f)`                                          
    f : `numpy.ndarray`                                                 
        frequency array                                                
    """
    if isinstance(omg_ref, dict):
    #unpack params
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
        M = float(omg_ref['M'])
    except:
        print 'WARNING: Chirp mass not set, setting to 8.7'
        pass
    try:
        lamb = float(omg_ref['lamb'])
    except:
        print 'WARNING: Mass fraction not set, setting it to 0.01'
        pass
    try:
        eta = float(omg_ref['eta'])
    except:
        print 'WARNING: Mass ratio not set, setting it to 0.1'
        pass
    try:
        zi = float(omg_ref['zi'])
    except:
        print 'WARNING: Spin Parameter not set, setting to 0.1'
        pass
        
    f = np.arange(flow, fhigh+df, df)
    omgwf = np.zeros(f.size)
    
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

    mod1 = lambda v_:(1 + (-323.0/224.0+ 451.0*eta/168.0)*v_**2 + zi*(27.0/8.0 - 11.0*eta/6.0)*v_**3)**2
    mod2 = lambda v_: (1+ (1.4545*zi - 1.8897)*v_ + (-1.8153*zi + 1.6557)*v_**2)**2

    w1 =lambda v1, mod1, mod2:(1 / v1) * (mod1 / mod2)
    w2 = lambda v1, v2, mod1:(1 / v1)*(1/v2**(4.0/3.0))*mod1
