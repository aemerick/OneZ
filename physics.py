__author__ = "aemerick <emerick@astro.columbia.edu>"

# --- external ---
import numpy as np

# --- internal ---
from constants import CONST as const
import config      as config


# helper functions for computing physics models

def s99_wind_velocity(L, M, T, Z):
    """
    Starburt99 stellar wind velocity model which computes
    the stellar wind as a four parameter function. Wind is
    returned in cm/s
    """

    v_wind = 1.23 - 0.30 * np.log10(L/const.Lsun) + 0.55*np.log10(M) +\
                    0.64 * np.log10(T) + 0.13 * np.log10(Z / const.Zsolar_s99)

    return 10.0**(v_wind) * 1.0E5 # km/s -> cm/s

def s99_wind_mdot(L, M, T, Z):
    """
    STARBURST99 stellar wind mass ejecta rate model which
    computes the stellar wind as a four parameter function. Mdot
    is returned in Msun / s
    """

    Mdot = -24.06 + 2.45*np.log10(L/const.Lsun) - 1.10 * np.log10(M) +\
                    1.31*np.log10(T) + 0.80 * np.log10(Z/const.Zsolar_s99)

    return 10.0**(Mdot) / const.yr_to_s

def SNIa_yields( elements , return_dict = False):
    """
    Wrapper around dictionary of SNIa yields from 
    Thieleman et. al. 1986 (Table 5). All isotopes for each element
    are summed together. If return_dict is true, returns yield
    dictionary instead, and 'elements' is ignored. 'elements' can be
    a single string or list of strings, where strings are atomic symbols
    """

    # dict of wSNIa values
    
    yields_dict ={'m_tot'   : 1.2447714757,
                  'm_metal' : 1.2447714757,
                  'C'       : 5.0E-2 + 4.5E-13,
                  'N'       : 2.7E-9 + 4.4E-9,
                  'O'       : 1.3E-1 + 1.1E-10 + 1.7E-12,
                  'F'       : 2.5E-13,
                  'Ne'      : 1.8E-3 + 1.1E-8  + 2.5E-3,
                  'Na'      : 1.8E-6,
                  'Mg'      : 1.6E-6 + 5.8E-6 + 4.0E-6,
                  'Al'      : 4.4E-4,
                  'Si'      : 1.5E-1 + 3.0E-4 + 3.4E-3,
                  'P'       : 1.4E-4,
                  'S'       : 8.2E-2 + 7.2E-4 + 1.5E-3 + 2.5E-8,
                  'Cl'      : 1.2E-4 + 2.8E-5,
                  'Ar'      : 1.7E-2 + 1.2E-3,
                  'K'       : 9.9E-5 + 6.6E-6,
                  'Ca'      : 1.5E-2 + 3.6E-5 + 4.2E-8 + 1.8E-5 + 1.3E-9 + 5.7E-12,
                  'Sc'      : 1.6E-7,
                  'Ti'      : 1.9E-5 + 3.1E-7 + 2.0E-4 + 9.3E-6 + 1.6E-6,
                  'V'       : 5.0E-9 + 2.8E-5,
                  'Cr'      : 2.3E-4 + 5.2E-3 + 6.6E-4 + 3.8E-5,
                  'Mn'      : 6.7E-3,
                  'Fe'      : 9.0E-2 + 6.3E-1 + 2.2E-2 + 2.5E-4,
                  'Co'      : 7.3E-4,
                  'Ni'      : 1.3E-2 + 1.4E-2 + 2.4E-4 + 5.1E-3 + 2.6E-7,
                  'Cu'      : 2.0E-6 + 8.5E-6,
                  'Zn'      : 1.3E-5 + 1.9E-5 + 8.2E-8 + 3.5E-7 + 1.0E-9,
                  'Ga'      : 1.0E-7 + 6.1E-9,
                  'Ge'      : 8.4E-7 + 5.7E-8 + 8.1E-11 + 1.8E-8}

    zero_elements = ['H','He','Li','Be','B','As','Se','Br','Kr','Rb','Sr','Y','Zr',
                     'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I',
                     'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',
                     'Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au',
                     'Hg','Tl','Pb','Bi']                 

    for e in zero_elements:
        yields_dict[e] = 0.0

    if return_dict:
        return yields_dict


    if isinstance(elements, basestring):
        return yields_dict[elements]
    else:
        return np.asarray([ yields_dict[x] for x in elements ])


def SNIa_probability(t, t_form, lifetime, DTD_slope = 1.0, NSNIa = 0.043,
                        z = 0.0):
    """
    Delay time distribution model to calculate dP/dt for a given
    white dwarf to explode as a Type Ia supernova as a function of 
    the time since the formation of its main sequence star projenitor.
    Parameters to set are the slope of the DTD and NSNIa, or the percent
    of WD's that explode as a Type Ia supernova within a hubble time. This
    number is observationally informed, but depends on one's choice of IMF
    and the mass range of MS stars that can form Type Ia white dwarf
    projenitors
    """

    dPdt = NSNIa

    if (DTD_slope == 1.0):
        dPdt /= np.log( (config.units.hubble_time(z) + t_form) / (t_form + lifetime ))
    else:
        dPdt *= (- DTD_slope + 1.0)
        dPdt /= ( (config.units.hubble_time(z) + t_form)**(-DTD_slope + 1.0) -\
                  (t_form + lifetime)**(-DTD_slope+1.0))
    
    dPdt *= (t)**(-DTD_slope)

    return dPdt

def WD_lifetime(t, t_form, lifetime, DTD_slope = 1.0, NSNIa = 0.043, z = 0):
    """
    Delay time distribution model to c.alculate the exact time at which a given WD
    will explode as a Type Ia supernova. Time is given as its lifetime,
    t_explosion = lifetime + t_form, which will essentially be the lifetime of the
    progenitor MS star (an input) plus the time it will spend as a WD. t_explosion
    equals infinity if the star never explodes (as will happen ~90-95% of the time).
    """

    # set tabulated properties:
    npoints  = 1000
    min_time = np.log10(lifetime / 10.0)
    max_time = np.log10(config.units.hubble_time(z))
    dt       = (max_time - min_time) / (1.0 * (npoints - 1))

    time = 10.0**(min_time + dt * np.arange(npoints))

    # tabulate the probability using the SNIa_probability function
    tabulated_probability    = np.zeros(npoints)
    tabulated_probability[0] = SNIa_probability(time[0] + t_form + lifetime,
                                                t_form, lifetime, DTD_slope, NSNIa, z)

    tabulated  = SNIa_probability(time + t_form + lifetime, t_form, lifetime, DTD_slope, NSNIa, z)

    f_a = tabulated[:-1]
    f_b = tabulated[1:]
    f_ab = SNIa_probability( 0.5*(time[1:] + time[:-1]) + t_form + lifetime,
                             t_form, lifetime, DTD_slope, NSNIa, z )


    tabulated_probability[1:] = (1.0/6.0) * (time[1:] - time[0:-1])*(f_a + 4.0*f_ab + f_b)
    tabulated_probability     = np.cumsum(tabulated_probability)

    rnum = np.random.random()

    if ( rnum < tabulated_probability[0] ):
        # very highly unlikely, explode right away
        # print warning since this is so unlikely
        print "WARNING: Type Ia going off immediately after star's death"
        WD_lifetime = time[0]
    elif (rnum > tabulated_probability[-1]):
        # never explode
        WD_lifetime = 1000.0 * config.units.hubble_time(z)

    else:

        WD_lifetime = time[ np.abs(tabulated_probability - rnum).argmin() ]

    return WD_lifetime

def white_dwarf_mass(M):
    """
    Initial to final mass function to return white dwarf mass 
    as a function of the mass of its main sequence star projenitor.
    IFMF taken from Salaris et. al. 2009
    """

    if M < 4.0:
        wd_mass = 0.134 * M + 0.331
    else:
        wd_mass = 0.047 * M + 0.679


    return wd_mass
