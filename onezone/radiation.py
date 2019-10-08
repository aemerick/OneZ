# --- external ---
import numpy as np

# --- internal ---
from .constants import CONST as const


def compute_blackbody_q1(T):
    """
    Computes flux of photons at HeI line for black body of 
    temperature T
    """

    # convert to unitless energy
    x = (const.E_HeI / const.eV_erg) / (const.k_boltz * T)
    q1 = photon_radiance(x)

    # normalization constants
    A =  2.0 * const.k_boltz**3 * T**3 / (const.h**3 *const.c**2)
   
    return A * q1

def compute_blackbody_q0(T):
    """
    Computes H ionizing photon flux for a black body of temperature T
    """

    # convert to a unitless energy to compute radiance
    x = (const.E_HI / const.eV_erg) / (const.k_boltz * T)
    q0 = photon_radiance(x)

    # constants to convert back to total count
    A=  2.0 * const.k_boltz**3 * T**3 / (const.h**3 *const.c**2)

    return A * q0 


def fuv_flux_blackbody(T):
    """
    Computes black body flux over the FUV band for a 
    black body of temperature T
    """

    # convert to unitless energies
    x2 = (const.E_HI / const.eV_erg) / (const.k_boltz * T)
    x1 = (6.0        / const.eV_erg) / (const.k_boltz * T)

    # normalization const
    A = 2.0 * const.k_boltz**4 * T**4 / (const.h**3 * const.c**2)

    # flux and renorm
    fuv_flux = A * black_body_flux(x1,x2)

    return fuv_flux

def LW_flux_blackbody(T):
    """
    Computes black body flux over the LW band for a 
    black body of temperature T
    """

    # convert to unitless energies
    x2 = (const.E_HI / const.eV_erg) / (const.k_boltz * T)
    x1 = (11.2       / const.eV_erg) / (const.k_boltz * T)

    # normalization const
    A = 2.0 * const.k_boltz**4 * T**4 / (const.h**3 * const.c**2)

    # flux and renorm
    LW_flux = A * black_body_flux(x1,x2)

    return LW_flux


# AE: To do rename this function to black_body_flux and
#     have this be the main function to use. Rename the other
#     one to something else, and rewrite everything to be
#     called via this. (a little easier to use rather
#     than having to precompute the unitless numbers)
def BB_flux(E1, E2, T):
    """
    Wrapper to compute black body flux given minimum and 
    maximum energies (in eV) and a black body temperature.
    """

    x1 = E1 / (const.k_boltz * T) / const.eV_erg
    x2 = E2 / (const.k_boltz * T) / const.eV_erg

    A = 2.0 * const.k_boltz**4 * T**4 / (const.h**3 * const.c**2)

    return A * black_body_flux(x1, x2)

def black_body_flux(x1,x2):
    """ 
    Compute the black body flux between given unitless energy range
    x = (E_photon / kT) using the series approximation to compute
    the two one-sided integrals. The returned value is unitless
    and needs to be scaled by :
     2 (kT)^4 / (h^3 c^2) to get units of   energy / area / steradian
    """

    return one_sided_black_body_flux(x1) - one_sided_black_body_flux(x2)

def one_sided_black_body_flux(x):
    """
    Compute the one sided black body flux intergral between
    x = (E_photon / kT) using the series approximation to compute
    The returned value is unitless and needs to be scaled by :
     2 (kT)^4 / (h^3 c^2) to get units of   energy / area / steradian
    """

    max_iter = 513
    min_iter = 4
    tolerance = 1.0E-10

    difference = 1.0

    sum = 0.0; old_sum = 0.0
    i = 1

    while((difference > tolerance and i < max_iter) or i < min_iter):
        old_sum = sum
        sum += (x*x*x/(1.0*i) + 3.0*x*x/(1.0*i*i) + 6.0*x/(1.0*i*i*i) + 6.0/(1.0*i*i*i*i))*np.exp(-i*x)
        difference  = sum - old_sum
        i = i + 1

    return sum


def photon_radiance(x):
    """
    Integral over black body spectrum to compute the number of photons
    with energies above the given unitless energy x, where x = Eh / kT. 
    Uses series approximation of integral
    """
    max_iter = 513
    min_iter = 4
    tolerance = 1.0E-10

    difference = 1.0E10

    sum = 0.0; old_sum = 0.0
    i = 1
    while((difference > tolerance and i < max_iter) or i < min_iter):
        old_sum = sum
        sum += ( x*x/(1.0*i) + 2.0*x/(1.0*i*i) + 2.0/(1.0*i*i*i))*np.exp(-i*x)
        difference = sum - old_sum
        i = i + 1

    return sum

def average_energy(E_i, T):
    """
    Given an energy in erg, computes the average energy of photons 
    above E_i for a black body of temperature T using a series approx of
    the integral over the black body spectrum
    """

    max_iter = 513
    min_iter = 4
    tolerance = 1.0E-10
    x = E_i / (const.k_boltz * T)


    difference = 1.0E10;
    sum = old_sum = 0.0;
    i = 1;
    while((difference > tolerance and i < max_iter) or i < min_iter):
        old_sum = sum
        sum += ( x*x*x/(1.0* i) + 3.0*x*x/(1.0* i*i)
                     + 6.0*x/(1.0* i*i*i) + 6.0 / (1.0*i*i*i*i)) * np.exp(-i*x)
        difference = sum - old_sum
        i = i + 1

    u_dens_sum = sum

    sum = 0.0; old_sum = 0.0; i = 1; difference = 1.0E10

    while( (difference > tolerance and i < max_iter) or i < min_iter):
        old_sum = sum
        sum += ( x*x/(1.0*i) + 2.0*x/(1.0*i*i) + 2.0/(1.0*i**3))*np.exp(-i*x)
        difference = sum - old_sum
        i = i + 1

    return (const.k_boltz * T)*(u_dens_sum / sum)

