# --- external ---
import numpy as np

# --- internal ---
from constants import CONST as const


def compute_blackbody_q1(T):
    x = (const.E_HeI / const.eV_erg) / (const.k_boltz * T)
    q1 = photon_radiance(x)

    A =  2.0 * const.k_boltz**3 * T**3 / (const.h**3 *const.c**2)
   
    return A * q1

def compute_blackbody_q0(T):

    x = (const.E_HI / const.eV_erg) / (const.k_boltz * T)
    q0 = photon_radiance(x)

    A=  2.0 * const.k_boltz**3 * T**3 / (const.h**3 *const.c**2)
    q0 = q0 * A

    return A * q0 


def fuv_flux_blackbody(T):
    x2 = (const.E_HI / const.eV_erg) / (const.k_boltz * T)
    x1 = (6.0        / const.eV_erg) / (const.k_boltz * T)


    A = 2.0 * const.k_boltz**4 * T**4 / (const.h**3 * const.c**2)

    fuv_flux = A * black_body_flux(x1,x2)

    return fuv_flux


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

