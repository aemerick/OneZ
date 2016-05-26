__author__ = "aemerick <emerick@astro.columbia.edu>"

# --- external ---
import numpy as np

# --- internal ---
from constants import CONST as const

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
