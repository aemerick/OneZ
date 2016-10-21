import numpy as np
from astropy import units as u
from onezone import constants as const


def normalize_abundance_ratio(x1, x2, input_type = 'abundance'):
    """
    Normalize abundance ratio to solar. x1 and x2 are tuples
    containing either element atomic number and abundance, or
    atomic symbol and abundance (abundance = number of particles).
    E.g.:
       normalize_abundance_ratio( ('Fe',1), ('H', 100) )
    or
       normalize_abundance_ratio( (26, 1), (1, 100) )

    Returns [x1/x2] where
       [x1/x2] = log10(x1/x2) - log10(x1_sun / x2_sun)

    Optionally, x1 and x2 can be mass (in cgs, unless astropy units
    are used) and 
    """

    # convert to abundance first
    x1_abund = x1[1]
    x2_abund = x2[1]

    if input_type == 'mass' or hasattr(x1[1],'unit'):

        # if we have units, use them to convert to grams
        # and strip units
        # otherwise assume grams
        if hasattr(x1[1], 'unit') and hasattr(x2[1], 'unit'):
            x1_abund = x1[1].to(u.g).value
            x2_abund = x2[1].to(u.g).value


        x1_abund = elemental_abundance(x1[0], x1[1])
        x2_abund = elemental_abundance(x2[0], x2[1])

    x1_solar = const.CONST.solar_abundance[x1[0]]
    x2_solar = const.CONST.solar_abundance[x2[0]]

    aratio = np.log10(x1_abund / x2_abund) - np.log10( x1_solar / x2_solar)

    return aratio

def elemental_abundance(element, mass):

    n = mass / (const.CONST.molecular_weight[element] * const.CONST.amu)

    return n

