import numpy as np
from astropy import units as u
from onezone import constants as const
from onezone import config as config

from onezone.plots import onezone_plot_tools as ptools


def renormalize_abundances(x, e1, e2, input_norm = 'solar', output_norm = None):

    if input_norm == 'solar':
        e1_solar = const.CONST.solar_abundance[e1]
        e2_solar = const.CONST.solar_abundance[e2]

        x = 1.0*x + (e1_solar - e2_solar)

    elif input_norm == 'gas':

        if e1 == 'O' and e2 == 'H':
            x = 1.0*x - 12.0



    if output_norm == 'gas':

        if e1 == 'O' and e2 == 'H':
            x = 1.0*x + 12.0
        # else, do nothing and just return ratio

    elif output_norm == 'solar':
        e1_solar = const.CONST.solar_abundance[e1]
        e2_solar = const.CONST.solar_abundance[e2]

        x = 1.0*x - (e1_solar - e2_solar)

    elif (not (output_norm is None)) and np.size(output_norm) == 1:

        x = 1.0*x + output_norm

    elif not (output_norm is None):
        print "No output options currently exist aside from 'solar', 'gas', or None"
        raise ValueError

    return x



def abundance_ratio(x1, x2, input_type = 'abundance', normalize = 'solar'):
    """
    Return an abundance ratio with optional normalization.
    x1 and x2 are tuples
    containing either element atomic number and abundance, or
    atomic symbol and abundance (abundance = number of particles).
    E.g.:
       normalize_abundance_ratio( ('Fe',1), ('H', 100) )
    or
       normalize_abundance_ratio( (26, 1), (1, 100) )

    Optionally, x1 and x2 can be mass (in cgs, unless astropy units
    are used) or actual abundances.

    Normalization options:
       "solar" or True: Normalizes to standard stellar abundance
             notation, where returned ratio is given as:

             [x1/x2] = log10(x1/x2) - log10(x1_sun / x2_sun)

       None or False: No normalization. Just returns:
                      log10(x1/x2)

       int / float  : Supplying a single value will return:
                      log10(x1/x2) + normalize
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

    aratio = np.log10(x1_abund / x2_abund)

    norm = 0.0

    if normalize is None or (normalize == False):
        norm = 0.0

    elif (normalize == 'solar') or (normalize == True):
        x1_solar = const.CONST.solar_abundance[x1[0]]
        x2_solar = const.CONST.solar_abundance[x2[0]]

        norm =  -1.0*(x1_solar - x2_solar) # np.log10( x1_solar / x2_solar)

    elif normalize == 'gas':

        if (((x1[0] == 'O') and (x2[0] == 'H')) or\
            ((x1[0] == 8)   and (x2[0] == 1))):

            norm = 12.0

    elif len(normalize) == 1:
        # assume a number
        norm = normalize


    return aratio + norm


def elemental_abundance(element, mass):

    n = mass / (const.CONST.molecular_weight[element] * const.CONST.amu)

    return n

def IMF_weighted_quantity(Mavg, imf = None, quantity = None, Z = 0.0043,
                          M   = None,
                          n_sample = 10000):
    """
    Given an IMF, which is a function that accepts a mass and 
    returns the value of the IMF at at point, this function
    computes IMF weighted properties of stars.

    IMF : function or string
         Pass a function that takes a single argument (mass) 
         and returns value of IMF at that mass OR a name of an
         IMF, in which case a default will be loaded (not functional).
         If none provided, the config.zone IMF settings are loaded.
         Default: None
    """

    if imf is None:
        imf = config.zone.imf

    elif not hasattr(imf, '__call__'):
        # then like a string, load appropriate IMF based on name
        print "Ability to use this using IMF string name not yet implemented"
        raise RuntimeError

    available_quantities = ['Q0','Q1','E0','E1','luminosity','lifetime',
                            'age_agb','Teff','L_FUV','L_LW', 'ones']

    if (quantity is None) or \
       ((not isinstance(quantity, basestring)) and M is None):
        print """ Must provide a quantity to IMF average
                  this is either a one argument function of mass,
                  a string name of any of the below listed names,
                  or a list of values. If values are passed, one must
                  also provide a corresponding list of masses """


        print available_quantities

    elif isinstance(quantity, basestring):
        s, m, z = ptools.star_sample( n_sample, [1.0, 100.0], Z)
        s = np.reshape(s, (n_sample,))
        m = np.reshape(m, (n_sample,))
        if quantity == 'ones':
            quantity_value = np.ones(np.size(m))
        else:
            quantity_value = np.array([item.properties[quantity] for item in s])
        
        select  = (m >= Mavg[0]) * (m < Mavg[1])
        denom   = imf.imf(m)
        num     = quantity_value[select] * denom[select]

        IMF_average = np.trapz(num) / np.trapz(denom)
    
    return IMF_average
        
