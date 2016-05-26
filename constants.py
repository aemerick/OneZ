__author__ = "aemerick <emerick@astro.columbia.edu>"

class constants:
    """ 
    Helpful contants. In cgs or cgs conversions 
    except ionization energies (eV)
    """

    def __init__(self):
        self.eV_erg  = 6.24150934326E11
        self.k_boltz = 1.380658E-16
        self.c       = 2.99792458E10
        self.h       = 6.6260755E-27
        self.G       = 6.6743E-8
        self.Msun = 1.998E33 # solar masses in cgs
        self.Rsun = 69.63E9  # solar radius in cgs
        self.Lsun = 3.9E33   # solar luminosity in cgs
        self.Tsun = 5777.0   # solar temperature in cgs
        self.tau_sun = 10.0E9 * 3.1536E7 # solar lifetime in cgs
        self.E_HI    = 13.6 # eV
        self.E_HeI   = 24.587

        self.Zsolar_ostar  = 0.01700 # Grevesse & Sauval 1998 - used in OSTAR2002
        self.Zsolar_parsec = 0.01524     # Caffau et. al. 2009 / 2011 - used in PARSEC SE code

        self.Zsolar_s99    = 0.02 # starburst 99 

        self.black_body_q0  = [2.89,  0.1]
        self.black_body_q1  = [5.20, 0.01]
        self.black_body_fuv = [1.0E-4, 1.0/4.35E4]
      
        self.yr_to_s = 3.16224E7


        return None


CONST = constants()
