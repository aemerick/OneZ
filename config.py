__author__ = "aemerick <emerick@astro.columbia.edu>"

# --- external ---
from collections import OrderedDict
import numpy as np

# --- internal ---
from .constants import CONST as const
from . import imf as imf

import os

#
# --------- Superclass for all parameters -------
#
class _parameters(object):

    def __init__(self):
        pass

    def help(self):
        print(self.__doc__)

    def reset_parameters_to_default(self):
        self.__init__()

class _globals(_parameters):
    """
    Glbal vals
    """

    def __init__(self):

        self.time = 0.0

global_values = _globals()

#
# ----------- Units ----------------------
#
class _units(_parameters):
    """
    Units:

        Set the conversions between code units and 'base'
        quantities. For the sake of the onezone code,
        the base unit of time is in seconds and the base unit
        of mass is in solar masses. Therefore, in order to
        set the time units to Myr (default) and mass units
        to solar masses (default), one would just do:

            >>> units.time = 3.1536E13
            >>> units.mass = 1.0

        Output masses are always in solar, while everything else (e.g.
        luminosity) is in cgs, regardless of what the code units are
        set to.
    """

    def __init__(self):
        self.time      = const.yr_to_s * 1.0E6
        self.mass      = 1.0

        self.omega_b = 0.0463
        self.omega_c = 0.2330
        self.omega_v = 0.7210
        self.omega_r = 1.0E-4

        self.omega_m = self.omega_c + self.omega_b
        self.omega   = self.omega_m + self.omega_v + self.omega_r

#        self.omega_matter = 0.258 # WMAP 5 year
#        self.omega_lambda = 0.742 # WMAP 5 year
        self.H_o          = 69.32 # km / s / Mpc



        return

    def H_a(self, a):
        H = self.H_o * np.sqrt((self.omega_v + self.omega_m*a**(-3) +\
                            self.omega_r*a**(-4) - (self.omega-1.0)*a**(-2)))
        return H

    def hubble_time(self, z):
        a = 1.0 / (1.0 + z)
        return (1.0 / (self.H_a(a) * (const.km) / (const.Mpc))) / (self.time)

units = _units()
#
# -------- Global Zone Parameters -------
#
class _zone_parameters(_parameters):
    """
    Zone Parameters:

        The below is a list of all parameters that are set
        by the zone, including default values. Parameters
        are listed in order of : 1) required, 2) not required but
        highly suggested, or 3) optional

        Required Parameters:

        initial_gas_mass (float) : initial gas mass in solar masses
        initial_dark_matter_mass (float) : DM mass of halo in solar
        initial_metallicity (float)      : initial gas metal fraction
        dt (float)                       : constant timestep size in code time
        t_final (float)                  : time to end simulation


        Suggested Parameters:

        imf (function) : IMF function used to generate stars. Currently
            function must be one of the imf functional forms defined
            in the ``imf'' module, but user can supply their own.
            Default imf.salpeter()

        star_formation_method (int) : switch between star formation
            schemes:

            1) constant, uniform SFR throughout evolution
            2) SFR computed based on gas mass and input SFR rate efficiency (cosmological)
            3) SFH table provided using SFH_filename parameter where
               either two columns are provided, time and SFR, or
               time and stellar mass. Column headers must be named
               appropriately as ("t" or "SFR" or "mass"). Time is assumed to be in Myr
               and SFR is assumed to be in Msun / yr (regardless of code units).

         use_SF_mass_reservoir (bool , optional) : One of two ways to deal with low
             SFR's to ensure accurate sampling of the IMF (see the second
             below). Make sure you understand these parameters and their
             implications before setting -- know also that they may need to
             be set in certain situations. This parameter turns on the
             reservoir method whereby M_sf = dt * SFR is added to a
             reservoir each timestep. If the reservoir exceeds the
             mass threshold ``SF_mass_reservoir_size'', SF occurs in that timestep
             using up all of the reservoir. This may lead to bursty, intermittent SF
             depnding on size of resivoir. Default = False

         use_stochastic_mass_sampling (bool, optional) : Second and preferred method to deal
             with low SFR's. Prevent normal SF if M_sf = dt * SFR is below some threshold
             value, ``stochastic_sample_mass''. Instead of shutting off SF that
             timestep completely, however, and instead of accumulating mass in a
             reseroir, compute the probablility that a chunk of gas of mass
             ``stochastic_sample_mass'' is made into stars that timestep as
             P = M_sf / stochastic_sample_mass . That chunk is then formed
             completely into stars using a random number draw. Default is True

         SF_mass_reservoir_size (float, optional) : Size of accumulation reservoir used with
         ``use_SF_mass_reservoir''. Default is 1000.0

         stochastic_sample_mass (float, optional) : Size of mass chunk allowed to form
             stochastically when SFR is low. Be careful setting this to too large
             a value. Not recommended to set below ~200.0 Msun depending
             on one's choice of maximum star particle mass. Default is 250.0.

         inflow_factor  (float, optional) : Sets the mass inflow rate as a function of
             the star formation rate. Default 0.05

         mass_loading_factor (float, optional) : Sets the mass outlflow rate as a function of
             the star formation rate. Default 0.1

         SFR_efficiency (float, optional) : For cosmologically derived SFR's, sets the
             star formation rate efficiency of converging gas to stars in a free fall time
             Default is 0.01

         Optional:

         t_o     (float, optional) : initial time. Default is 0.0
         t_final (float, optional) : simulation end time. Default is 10 Gyr

    """

    def __init__(self):
        self.initial_gas_mass         = 0.0
        self.initial_dark_matter_mass = 0.0
        self.initial_metallicity      = 0.0
        self.initial_stellar_mass     = 0.0
        self.species_to_track         = OrderedDict()
        self.initial_abundances       = None
        self.track_massive_star_ejecta_mass = 25.0

        self.imf                     = imf.salpeter()
        self.M_min                   = self.imf.M_min
        self.M_max                   = self.imf.M_max
        self.alpha                   = self.imf.alpha

        self.star_formation_method    = 1           # 0, 1, 2, 3, 4
        self.SFR_filename             = "SFR.in"    # t in Myr,   SFR in Msun / yr
        self.constant_SFR             = 10.0        # code mass / code time
        self.outflow_filename         = "mass_outflow.in"

        self.cosmological_evolution   = False       # on or off
        self.initial_redshift         = 0.0       #
        self.final_redshift           = 0.0
        self.current_redshift         = self.initial_redshift

        self.use_SF_mass_reservoir     = False
        self.SF_mass_reservoir_size    = 1000.0

        self.use_stochastic_mass_sampling = True
        self.stochastic_sample_mass       = 200.0

        # - inflow, outflow, and efficiency parameters
        self.mass_outflow_method      = 0               # 1 is use below, 2 is read from file
                                                        # 3: use mass loading factor for total, H, and He
                                                        #    otherwise, use constant fraction of production (read from file)
                                                        #
        self.inflow_factor            = 1.03            # ratio of inflow to outflow
        self.mass_loading_factor      = 15.7            # constant if non-cosmological
                                                        # value at z = 0 if cosmological
        self.mass_loading_index       = 3.32            # index to redshift power law
        self.SFR_efficiency           = 6.10E-4         #  1 / Myr
        self.SFR_dyn_efficiency       = 0.0950          # unitless (sf / f_dyn)

        self.t_o                      = 0.0             # Myr
        self.t_final                  = 1.0E4           # Myr
        self.dt                       = 1.0             # Myr
        self.max_dt                   = 100.0           # Myr (max dt to use when adaptive)
        self.adaptive_timestep        = True
        self.timestep_safety_factor   = 4


        self._maximum_stars      = None
        self.optimize            = True

        self.constant_metallicity = False               # Use fixed Z for gas
        self.minimum_star_particle_mass = -1            # when > 0, adds all stars < this into one bin per timestep

        # assert time units here

    @property
    def maximum_stars(self):

        if self._maximum_stars == None:
            return self._maximum_stars
        else:
            return int(self._maximum_stars)

    @maximum_stars.setter
    def maximum_stars(self, val):
        self._maximum_stars = val
        return

    @property
    def M_min(self):
        return self.__M_min

    @property
    def M_max(self):
        return self.__M_max

    @property
    def alpha(self):
        return self.__alpha

    @M_min.setter
    def M_min(self, value):
        self.__M_min = value
        self.imf.M_min = self.__M_min
        return

    @M_max.setter
    def M_max(self, value):
        self.__M_max = value
        self.imf.M_max = self.__M_max
        return

    @alpha.setter
    def alpha(self, value):
        self.__alpha = value
        self.imf.alpha = self.__alpha
        return


zone = _zone_parameters()

#
# ----------------- Stars and Stellar Evolution ------------------
#
class _star_particle_parameters(_parameters):
    """
    Star and Stellar Physics Parameters:

        The below is a list of all parameters that are set to be
        used in evolving stars and controlling the underlying stellar
        physics properties.

        SNII_mass_threshold (float) : Lower mass limit for stars to
            explode as a Type II supernovae at the end of their life.
            Default is 8.0

        SNIa_candidate_mass_bounds (list or array of floats) : Size
            two (lower and upper bound) boundaries for mass range where
            stars turn into WD's that are candidates for exploding as
            Type 1a. Default [3.0, 8.0]

        DTD_slope (float) : Slope of the delay time distribution (DTD)
            model used to compute probability of SNIa candidates
            exploding as SNIa in a given timestep. Slope is beta,
            where probability goes as t^(-beta). Default 1.0

        NSNIa (float) : Fraction of SNIa candidates that will explode
            as Type Ia supernovae in a hubble time. Default 0.043

    """

    def __init__(self):

        self.SNII_mass_threshold           = 8.0
        self.SNIa_candidate_mass_bounds    = [3.0, 8.0]

        self.DTD_slope                     = 1.0
        self.NSNIa                         = 0.043

        self.use_snII                      = True
        self.use_snIa                      = True

        self.use_stellar_winds             = True
        self.use_AGB_wind_phase            = True
        self.AGB_wind_phase_mass_threshold = 8.0
        self.AGB_wind_velocity             = 20.0      # km / s

        self.direct_collapse_mass_threshold = 25.0     # top of NuGrid data set
        self.extrapolate_snII_yields        = False
        self.use_massive_star_yields        = True


        self.normalize_black_body_to_OSTAR = True
        self.black_body_correction_mass    = 20.0
        self.black_body_q0_factors         = const.black_body_q0
        self.black_body_q1_factors         = const.black_body_q1
        self.black_body_FUV_factors        = const.black_body_fuv
        self.black_body_LW_factors         = const.black_body_LW


stars = _star_particle_parameters()
#
# ----------------- Input and Output --------------
#
class _io_parameters(_parameters):

    def __init__(self):
        self.dump_output_basename     = 'dump'
        self.dt_dump                  = 0.0
        self.cycle_dump               = 0

        self.pickle_output_basename   = 'pickle'
        self.dt_pickle                = 0
        self.cycle_pickle             = 0

        self.summary_output_filename  = 'summary_output.txt'
        self.dt_summary               = 0.0
        self.cycle_summary            = 0

        self._abundance_output_filename = None # 'abundances.dat'

        self.radiation_binned_output  = 0 # bin rad in mass bins - expensive


    @property
    def abundance_output_filename(self):
        return self._abundance_output_filename

    @abundance_output_filename.setter
    def abundance_output_filename(self, value):
        self._abundance_output_filename = value

        self._initialize_abundance_output()

        return

    def _clean_up(self):

        self._abundance_output_file.close()

        return

    def _initialize_abundance_output(self):

        if self._abundance_output_filename is None:
            return

        elif os.path.isfile(self._abundance_output_filename):
            self._abundance_output_file = open(self._abundance_output_filename,'a')
        else:
            self._abundance_output_file = open(self._abundance_output_filename,'w')

            self._abundance_output_file.write("# t id M Z lifetime")
            for e in zone.species_to_track:
                self._abundance_output_file.write(" " + (e))
            self._abundance_output_file.write("\n")

        return

io  = _io_parameters()

#
# ------------------ Data Table ----------------------
#
class _data_table(_parameters):

    def __init__(self):

        self.yields_mass_limits = [1.0, 25.0]

data = _data_table()
#
# ------------- Helper Functions -------------
#

def information():
    """
    Welcome to the configuration parameters for the onezone
    chemical enrichment model. Parameters are classified by
    whether or not they belong to the more global onezone gas
    reservoir, the stars (and stellar physics) itself, or
    input/output. The parameters can be accessed and modified
    as attributes of the following objects:

        Zone            :    zone
        Stars           :    stars
        Input / Output  :    io

    More information about these parameters can be found by
    calling the help method on a given object, e.g.:
         zone.help()
    which will print an explanation about each parameter, whether or
    not it requires user to set (i.e. will fail if defaults are kept),
    and its default value.

    Call the 'help' function in this module to see a description of
    all parameters at once.
    """

    print(information.__doc__)

def help():
    """
    Print the docstrings for all config parameter types.
    """

    for obj in [zone, stars, io]:
        obj.help()

    return

def reset_all_parameters():
    """
    Reset all parameter types to their default values
    """

    for obj in [zone, stars, io]:
        obj.reset_parameters_to_default()

    print("ALL parameters reset to default")

    return
