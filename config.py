__author__ = "aemerick <emerick@astro.columbia.edu>"

# --- external ---
from collections import OrderedDict

# --- internal ---
from constants import CONST as const
import imf as imf


#
# --------- Superclass for all parameters -------
#
class _parameters:

    def __init__(self):
        pass

    def help(self):
        print self.__doc__

    def reset_parameters_to_default(self):
        self.__init__()

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

        Suggested Parameters:

        imf (function) : IMF function used to generate stars. Currently
            function must be one of the imf functional forms defined
            in the ``imf'' module, but user can supply their own.
            Default imf.salpeter()

        star_formation_method (int) : switch between star formation
            schemes:

            1) constant, uniform SFR throughout evolution
            2) SFR computed cosmologically
            3) SFH table provided using SFH_filename parameter where
               either two columns are provided, time and SFR, or 
               time and stellar mass. Column headers must be named
               appropriately as ("time" or "SFR" or "mass").

         dt (float) : time step size. Default 1.0 code units (Myr)
    """

    def __init__(self):
        self.initial_gas_mass         = 0.0
        self.initial_dark_matter_mass = 0.0
        self.initial_metallicity      = 0.0
        self.species_to_track         = OrderedDict()
        self.initial_abundances       = None

        self.imf                      = imf.salpeter()
        self.star_formation_method    = 1          # 
        self.SFH_filename             = None
        self.SFR                      = 1.0        # code mass / code time
        self.SF_accumulation_mass     = 100.0

        # - inflow, outflow, and efficiency parameters
        self.inflow_factor            = 0.05
        self.mass_loading_factor      = 0.1
        self.SFR_efficiency           = 0.01


        self._time_units              = const.yr_to_s * 1.0E6 
        self.t_o                      = 0.0             # Myr
        self.dt                       = 1.0             # Myr

    @property
    def time_units(self):
        return self._time_units

    @time_units.setter
    def time_units(self, value):
        self._time_units = value

        # assert time units here

zone = _zone_parameters()

#
# ----------------- Stars and Stellar Evolution ------------------
#
class _star_particle_parameters(_parameters):


    def __init__(self):
    
        self.SNII_mass_threshold           = 8.0    
        self.SNIa_candidate_mass_bounds    = [3.0, 8.0]

        self.DTD_slope                     = 1.0
        self.NSNIa                         = 0.043

        self.use_AGB_wind_phase            = True
        self.AGB_wind_phase_mass_threshold = 8.0
        

        self.normalize_black_body_to_OSTAR = True
        self.black_body_correction_mass    = 20.0
        self.black_body_q0_factors         = const.black_body_q0
        self.black_body_q1_factors         = const.black_body_q1
        self.black_body_FUV_factors        = const.black_body_fuv


stars = _star_particle_parameters()
#
# ----------------- Input and Output --------------
#
class _io_parameters(_parameters):

    def __init__(self):
        self.dump_output_basename     = 'dump'
        self.dt_dump                  = 0.0
        self.cycle_dump               = 0

        self.summary_output_filename  = 'summary_output.txt' 
        self.dt_summary               = 0.0
        self.cycle_summary            = 0

io  = _io_parameters()

#
# ------------- Helper Functions -------------
#

def information():
    """
    Welcome to the configuration parameters for the onezone
    chemical enrichment model. Parameters are classified by
    whether or not they belong to the more global onezone gas
    resivoir, the stars (and stellar physics) itself, or 
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

    print information.__doc__

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

    print "ALL parameters reset to default"

    return
