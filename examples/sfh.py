"""
    Author : Andrew Emerick


    Runs a single cluster initialzed at t = 0 through
    its full evolution, tracking the abundances over time.

    Most of this script sets up the configuration parameters
    (some of which are defaults anyway, but given here again).
   
    For 2 x 10^5 solar masses of stars, this should take around
    10 - 20 minutes depending on the computer and desired output
    frequency
"""

from onezone import zone
from onezone import imf
from onezone import config

#
# Final time in Myr
#
config.zone.t_final = 14000.0

#
# species to trace (m_tot and m_metal must be listed)
# All elements through Bismuth are available as single isotopes
#
config.zone.species_to_track = ['m_tot', 'm_metal', 'H', 'He', 'C', 'N','O']

#
# I/o spacing for full dumps (pickles) and summary output
# in Myr
#
config.io.dt_dump            = 14000.0
config.io.dt_summary         = 2.0
config.io.cycle_dump         = 0
config.io.cycle_summary      = 0

#
# dt = 1/32 of smallest stellar lifetime
#
config.zone.adaptive_timestep      = True
config.zone.timestep_safety_factor = 32


#
# Gas properties
#
config.zone.initial_gas_mass         = 2.0E6
config.zone.initial_metallicity      = 0.0004 # lowest posssible metallicity value
config.zone.initial_stellar_mass     = 0.0    # 8.0E5 * 2.5   ### 2.0E5


#
# Set up IMF
#
config.zone.imf = imf.salpeter(M_min = 1.0, M_max = 100.0, alpha = 1.35)

#
# Include any inflow or outflow.
#    inflow is factor of outflow
#    outflow is factor of SFR
#
config.zone.inflow_factor       = 0.0
config.zone.mass_loading_factor = 100.0

#
# Form stars at beginning of simulation
# without any further SF
#
config.zone.star_formation_method = 2
config.zone.SFH_filename          = "fiducial_SFR.dat"

#
# Physics to include in tracking yields
#
config.stars.use_stellar_winds = True
config.stars.use_snII          = True
config.stars.use_snIa          = True

# --------------------------------------------------------------



#
# Initialize zone, set initial abundances to primordial
# (minus the initial metallicity) and evolve
#

sim = zone.Zone()

sim.set_initial_abundances(config.zone.species_to_track)

sim.evolve()
