#
# Example useage for a parallel run
#
import numpy as np
from joblib import Parallel, delayed
import multiprocessing
import sys, glob

# onezone imports
from onezone import zone, imf, config
from onezone.constants import CONST as const

#------------------------------------
#
# Fill in onezone parameters here
#    example is VERY simple
config.units.time            = 1.0E6 * const.yr_to_s

config.zone.initial_gas_mass = 1.0E6
config.zone.constant_SFR     = 1.0E-4 / const.yr_to_s * config.units.time
config.zone.t_final          = 1000.0   # 1000 Myr
config.zone.star_formation_method = 1   # constant SFR

config.zone.max_dt                = 5.0 # Myr
config.zone.adaptive_timestep     = True

config.zone.species_to_track = ['m_tot','m_metal','H','He','C','N','O','Mg','Fe','Ba']

config.io.dt_summary             = 10.0
config.io.cycle_summary          = 0
config.io.dt_dump                = 500.0

#
# -----------------------------------
#

# ------------------- Set up and run a parallel run -------

run_basename = "run%0004i"
summary_basename = "_summary_output.txt"
istart = len(glob.glob("*_summary_output.txt"))

def run(i):
  # runs one instance of the one zone model
  np.random.seed(i) # set different seed on each process
  config.io.summary_output_filename = run_basename%(istart+i) + summary_basename
  config.io.dump_output_basename    = run_basename%(istart+i) + '_dump'
  config.io.abundance_output_filename = run_basename%(istart+i) + '_abundances.dat'

  sim = zone.Zone()
  sim.set_initial_abundances(config.zone.species_to_track)

  sim.evolve()

  del(sim)
  return


if __name__ == "__main__":

    # default run 1 simulation per core available
    n_jobs        = multiprocessing.cpu_count()
    n_simulations = n_jobs*1

    if len(sys.argv) > 1:
        n_simulations = int(sys.argv[1]))
    if len(sys.argv) > 2:
        n_jobs = int(sys.argv[2])

    Parallel(n_jobs=n_jobs)(delayed(run)(i) for i in np.arange(n_simulations))

