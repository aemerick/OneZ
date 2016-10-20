"""
   
    Author : A. Emerick
    Date   : May 2016
    
    Purpose:

"""

__author__ = "aemerick <emerick@astro.columbia.edu>"

# external
import numpy as np
from collections import OrderedDict

try:
    import cPickle as pickle
except:
    print "WARNING: cPickle unavailable, using pickle - data dumps will be slow"
    import pickle

# internal
import imf  as imf
import star as star
import config as config
from constants import CONST as const



class Zone:
    """
    Zone Class

    Single zone gas reservoir attached to star formation and 
    chemical enrichment models. Parameters for a given simulation
    can be set using config. Zone initializes a simulation using
    initial conditions parameters set in config. 

    Once initialized and abundances are set, simulation can be 
    evolved:
        >>> sim = Zone()
        >>> sim.set_initial_abundances(list_of_element_names)
        >>> sim.evolve()
 
    I/O controlled by config parameters, but can be done manually
    with a full (pickle) dump or an output of summary statistics:

        >>> sim.write_full_dump()
        >>> sim.write_summary_output()
    """


    def __init__(self):

        #
        # important values and stars
        #
        self.M_gas     = config.zone.initial_gas_mass
        self.M_DM      = config.zone.initial_dark_matter_mass
        self.all_stars = star.StarList()
        self.Z         = config.zone.initial_metallicity
        self._M_sf_reservoir = 0.0

        self.initial_abundances = config.zone.initial_abundances
        self.species_masses     = OrderedDict()

        #
        # some private things
        #
        self._t_last_dump           = 0.0
        self._t_last_summary        = 0.0
        self._cycle_number          = 0
        self._cycle_last_dump       = 0
        self._cycle_last_summary    = 0
        self._output_number         = 0
        self._summary_output_number = 0

        self.t  = config.zone.t_o
        self.dt = config.zone.dt

        self._summary_data = {}
        self.Mdot_ej = 0.0
        self.Mdot_DM = 0.0
        self.Mdot_ej_masses = OrderedDict()
        self.SN_ej_masses   = OrderedDict()

        self.N_SNIa = 0
        self.N_SNII = 0

        # 
        # Create stars if starting with initial cluster
        #
        for e in config.zone.species_to_track:
            self.species_masses[e] = 0.0

        if (config.zone.initial_stellar_mass > 0.0):
            self._make_new_stars( M_sf = config.zone.initial_stellar_mass )

        self._compute_dt()

        self._update_globals()

        return 

    def set_initial_abundances(self, elements, abundances = None):
        """
        Set initial abundances using a list of desired elements and 
        associated abundances. If no abundances are provided, all are
        set to zero except H and He. Abundances dict does not have to be 
        complete
        """
        self.initial_abundances = OrderedDict()

        if abundances == None:
            abundances = {'empty' : 0.0}

        for e in elements:
            if e in abundances.iterkeys():
                self.initial_abundances[e] = abundances[e]
            elif e == 'H':
                self.initial_abundances[e] = 0.75*(1.0 - self.Z)
            elif e == 'He':
                self.initial_abundances[e] = 0.25*(1.0 - self.Z)
            elif e == 'm_metal':
                self.initial_abundances[e] = self.Z
            elif e == 'm_tot':
                self.initial_abundances[e] = 1.0
            else:
                self.initial_abundances[e] = 0.0

        for e in self.initial_abundances.iterkeys():
            self.species_masses[e] = self.M_gas * self.initial_abundances[e]


        #
        # One day, set this as list with second list of conditionals
        # so one can (at runtime) arbitrarily track on any condition
        #
        self.special_mass_accumulator = OrderedDict()
        self.special_mass_accumulator['m_massive'] = 0.0

        return None
        
    def _accumulate_new_sn(self):
        """
        Looks through all stars, checking to see if any new SN or SNIa
        occured in current timestep. Adds these to the counters
        """


        star_type = list(self.all_stars.property_asarray('type'))

        #if np.size(star_type) > 1:
        self.N_SNIa += star_type.count('new_SNIa_remnant')
        self.N_SNII += star_type.count('new_remnant')
        
        return

    def evolve(self):
        """
           Evolves the system until the end time assigned in config,
           using a constant timestep. This handles the formation of new
           stars, gas inflow and outflow, enrichment from stars, and
           outputs.
        """

        while self.t < config.zone.t_final:
            self._compute_dt()

            #
            # Check if output conditions are met
            #
            self._check_output()

            #
            # I) Evolve stars, computing ejecta rates and abundances
            #
 
            self._evolve_stars()

            #
            # II) Sum up and store number of supernovae
            #
            self._accumulate_new_sn()

            #
            # III) Compute SFR and make new stars
            #
            self._compute_sfr()
            self.M_sf = self._make_new_stars()

            #
            # IV) Compute inflow and outflow
            #
            self._compute_outflow()
            self._compute_inflow()
            self._compute_mdot_dm()

            #
            # V) Add/remove gas from zone due to inflow,
            #    outflow, SF, and stellar ejecta
            #
            abundances = self.abundances

            temp1 = abundances['m_metal'] * 1.0

            new_gas_mass =  self.M_gas + (self.Mdot_in + self.Mdot_ej_masses['m_tot'] -\
                           self.Mdot_out) * self.dt - self.M_sf +\
                           self.SN_ej_masses['m_tot']

            
            #
            # VI) Check if reservoir is empty
            # 
            if self.M_gas <= 0:
                self.M_gas = 0.0
                print "Gas in zone depleted. Ending simulation"
                break

            # 
            # VII) Compute increase / decrease of individual abundances
            #
            for e in self.species_masses.iterkeys():
                self.species_masses[e] = self.species_masses[e] +  (self.Mdot_in  * self.Mdot_in_abundances(e) +\
                                           self.Mdot_ej_masses[e] -\
                                           self.Mdot_out * abundances[e]) * self.dt -\
                                           self.M_sf * abundances[e] + self.SN_ej_masses[e]
            self.M_gas = new_gas_mass
            self.M_DM  = self.M_DM + self.Mdot_DM * self.dt

            #
            # VII) i) ensure metallicity is consistent with new abundances
            #
            self._update_metallicity()

            #
            # VIII) End of evolution, increment counters
            #
            self.t += self.dt
            self._cycle_number += 1
            self._update_globals()


        #
        # At end of simulation, force summary and dump 
        # outputs
        #
        self._check_output(force=True)

        return

    def _update_globals(self):

        config.global_values.time = self.t

        return

    @property
    def current_redshift(self):

        t_h = config.units.hubble_time
        z   = (2.0 * t_h / (3.0*self.t))**(2.0/3.0) - 1.0

        return z

    def _update_metallicity(self):
 
        self.Z = self.species_masses['m_metal'] / self.M_gas

        return

    @property
    def abundances(self):
        """
        Returns dictionary of abundances
        """

        abund = {}

        for x in self.species_masses.iterkeys():
            abund[x] = self.species_masses[x] / self.M_gas


        return abund

    def Mdot_in_abundances(self, e):
        """
        Inflow abundances
        """

        if e == 'H':
            abund = 0.75
        elif e == 'He':
            abund = 0.25
        elif e == 'm_tot':
            abund = 1.0
        else:
            abund = 0.0

        return abund
    
    @property
    def N_stars(self):
        return self.all_stars.N_stars

    @property
    def M_stars(self):
        return np.sum(self.all_stars.M)

    def _compute_dt(self):
        
        if config.zone.adaptive_timestep:
            lifetimes = self.all_stars.property_asarray('lifetime','star')

            if np.size(lifetimes) > 1:
                min_lifetime = np.min( lifetimes ) / (config.units.time)
                self.dt      = min_lifetime / (1.0 * config.zone.timestep_safety_factor)

        return

    def _evolve_stars(self):
        """
        Evolve stars using star list methods, and compute the total
        amount of ejecta from winds and supernova
        """

        #
        # zero mass accumulators before evolving
        # 
        for key in self.species_masses.iterkeys():
            self.Mdot_ej_masses[key] = 0.0
            self.SN_ej_masses[key]   = 0.0
 
        #
        # advance each star one timestep
        # optional mass arguments mean stars add 
        # to ejecta bins during evolution (winds and SN)
        # to limit number of loops through star list
        #
        self.all_stars.evolve(self.t, self.dt, ej_masses    = self.Mdot_ej_masses,  
                                               sn_masses    = self.SN_ej_masses,
                                               special_accumulator = self.special_mass_accumulator)

        self.Mdot_ej = self.Mdot_ej_masses['m_tot'] * config.units.time
     
        for e in self.species_masses.iterkeys():
            self.Mdot_ej_masses[e] *= config.units.time

        return
        #
        # set total dM/dt from all stellar winds
        #
      #  self.Mdot_ej = 0.0
     #   mass_loss_rate = self.all_stars.property_asarray('Mdot_ej') * config.units.time
    #    self.Mdot_ej = np.sum( mass_loss_rate )

        #
        # set dM/dt for all species
        # add up SN ejecta
        #
   #     i = 0
   #     for e in self.species_masses.iterkeys():
            #
            # add wind ejecta abundaces from "stars" and stars that may have formed SN
            # but were alive for part of timestep.
            #
    #        self.Mdot_ej_masses[e]  = np.sum(self.all_stars.species_asarray('Mdot_ej_' + e, 'star') * self.all_stars.property_asarray('Mdot_ej','star'))
   #         self.Mdot_ej_masses[e] += np.sum(self.all_stars.species_asarray('Mdot_ej_' + e, 'new_WD')      * self.all_stars.property_asarray('Mdot_ej', 'new_WD'))
  #          self.Mdot_ej_masses[e] += np.sum(self.all_stars.species_asarray('Mdot_ej_' + e, 'new_remnant') * self.all_stars.property_asarray('Mdot_ej','new_remnant'))
 #           self.Mdot_ej_masses[e] *= config.units.time
#
#            self.SN_ej_masses[e]   = np.sum(self.all_stars.species_asarray('SN_ej_' + e, 'new_SNIa_remnant'))
#            self.SN_ej_masses[e]  += np.sum(self.all_stars.species_asarray('SN_ej_' + e, 'new_remnant'))
#            i = i + 1
#
#
#
#        return

    def _make_new_stars(self, M_sf = -1):
        """
        Sample IMF to make new stars. Includes methods to 
        handle low SFR's, i.e. when SFR * dt < maximum star mass
        """
        #
        # compute the amount of gas to convert
        # into stars this timestep. Otherwise, make fixed
        # mass of stars
        #
        if (M_sf < 0.0):
            M_sf = self.dt * self.Mdot_sf

        if config.zone.use_SF_mass_reservoir and M_sf > 0.0:
            # 
            # Accumulate mass into "reservoir" and wait until
            # this is surpassed to form stars
            #

            self._M_sf_reservoir += M_sf

            if (self._M_sf_reservoir > config.zone.SF_mass_reservoir_size):

                # sample from IMF until M_sf is reached
                M_sf = self._M_sf_reservoir
                self._M_sf_reservoir = 0.0 # reset counter
            else:
                M_sf = 0.0

        elif (config.zone.use_stochastic_mass_sampling and\
                M_sf < config.zone.stochastic_sample_mass) and M_sf > 0.0:
            #
            # Prevent undersampling the IMF by requiring minimum mass
            # threshold. Allow SF to happen stochastically when M_sf is
            # below this threshold
            #

            probability_of_sf = M_sf / config.zone.stochastic_sample_mass

            if (probability_of_sf > np.random.random()):
                M_sf = config.zone.stochastic_sample_mass
            else:
                M_sf = 0.0

        #
        # Make new stars if M_gas->star > 0
        # 
        if M_sf > 0.0:

            # sample from IMF and sum sampled stars
            # to get actual star formation mass 
            star_masses = config.zone.imf.sample(M = M_sf)
            M_sf = np.sum(star_masses)

            # add each new star to the star list
            for m in star_masses:
#                print self.abundances
                self.all_stars.add_new_star( star.Star(m, self.Z, self.abundances,
                                                 tform=self.t,id=self._assign_particle_id()))


        return M_sf


    def _assign_particle_id(self):
        """
        Generates unique, consecutive ID for particle
        """

        if (not hasattr(self, '_global_id_counter')):
            self._global_id_counter = 0

        num = self._global_id_counter

        self._global_id_counter += 1
        return num

    def _compute_mdot_dm(self):
        """
        For cosmological simulations, computes growth of DM halo
        """

        self.Mdot_DM = 0.0
        
        if not config.zone.cosmological_evolution:
            return

        self.Mdot_DM = 46.1 * ( self.M_DM / 1.0E12 )**(1.1) *\
                       (1.0 + 1.11 * self.current_redshift) *\
              np.sqrt( config.units.omega_matter * (1.0 + self.current_redshift)**3 + config.units.omega_lambda)

        self.Mdot_DM *= const.yr_to_s

        self.Mdot_DM /= config.units.time

        return

    def _compute_inflow(self):
        """
        Compute inflow rate, as a function of outflow
        """

        self.Mdot_in  = config.zone.inflow_factor * self.Mdot_out
        return

    def _compute_outflow(self):
        """
        Compute outflow rate, as a function of SFR
        """

        # If either of the discrete SF sampling methods are used,
        # outflow should be determined by mass of stars formed, not
        # rate

        factor = config.zone.mass_loading_factor

        if config.zone.cosmological_evolution:
            factor = factor * (1.0 + self.current_redshift)**(-config.zone.mass_loading_index/2.0)

        if config.zone.use_SF_mass_reservoir or config.zone.use_stochastic_mass_sampling:
            self.Mdot_out = config.zone.mass_loading_factor * self.M_sf / self.dt
        else:
            self.Mdot_out = config.zone.mass_loading_factor * self.Mdot_sf

        return

    @property
    def t_dyn(self):
        if not config.zone.cosmological_evolution:
            print "Error: Cannot compute cosmological dynamical time with non-cosmological simulation"
            raise NotImplementedError

        return 0.1 * config.units.hubble_time * (1.0 + self.current_redshift)**(-3.0/2.0)

    def _compute_sfr(self):
        """
        Compute SFR using method set in config
        """

        if config.zone.star_formation_method == 0:
            # no star formation
            self.Mdot_sf = 0.0

        elif config.zone.star_formation_method == 1:
            # constant, user supplied SFR
            self.Mdot_sf = config.zone.constant_SFR

        elif config.zone.star_formation_method == 2 :
            self.Mdot_sf = config.zone.SFR_efficiency * self.M_gas

        elif config.zone.star_formation_method == 3:
            self.Mdot_sf = self.zone.SFR_dyn_efficiency * self.M_gas / self.t_dyn
 

        elif config.zone.star_formation_method == 4 :
            print "SFH from file not yet implemented"
            raise NotImplementedError

    def _check_output(self, force = False):
        """
        Checks output conditions, output if any are met
        """

        if force:
            self.write_full_dump()
            self.write_summary_output()
            return

        #
        # Otherwise output based on user supplied conditions
        # 

        #
        # check for full write out
        #
        if( (self.t - self._t_last_dump) >= config.io.dt_dump and\
              config.io.dt_dump > 0 ):
            self._t_last_dump = self.t
            self.write_full_dump()

        if( self._cycle_number == 0 or\
           ((self._cycle_number - self._cycle_last_dump) >= config.io.cycle_dump )\
           and config.io.cycle_dump > 0 ): 
            self._cycle_last_dump = self._cycle_number
            self.write_full_dump()

        #
        # now check for partial (summary) writes
        #
        if(  self._cycle_number == 0 or\
            (self._cycle_number - self._cycle_last_summary >= config.io.cycle_summary)\
            and config.io.cycle_summary > 0):

            self._cycle_last_summary = self._cycle_number
            self.write_summary_output()

        if( ((self.t - self._t_last_summary) > config.io.dt_summary) and
            config.io.dt_summary > 0 ):
            self._t_last_summary = self.t
            self.write_summary_output()

        return

        
    def write_full_dump(self):
        """
        Pickle current simulation 
        """

        name = config.io.dump_output_basename + "_%00004i"%(self._output_number)

        print "Writing full dump output as " + name + " at time t = %4.4f"%(self.t)

        pickle.dump( self , open(name, "wb"))

        self._output_number += 1

        return


    def _accumulate_summary_data(self):
        """
        Set up list of all of the data to output. Sum over particle properties
        to do this.
        """

        self._summary_data = OrderedDict()

        self._summary_data['t']       = self.t

        self._summary_data['M_gas']   = self.M_gas
        self._summary_data['M_DM']    = self.M_DM
        self._summary_data['M_star']  = np.sum(self.all_stars.M())       

        self._summary_data['Z_gas']   = self.Z
        self._summary_data['Z_star']  = np.average( self.all_stars.Z() )

        self._summary_data['N_star']  = self.N_stars

        self._summary_data['N_SNIa']  = self.N_SNIa
        self._summary_data['N_SNII']  = self.N_SNII

        sum_names = ['Mdot_ej', 'L_FUV', 'L_LW', 'Q0', 'Q1']
        for n in sum_names:
            self._summary_data[n] = np.sum(self.all_stars.property_asarray(n))

        self._summary_data['L_bol']  = np.sum(self.all_stars.property_asarray('luminosity'))
        self._summary_data['L_wind'] = np.sum(self.all_stars.property_asarray('mechanical_luminosity'))

        self._summary_data['L_Q0'] = np.sum(self.all_stars.property_asarray('Q0') * self.all_stars.property_asarray('E0'))
        self._summary_data['L_Q1'] = np.sum(self.all_stars.property_asarray('Q1') * self.all_stars.property_asarray('E1'))

        # now do all of the abundances
        
        for e in self.abundances.iterkeys():
            self._summary_data[e] = self.abundances[e]
            self._summary_data[e + '_mass'] = self.species_masses[e]

        for key in self.special_mass_accumulator.iterkeys():
            self._summary_data[key] = self.special_mass_accumulator[key]

        if config.io.radiation_binned_output:
            condition_1 = {'mass': lambda x : (x.M >= 1.0)  * (x.M < 8.0) *  (x.properties['type'] == 'star')} 
            condition_2 = {'mass': lambda x : (x.M >= 8.0)  * (x.M < 16.0) * (x.properties['type'] == 'star')}
            condition_3 = {'mass': lambda x : (x.M >= 16.0) * (x.M < 24.0) * (x.properties['type'] == 'star')}
            condition_4 = {'mass': lambda x : (x.M >= 24.0) * (x.M < 1000.0) * (x.properties['type'] == 'star')}

            self._summary_data['low_mass_LQ0']   = np.sum(self.all_stars.property_asarray('Q0', subset_condition = condition_1) *\
                                                          self.all_stars.property_asarray('E0', subset_condition = condition_1))
            self._summary_data['int_mass_LQ0']   = np.sum(self.all_stars.property_asarray('Q0', subset_condition = condition_2) *\
                                                          self.all_stars.property_asarray('E0', subset_condition = condition_2))
            self._summary_data['high_mass_LQ0']  = np.sum(self.all_stars.property_asarray('Q0', subset_condition = condition_3) *\
                                                          self.all_stars.property_asarray('E0', subset_condition = condition_3))
            self._summary_data['vhigh_mass_LQ0'] = np.sum(self.all_stars.property_asarray('Q0', subset_condition = condition_4) *\
                                                          self.all_stars.property_asarray('E0', subset_condition = condition_4))

            self._summary_data['low_mass_LQ1']   = np.sum(self.all_stars.property_asarray('Q1', subset_condition = condition_1) *\
                                                          self.all_stars.property_asarray('E1', subset_condition = condition_1))
            self._summary_data['int_mass_LQ1']   = np.sum(self.all_stars.property_asarray('Q1', subset_condition = condition_2) *\
                                                          self.all_stars.property_asarray('E1', subset_condition = condition_2))
            self._summary_data['high_mass_LQ1']  = np.sum(self.all_stars.property_asarray('Q1', subset_condition = condition_3) *\
                                                          self.all_stars.property_asarray('E1', subset_condition = condition_3))
            self._summary_data['vhigh_mass_LQ1'] = np.sum(self.all_stars.property_asarray('Q1', subset_condition = condition_4) *\
                                                          self.all_stars.property_asarray('E1', subset_condition = condition_4))

            FUV_1 = self.all_stars.property_asarray('L_FUV', subset_condition = condition_1)
            FUV_2 = self.all_stars.property_asarray('L_FUV', subset_condition = condition_2)
            FUV_3 = self.all_stars.property_asarray('L_FUV', subset_condition = condition_3)
            FUV_4 = self.all_stars.property_asarray('L_FUV', subset_condition = condition_4)

            self._summary_data['low_mass_LFUV']   = np.sum(FUV_1)
            self._summary_data['int_mass_LFUV']   = np.sum(FUV_2)
            self._summary_data['high_mass_LFUV']  = np.sum(FUV_3)
            self._summary_data['vhigh_mass_LFUV'] = np.sum(FUV_4)

            LW_1 = self.all_stars.property_asarray('L_LW', subset_condition = condition_1)
            LW_2 = self.all_stars.property_asarray('L_LW', subset_condition = condition_2)
            LW_3 = self.all_stars.property_asarray('L_LW', subset_condition = condition_3)
            LW_4 = self.all_stars.property_asarray('L_LW', subset_condition = condition_4)

            self._summary_data['low_mass_LLW']   = np.sum(LW_1)
            self._summary_data['int_mass_LLW']   = np.sum(LW_2)
            self._summary_data['high_mass_LLW']  = np.sum(LW_3)
            self._summary_data['vhigh_mass_LLW'] = np.sum(LW_4)

            self._summary_data['low_mass_count'] = np.size(FUV_1)
            self._summary_data['int_mass_count'] = np.size(FUV_2)
            self._summary_data['high_mass_count'] = np.size(FUV_3)
            self._summary_data['vhigh_mass_count'] = np.size(FUV_4)

            


        return 

    def write_summary_output(self):
        """
        Write out summary output by appending to ASCII file. Filename is overwritten
        at first write out, so be careful.
        """

        self._accumulate_summary_data()

        ncol = np.size(self._summary_data.keys())

        if self._summary_output_number == 0: # print the header only once
            header = " " + " ".join(self._summary_data.keys()) + "\n"

            f = open(config.io.summary_output_filename,'w')
            f.write(header)
        else:
            f = open(config.io.summary_output_filename, 'a')


        fmt = "%5.5E "*ncol + "\n"
        f.write(fmt%(  tuple(self._summary_data.values()) ))

        self._summary_output_number += 1
        self._summary_data.clear()
        f.close()


 
