"""
   
    Author : A. Emerick
    Date   : May 2016
    
    Purpose:

"""

__author__ = "aemerick <emerick@astro.columbia.edu>"

# external
import numpy as np
from collections import OrderedDict
import pickle 

# internal
import imf  as imf
import star as star
from constants import CONST as const

class Zone:

    def __init__(self, M_gas = None, M_DM = None):

        # important values and stars
        self.M_gas = M_gas
        self.M_DM  = M_DM
        self.all_stars = star.StarList()
        self.Z         = 0.001

        self.initial_abundances = OrderedDict()
        self.species_masses     = OrderedDict()

        self.imf = imf.salpeter()

        # some private things
        self._t_last_dump        = 0.0
        self._t_last_summary     = 0.0
        self._cycle_number       = 0
        self._cycle_last_dump    = 0
        self._cycle_last_summary = 0
        self._output_number      = 0
        self._summary_output_number = 0
        self._summary_filename = 'summary_output.txt'


        self.t  = 0.0
        self.dt = 1.0 # Myr

        self.parameters = {}

        self._set_default_parameters()
        self._summary_data = {}
        self.Mdot_ej = 0.0
        self.Mdot_ej_masses = OrderedDict()
        self.SN_ej_abundances   = OrderedDict()

        self.N_SNIa = 0
        self.N_SNII = 0

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
            if e in abundances.keys():
                self.initial_abundances[e] = abundances[e]
            elif e == 'H':
                self.initial_abundances[e] = 0.75
            elif e == 'He':
                self.initial_abundances[e] = 0.25
            elif e == 'm_metal':
                self.initial_abundances[e] = self.Z
            elif e == 'm_tot':
                self.initial_abundances[e] = 1.0
            else:
                self.initial_abundances[e] = 0.0


        for e in self.initial_abundances.keys():
            self.species_masses[e] = self.M_gas * self.initial_abundances[e]

        return None
        
    def _accumulate_new_sn(self):
        star_type = self.all_stars.property_asarray('type')
        
        if np.size(star_type) > 1:
            self.N_SNII += len( star_type[ star_type == 'new_SNIa_remnant'] )
            self.N_SNIa += len( star_type[ star_type == 'new_remnant']      )
        
        return

    def evolve(self, tend):
       """
          Evolves the system
       """

       self.tend         = tend

       self._check_output() # do initial output

       while self.t < self.tend:
           # step one, evolve stars and compute m_ej
 
           self._evolve_stars()

           #
           # Count number of core collapse and SNIa 
           #
           self._accumulate_new_sn()

           # step two, compute m_sf
           self._compute_sfr()
           self._compute_outflow()
           self._compute_inflow()

           # form new stars here
           M_sf = self._make_new_stars()
           self.M_sf = M_sf

           #
           # step three compute total inflow, outflow,
           #
 
           self.M_gas += (self.Mdot_in + self.Mdot_ej -\
                          self.Mdot_out) * self.dt - M_sf

           if self.M_gas < 0:
               self.M_gas = 0.0

           # 
           # do the same for abundances
           #
           abundances = self.abundances
           for e in self.species_masses.keys():
               self.species_masses[e] += (self.Mdot_in  * self.Mdot_in_abundances(e) +\
                                          self.Mdot_ej_masses[e] +\
                                          self.Mdot_out * abundances[e]) * self.dt -\
                                          self.M_sf * abundances[e] + self.SN_ej_abundances[e]

           self._update_metallicity()
           # compute abundances

           self.t += self.dt
           self._cycle_number += 1      
 
           # check output condition
           self._check_output()

    def _update_metallicity(self):
    
        self.Z = self.species_masses['m_metal'] / self.M_gas

        return

    @property
    def abundances(self):
        abund = {}

        for x in self.species_masses.keys():
            abund[x] = self.species_masses[x] / self.M_gas

        return abund

    def Mdot_in_abundances(self, e):
        if e == 'H':
            abund = 0.75
        elif e == 'He':
            abund = 0.25
        else:
            abund = 0.0

        return abund
    
    @property
    def N_stars(self):
        return np.size(self.all_stars.stars)

    @property
    def M_stars(self):
        self._M_stars
        return np.sum(self.all_stars.M)

    def _evolve_stars(self):

        # advance each star one timestep using
        self.all_stars.evolve(self.t, self.dt)

        self.Mdot_ej = 0.0
        mass_loss_rate = self.all_stars.property_asarray('Mdot_ej') * 1.0E6 * const.yr_to_s
        self.Mdot_ej = np.sum( mass_loss_rate )

        i = 0
        


        for e in self.species_masses.keys():
            #
            # add wind ejecta abundaces from "stars" and stars that may have formed SN
            # but were alive for part of timestep.
            #
            self.Mdot_ej_masses[e]  = np.sum(self.all_stars.species_asarray('Mdot_ej_' + e, 'star') * self.all_stars.property_asarray('Mdot_ej','star'))
            self.Mdot_ej_masses[e] += np.sum(self.all_stars.species_asarray('Mdot_ej_' + e, 'new_WD')      * self.all_stars.property_asarray('Mdot_ej', 'new_WD'))
            self.Mdot_ej_masses[e] += np.sum(self.all_stars.species_asarray('Mdot_ej_' + e, 'new_remnant') * self.all_stars.property_asarray('Mdot_ej','new_remnant'))
            self.Mdot_ej_masses[e] *= 1.0E6 * const.yr_to_s

            self.SN_ej_abundances[e]   = np.sum(self.all_stars.species_asarray('SN_ej_' + e, 'new_SNIa_remnant'))
            self.SN_ej_abundances[e]  += np.sum(self.all_stars.species_asarray('SN_ej_' + e, 'new_remnant'))
            i = i + 1



    def _make_new_stars(self):
        
        # given SFR and gas resivoir
        M_sf = self.dt * self.Mdot_sf
        # sample from IMF until M_sf is reached
        star_masses = self.imf.sample(M = M_sf)

        M_sf = np.sum(star_masses)

        for m in star_masses:
            self.all_stars.append( star.Star(m, self.Z, self.abundances,
                                             tform=self.t, id=self._assign_particle_id()))

        return M_sf

    def _set_default_parameters(self):
        self.parameters['inflow_factor']  = 0.05
        self.parameters['mass_loading']   = 0.1
        self.parameters['sfr_efficiency'] = 0.01

    def _assign_particle_id(self):
        if (not hasattr(self, '_global_id_counter')):
            self._global_id_counter = 0

        num = self._global_id_counter

        self._global_id_counter += 1
        return num

    def _compute_inflow(self):
        self.Mdot_in  = self.parameters['inflow_factor'] * self.Mdot_out

    def _compute_outflow(self):
        self.Mdot_out = self.parameters['mass_loading'] * self.Mdot_sf

    def _compute_sfr(self):
        self.Mdot_sf  = self.parameters['sfr_efficiency'] * self.M_gas / (1.0) # need to fix!!! 1 = 1 Myr

    def _check_output(self):

        # check for full write out
        if( (self.t - self._t_last_dump) >= self.parameters['dt_dump'] and
            self.parameters['dt_dump'] > 0 ):
            self._t_last_dump = self.t
            self._write_full_dump()

        if( self._cycle_number == 0 or\
           ((self._cycle_number - self._cycle_last_dump) >= self.parameters['cycle_dump'])\
           and self.parameters['cycle_dump'] > 0 ): 
            self._cycle_last_dump = self._cycle_number
            self._write_full_dump()

        # now check for partial writes
        if(  self._cycle_number == 0 or\
            (self._cycle_number - self._cycle_last_summary >= self.parameters['cycle_summary'])\
            and self.parameters['cycle_summary'] > 0):

            self._cycle_last_summary = self._cycle_number
            self._write_summary_output()

        if( ((self.t - self._t_last_summary) > self.parameters['dt_summary']) and
            self.parameters['dt_summary'] > 0 ):
            self._t_last_summary = self.t
            self._write_summary_output()


        
    def _write_full_dump(self):
        # for now just pickle, but need to get fancy later

        name = "output_name_%00004i"%(self._output_number)

        print "Writing full dump output as " + name + " at time t = %4.4f"%(self.t)

        pickle.dump( self , open(name, "wb"))

        self._output_number += 1


    def _accumulate_summary_data(self):

        self._summary_data = OrderedDict()


        self._summary_data['t']       = self.t

        self._summary_data['M_gas']   = self.M_gas
        self._summary_data['M_DM']    = self.M_DM
        self._summary_data['M_star']  = np.sum(self.all_stars.M())       

        self._summary_data['Z_gas']   = self.Z
        self._summary_data['Z_star']  = np.average( self.all_stars.Z() )

        self._summary_data['N_star']  = len(self.all_stars.stars)

        self._summary_data['N_SNIa']  = self.N_SNIa
        self._summary_data['N_SNII']  = self.N_SNII

        sum_names = ['Mdot_ej', 'L_FUV', 'Q0', 'Q1']
        for n in sum_names:
            self._summary_data[n] = np.sum(self.all_stars.property_asarray(n))

        self._summary_data['L_Q0'] = np.sum(self.all_stars.property_asarray('Q0') * self.all_stars.property_asarray('E0'))
        self._summary_data['L_Q1'] = np.sum(self.all_stars.property_asarray('Q1') * self.all_stars.property_asarray('E1'))

        # now do all of the abundances
        
        for e in self.abundances.keys():

            self._summary_data[e] = self.abundances[e]

        return 

    def _write_summary_output(self):

        # probably just do this as a line in a txt file
        # with colums:
        # t mgas mdm mstars nstars mdot_sf mdot_in mdot_out mdot_ej metallicity L_wind L_SN L_fuv L_q1 L_q1 species_fracs 
        # where the luminosities are summed over all star particles from stellar wind
        # supernova and radiation contributions

        self._accumulate_summary_data()

        ncol = np.size(self._summary_data.keys())

        

        if self._summary_output_number == 0:
            header = " " + " ".join(self._summary_data.keys()) + "\n"

            f = open(self._summary_filename,'w')
            f.write(header)
        else:
            f = open(self._summary_filename, 'a')


        fmt = "%5.5E "*ncol + "\n"
        f.write(fmt%(  tuple(self._summary_data.values()) ))

        self._summary_output_number += 1
        self._summary_data.clear()
        f.close()


 
