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
        self.M_gas = M_gas
        self.M_DM  = M_DM
        self.stars = []
        self.Z     = 0.001

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

    def evolve(self, tend):
       """
          Evolves the system
       """

       self.tend         = tend

       self._check_output() # do initial output

       while self.t < self.tend:
           # step one, evolve stars and compute m_ej
 
           self._evolve_stars()

           # step two, compute m_sf
           self._compute_sfr()
           self._compute_outflow()
           self._compute_inflow()

           # form new stars here
           M_sf = self._make_new_stars()

           # step three compute inflow, outflow,

           self.M_gas += (self.Mdot_in + self.Mdot_ej -\
                             self.Mdot_out) * self.dt - M_sf

           if self.M_gas < 0:
               self.M_gas = 0.0

           self.t += self.dt
           self._cycle_number += 1      
 
           # check output condition
           self._check_output()


    def _evolve_stars(self):

        # advance each star one timestep using
        map( lambda x : x.evolve(self.t, self.dt), self.stars)

        self.Mdot_ej = 0.0
        self.Mdot_ej = np.sum([x.Mdot_ej for x in self.stars])
        self.Mdot_ej = self.Mdot_ej * 1.0E6 * const.yr_to_s # msun / myr

    def _make_new_stars(self):
        
        # given SFR and gas resivoir
        M_sf = self.dt * self.Mdot_sf

        # sample from IMF until M_sf is reached
        star_masses = self.imf.sample(M = M_sf)

        M_sf = np.sum(star_masses)

        for m in star_masses:
            self.stars.append( star.Star(m, self.Z, tform = self.t,
                                               id=self._assign_particle_id()))

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

        star_parse = lambda name: [x.properties[name] for x in self.stars]
        sum_parse  = lambda name: np.sum(star_parse(name))

        self._summary_data['t']       = self.t

        self._summary_data['M_gas']   = self.M_gas
        self._summary_data['M_DM']    = self.M_DM

        self._summary_data['Z_gas']   = self.Z
        self._summary_data['Z_star']  = np.average( [x.Z for x in self.stars])

        M                             = [x.M for x in self.stars]
        self._summary_data['M_star']  = np.sum(M)
        self._summary_data['N_star']  = np.size(M)
        self._summary_data['Mdot_ej'] = np.sum([x.Mdot_ej for x in self.stars])


        self._summary_data['L_FUV'] = sum_parse('L_FUV')
        self._summary_data['Q0']    = sum_parse('Q0')
        self._summary_data['Q1']    = sum_parse('Q1')


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


 
