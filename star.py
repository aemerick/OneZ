
__author__ = "aemerick <emerick@astro.columbia.edu>"

# --- external ---
import numpy as np
import gc
from collections import OrderedDict

# --- internal ---
import data_tables as DT
import radiation   as rad
import physics     as phys

from constants import CONST as const

# ------- load the global tables -----------
SE_TABLE  = DT.StellarEvolutionData()
RAD_TABLE = DT.RadiationData()

SN_YIELD_TABLE = DT.StellarYieldsTable('SNII')
WIND_YIELD_TABLE = DT.StellarYieldsTable('wind')

class StarParticle:

    def __init__(self, M, Z, abundances=None, tform=0.0, id = 0):

        self.M   = M
        self.M_o = M
        self.Z   = Z

        self.tform = tform
        self.id    = id


        self.properties = {} # constant properties
        self.properties['type'] = 'star'

        if not abundances == None:
            self.wind_ejecta_abundances = OrderedDict()
            self.sn_ejecta_abundances = OrderedDict()

            for e in abundances.keys():
                self.wind_ejecta_abundances[e] = 0.0
                self.sn_ejecta_abundances[e]  = 0.0


        else:
            self.wind_ejecta_abundances = OrderedDict() # time varying
            self.sn_ejecta_abundances   = OrderedDict() # constant (only go SN once)

       

    def evolve(self, t, dt):
        pass

    def _assign_properties(self):
        """
        Compute stellar properties given M and Z
        """
        pass


class Star(StarParticle):
        
    def __init__(self, *args, **kwargs):
        StarParticle.__init__(self, *args, **kwargs)



        self._assign_properties()

    def evolve(self, t, dt):
        """
        Evolve the star
        """

        age = t - self.tform

        #
        # check and update Mdot_ej from stellar winds
        #
        self.Mdot_ej = 0.0
        self.stellar_wind_parameters(age, dt)
        self.Mdot_ej = self.properties['Mdot_wind']

        SN_mass_loss = 0.0

        if (age + dt > self.properties['lifetime'] / (1.0E6 * const.yr_to_s)):

            if 'new' in self.properties['type']:
                #
                # star changed types in previous timestep, update to "old"
                # and do nothing else
                #
                self.properties['type'] = self.properties['type'].replace('new_','')
                self._clear_properties()



            if self.M > 8.0:
                #
                # Core collapse supernova - change type and compute yields
                #
                self.set_SNII_properties()
                self.properties['type'] = 'new_remnant'
                SN_mass_loss = self.sn_ejecta_abundances['m_tot']
            else:
                #
                # Otherwise, form a white dwarf when dead and label as
                # candidate for future SNIa
                #
                self.M = phys.white_dwarf_mass(self.Mo)
                self.properties['type']               = 'new_WD'

                if self.M_o > 3.0 and self.M_o < 8.0:
                    self.properties['SNIa_candidate'] = True
                    self.properties['PSNIa']          = 0.0 # initialize prob of going SNIa
                else:
                    self.properties['SNIa_candidate'] = False

            # if this is a WD, need to check and see if it will explode            
            if self.properties['type'] == 'WD':
                if self.properties['SNIa_candidate']:
                    self.properties['PSNIa'] = dt * phys.SNIa_probability(self.M, t,
                                                                     self.t_form,
                                                                     self.properties['lifetime'])

                    if self.properties['PSNIa'] * dt > np.random.rand():

                        # go Type Ia supernova
                        self.properties['type'] = 'new_SNIa_remnant'
                        self.set_SNIa_properties()
                        SN_mass_loss = self.sn_ejecta_abundances['m_tot']

        #
        # Compute total mass lost through supernova and wind
        #
        M_loss = self.Mdot_ej * dt + SN_mass_loss
        
        # update mass
        self.M = self.M - M_loss

        if self.M <= 0.0:
            print "ERROR IN STAR: Negative or zero stellar mass"
            raise RuntimeError

        return None
    def set_SNIa_properties(self):
        # need to rename and combine this and functions below

        if len(self.wind_ejecta_abundances.keys()) > 0:
            yields = phys.SNIa_yields(self.wind_ejecta_abundances.keys())

            i = 0
            for e in self.sn_ejecta_abundances.keys():
                self.sn_ejecta_abundances[e] = yields[i]

        else:
            return NotImplementedError

        return yields

    def set_SNII_properties(self):

        if len(self.wind_ejecta_abundances.keys()) > 0:
            yields =  SN_YIELD_TABLE.interpolate([self.M_o, self.Z],
                                                  self.wind_ejecta_abundances.keys())

            i = 0
            for e in self.sn_ejecta_abundances.keys():
                self.sn_ejecta_abundances[e] = yields[i]

        else:
            raise NotImplementedError

    @property
    def mechanical_luminosity(self):
        return self.properties['Mdot_wind'] * const.Msun * self.properties['v_wind']**2

    def surface_gravity(self):
        return const.G * self.M * const.Msun / self.properties['R']**2

    def surface_area(self):
        return 4.0 * np.pi * self.properties['R']**2

    def _clear_properties(self):
        """
        zeroes certain properties after star dies
        """
        zeroed_properites = ['E0', 'E1', 'L_FUV', 'Q0', 'Q1', 'luminosity', 'v_wind',
                             'Mdot_wind']

        for p in zeroed_properties:
            self.properties[p] = 0.0

        return

    def ionizing_photons(self, photon_type):

        if not 'q' in photon_type:
            if photon_type == 'HI':
                photon_type = 'q0'
            elif photon_tupe == 'HeI':
                photon_type = 'q1'

        if self.properties.haskey(photon_type):
            return self.properties[photon_type] * self.surface_area()
        else:
            return 0.0

    def fuv_luminosity(self):
        if self.properties.haskey('FUV_flux'):
            return self.properties['FUV_flux'] * self.surface_area()
        else:
            return 0.0

    def stellar_wind_parameters(self, age, dt):

        if True:
            # need to compute wind velocities for all stars
            # check if star's wind is ON
            do_wind = True

            if self.M < 8.0:
                if age < self.properties['age_agb'] / ( 1.0E6 * const.yr_to_s):
                    do_wind = False
                    wind_lifetime = 0.0
                else:
                    wind_lifetime = (self.properties['lifetime'] - self.properties['age_agb'])
             
            else:
              
                wind_lifetime = self.properties['lifetime']

            if wind_lifetime < dt:
                wind_lifetime = dt

            if do_wind and age < self.properties['lifetime']:           
                Mdot   = self.properties['M_wind_total'] / wind_lifetime                
            else:
                Mdot   = 0.0

        else:

            Mdot  = phys.s99_wind_mdot(self.properties['luminosity'], self.M_o,
                                      self.properties['Teff'], self.Z)


        vwind = phys.s99_wind_velocity( self.properties['luminosity'], self.M_o,
                                        self.properties['Teff'], self.Z)


        self.properties['Mdot_wind'] = Mdot
        self.properties['v_wind']    = vwind

    def compute_stellar_wind_yields(self):
        yields = np.asarray(WIND_YIELD_TABLE.interpolate([self.M_o, self.Z],
                                                          self.wind_ejecta_abundances.keys()))


        return np.asarray(yields)


    def _assign_properties(self):

        p_list = ['luminosity', 'radius',
                  'lifetime'  , 'age_agb', 'L_FUV',
                  'Q1', 'Q2', 'E_Q1', 'E_Q2']

        # update SE table to handle arrays 
        
        L, T, R, lifetime, age_agb = SE_TABLE.interpolate([self.M,self.Z], ['L','Teff','R','lifetime','age_agb'])
        self.properties['luminosity']  = L * const.Lsun
        self.properties['Teff']        = T
        self.properties['R']           = R
        self.properties['lifetime']    = lifetime
        self.properties['age_agb']     = age_agb

        Q0, Q1, FUV = RAD_TABLE.interpolate([self.properties['Teff'],
                                             self.surface_gravity(),
                                             self.Z], ['q0','q1','FUV_flux'])
        
        E0  = rad.average_energy(const.E_HI/ const.eV_erg, self.properties['Teff'])
        E1  = rad.average_energy(const.E_HeI/const.eV_erg, self.properties['Teff'])


        use_blackbody = False
        for a in [Q0, Q1, FUV]:
            if a == 'offgrid':
                use_blackbody = True
                break;

        if use_blackbody:
            FUV = rad.fuv_flux_blackbody(self.properties['Teff'])
            Q0  = rad.compute_blackbody_q0(self.properties['Teff'])
            Q1  = rad.compute_blackbody_q1(self.properties['Teff'])

            if self.M_o < 20.0:
                corr_ind = 0
            else:
                corr_ind = 1

            Q0  *= const.black_body_q0[corr_ind]
            Q1  *= const.black_body_q1[corr_ind]
            FUV *= const.black_body_fuv[corr_ind]

    
        self.properties['Q0']    = Q0 * self.surface_area()
        self.properties['E0']    = E0
        self.properties['Q1']    = Q1 * self.surface_area()
        self.properties['E1']    = E1
        self.properties['L_FUV'] = FUV * self.surface_area()

        self.Mdot_ej = 0.0

        #
        # Interpolate and store wind and supernova abundances
        #

        yields = self.compute_stellar_wind_yields()
        
        i = 0 
        for e in self.wind_ejecta_abundances.keys():
            self.wind_ejecta_abundances[e] = yields[i]
            i = i + 1

        # convert to abundances
        self.properties['M_wind_total'] = self.wind_ejecta_abundances['m_tot']
        if self.properties['M_wind_total'] > 0.0:
            for e in self.wind_ejecta_abundances.keys():
                self.wind_ejecta_abundances[e] /= self.properties['M_wind_total']

        self.properties['Mdot_wind'] = 0.0
        self.properties['v_wind']    = 0.0

class StarList:

    def __init__(self, maximum_stars = None, conserve_memory = False):
                          #
        if maximum_stars != None and conserve_memory:
            self.stars = [None] * maximum_stars
        else:
            self.stars          = []

        if maximum_stars == None:
            maximum_stars = 1.0E99

        self._maximum_stars = maximum_stars

        self._conserve_memory = conserve_memory


    def evolve(self, t, dt):
        map( lambda x : x.evolve(t, dt), self.stars)
        return

    def append(self, new_star):
        # apparently a possible bug in appending objects to list with gc
        gc.disable()
        self.stars.append(new_star)
        gc.enable()
        return

    def property_asarray(self, name, star_type = 'all'):

        if len(self.stars) == 0:
            return 0.0
        
        if not star_type == 'all': # only look over a subselection of stars
            _star_subset = [x for x in self.stars if x.properties['type']== star_type]
        else:
            _star_subset = self.stars        

        if len(_star_subset) == 0:
            return np.asarray(0.0)


        if name == 'mass' or name == 'Mass' or name == 'M':
            array = np.asarray( [x.M for x in _star_subset])
        elif name == 'Z' or name == 'metallicity' or name == 'Metallicity':
            array = np.asarray( [x.Z for x in _star_subset])
        elif name == 'Mdot_ej':
            array = np.asarray( [x.Mdot_ej for x in _star_subset] )
        elif name == 'mechanical_luminosity':
            array = np.asarray( [x.mechanical_luminosity for x in _star_subset])
        elif name == 'id':
            array = np.asarray( [x.id for x in _star_subset])
        else:
            try:
                array = np.asarray( [x.properties[name] for x in _star_subset] )
            except KeyError:
                print name, "star property or value not understood for " + star_type + " stars"
                raise KeyError
        
        return array
        
    def property_names(self, mode='unique', star_type = 'all'):
        if len(self.stars) == 0:
            return None

        if not star_type == 'all':
            _star_subset = [x for x in self.stars if x.properties['type'] == star_type]
        else:
            _star_subset = self.stars

        all_keys = [x.properties.keys() for x in _star_subset]
        unique_keys = np.unique([item for sublist in all_keys for item in sublist])

        if mode == 'unique':   # get all unique keys from all stars in set
            return unique_keys 
        elif mode == 'shared': # get only properties shared among all stars in set
            i = 0
            return np.unique([element for element in unique_keys if element in all_keys[i] for i in np.arange(0,len(all_keys))])


    def species_asarray(self, name, star_type = 'all'):
        """
        Return either the ejecta rates or chemical tags for stars as array
        """

        if not star_type == 'all':
            _star_subset = [x for x in self.stars if x.properties['type'] == star_type]
        else:
            _star_subset = self.stars

        if len(_star_subset) == 0:
            return np.zeros(len(name))

        if 'Mdot' in name:
            if 'ej' in name:
                name = name.replace('Mdot_ej_','')
            else:
                name = name.replace('Mdot_','')

            func = lambda x , y : x.wind_ejecta_abundances[y]

        elif 'SN' in name:
            if 'ej' in name:
                name = name.replace('SN_ej_','')
            else:
                name = name.replace('SN_','')

            func = lambda x, y : x.sn_ejecta_abundances[y]

        else:
            func = lambda x, y : x.abundances[y]

        return np.asarray([ func(x, name) for x in self.stars])


    def Z(self):
        """
        List comprehension to return all metallicities as numpy array
        """
        return np.asarray([x.Z for x in self.stars])

    def M(self):
        """
        List comprehension to return all masses as np array
        """
        return np.asarray([x.M for x in self.stars])
    


