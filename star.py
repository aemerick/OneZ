
__author__ = "aemerick <emerick@astro.columbia.edu>"

# --- external ---
import numpy as np
import gc

# --- internal ---
import data_tables as DT
import radiation   as rad
import physics     as phys

from constants import CONST as const

# ------- load the global tables -----------
SE_TABLE  = DT.StellarEvolutionData()
RAD_TABLE = DT.RadiationData()

class StarParticle:

    def __init__(self, M, Z, tform=0.0, id = 0):

        self.M   = M
        self.M_o = M
        self.Z   = Z

        self.tform = tform
        self.id    = id


        self.properties = {} # constant properties
        self.properties['status'] = 'live'

        self.abundances        = {} # property of particle
        self.ejecta_abundances = {} # time varying

       

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

        #
        # check and update Mdot_ej from stellar winds
        #
        self.stellar_wind_parameters()
        self.Mdot_ej = self.properties['Mdot_wind']

        #
        # check age
        #
        age = t - self.tform

        if (age + dt > self.properties['lifetime']):
            self.properties['status'] = ['dead']

            
            
            #
            # conditional statements to turn particle into WD
            # or remnant
            #

    def mechanical_luminosity(self):
        return self.properties['Mdot_wind'] * const.Msun * self.properties['v_wind']**2

    def surface_gravity(self):
        return const.G * self.M * const.Msun / self.properties['R']**2

    def surface_area(self):
        return 4.0 * np.pi * self.properties['R']**2

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

    def stellar_wind_parameters(self, mode = None):

        Mdot  = phys.s99_wind_mdot(self.properties['luminosity'], self.M_o,
                                  self.properties['Teff'], self.Z)
        vwind = phys.s99_wind_velocity( self.properties['luminosity'], self.M_o,
                                        self.properties['Teff'], self.Z)

        self.properties['Mdot_wind'] = Mdot
        self.properties['v_wind']    = vwind



    def _assign_properties(self):

        p_list = ['luminosity', 'radius',
                  'lifetime'  , 'age_agb', 'L_FUV',
                  'Q1', 'Q2', 'E_Q1', 'E_Q2']

        # update SE table to handle arrays 

        

        self.properties['luminosity'],
        self.properties['Teff', self.properties['R'],
        self.properties['lifetime'], self.properties['age_agb'] =\
                SE_TABLE.interpolate([self.M,self.Z], ['L','Teff','R','lifetime','age_agb'])
        self.properties['luminosity'] *= const.Lsun

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
        gc.disabe()
        self.stars.append(new_star)
        gc.enable()
        return

    def property_asarray(self, name, star_type = 'all'):

        if not star_type == 'all': # only look over a subselection of stars
            _star_subset = [x for x in self.stars if x.startype == star_type]
        
        if name == 'mass' or name == 'Mass' or name == 'M':
            array = np.asarray( [x.M for x in _star_subset])
        elif name == 'Z' or name == 'metallicity' or name == 'Metallicity':
            array = np.asarray( [x.Z for x in _star_subset])
        elif name == 'Mdot_ej':
            array = np.asarray( [x.Mdot_ej for x in _star_subset] )
        elif name == 'id':
            array = np.asarray( [x.id for x in _star_subset])
        elif hasattr(_star_subset[0], name):
            array = np.asarray( [x.properties[name] for x in _star_subset] )
        else:
            print name, "star property or value not understood for " + star_type + " stars"
            raise KeyError
        
        return array
        

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
    


