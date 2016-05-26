# --- external ---
import numpy as np

# --- internal ---
import data_tables as DT

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


        self._assign_properties()

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

    def evolve(self, t, dt):
        """
        Evolve the star
        """

        #
        # check and update Mdot_ej from stellar winds
        #
        self.Mdot_ej = 0.0

        #
        # check age
        #
        age = t - self.tform

        if (age + dt > self.properties['lifetime']):
            self.particle_properties['status'] = ['dead']
            
            #
            # conditional statements to turn particle into WD
            # or remnant
            #

    @property
    def surface_gravity(self):
        return const.G * self.properties['M']/self.properties['R']**2

    @property
    def surface_area(self):
        return 4.0 * np.pi * self.properties['R']**2

    def _assign_properties(self):

        p_list = ['luminosity', 'radius',
                  'lifetime'  , 'l_fuv',
                  'l_q1', 'l_q2', 'e_q1', 'e_q2']

        # update SE table to handle arrays 

        self.properties['luminosity'] =\
                   SE_TABLE.interpolate([self.M, self.Z], 'L')
        self.properties['Teff'] =\
                   SE_TABLE.interpolate([self.M, self.Z], 'Teff')
        self.properties['R'] =\
                   SE_TABLE.interpolate([self.M, self.Z], 'R')
        self.properties['lifetime'] =\
                   SE_TABLE.interpolate([self.M, self.Z], 'lifetime')
        self.properties['age_agb'] =\
                   SE_table.interpolate([self.M, self.Z], 'age_agb')
                                                  
        
    
        self.Mdot_ej = 0.0



    def _stellar_evolution_grid_position()

#def WhiteDwarf(StarParticle):


    


