class StarParticle:

    def __init__(self, M, Z, tform=0.0, id = 0):

        self.M   = M
        self.M_o = M
        self.Z   = Z

        self.tform = tform
        self.id    = id

        self._assign_properties()

        self.particle_properties = {}
        self.particle_properties['status'] = ['live']

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

        if (age + dt > self.lifetime):
            self.particle_properties['status'] = ['dead']
            
            #
            # conditional statements to turn particle into WD
            # or remnant
            #


    def _assign_properties(self):

        p_list = ['luminosity', 'radius',
                  'lifetime'  , 'l_fuv',
                  'l_q1', 'l_q2', 'e_q1', 'e_q2']

 #       self._stellar_evolution_grid_position()
  #      self._radiation_properites_grid_position()

    
    
        self.L = 1.0
        self.R = 1.0
        self.lifetime = 2.5
        self.Mdot_ej = 0.0

    #    self.L        = self.compute_property("luminosity")
   #     self.R        = self.compute_property("radius")
     #   self.lifetime = self.compute_property("lifetime")



    def _stellar_evolution_grid_position()

#def WhiteDwarf(StarParticle):


    


