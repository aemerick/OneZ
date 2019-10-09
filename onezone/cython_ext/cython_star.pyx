
__author__ = "aemerick <emerick@astro.columbia.edu>"

# need to allow dimension switch in interpolation routines
# when one of the dims is exactly equal to one of the grid points
_interpolation_hack = 0.999999999 # do this for now

# --- external ---
import numpy as np
import gc
from collections import OrderedDict
import itertools


# --- internal ---
from onezone import data_tables as DT
from onezone import radiation   as rad
from onezone import physics     as phys
from onezone import config      as config


from onezone.constants import CONST as const

#
# ------- load the global data tables -----------
#
SE_TABLE  = DT.StellarEvolutionData()
RAD_TABLE = DT.RadiationData()

SN_YIELD_TABLE           = DT.StellarYieldsTable('SNII')
WIND_YIELD_TABLE         = DT.StellarYieldsTable('wind')
MASSIVE_STAR_YIELD_TABLE = DT.StellarYieldsTable('massive_star')
POPIII_YIELD_TABLE       = DT.StellarYieldsTable("popIII")

cdef class StarParticle:

    cdef double M, Z, age, t_now, M_o, tform
    cdef int id
    cdef dict properties, sn_ejecta_masses, wind_ejecta_abundances

    def __init__(self, M = None, Z = None, abundances={'m_tot':1.0}, tform=0.0, id = 0, M_o = None,
                      age = None, t_now = None):
        """
        Initialize star particle with mass and metallicity. Particle
        properties are assigned using input M and Z to interpolate
        over tables.

        Args:
            M (float): Star mass in solar masses
            Z (float): star metallicity in solar masses
            abundances (optional, dict): A dictionary of local chemical
                abundances in the SF region. If this is passed,
                star is assigned enrichment properties. Default to None.
            tform (optional, float): Time at formation. Default is 0.0.
            id (optional int): Optional unique id number to assign to particle
                if tracking many. Default is 0

        """
        if M is None or Z is None:
            raise ValueError("Must set values for mass and metallicity")

        self.M   = M
        if M_o is None:
            self.M_o = M
        else:
            self.M_o = M_o
        self.Z   = Z

        self.tform = tform
        self.id    = id

        if age is None and t_now is None:
            self.age = 0.0

        elif age is None:
            self.age = t_now - self.tform

        else:
            self.age = age

            if (t_now - self.age) != self.tform:
                _my_print("Supplied Particle age and formation time do not agree with current time")

            self.tform = t_now - self.age

        self.properties = {}

        self.wind_ejecta_abundances = {} #OrderedDict()
        self.sn_ejecta_masses       = {} #OrderedDict()

        if not abundances is None:
            for e in abundances.keys():
                self.wind_ejecta_abundances[e] = 0.0
                self.sn_ejecta_masses[e]  = 0.0


        return

    def evolve(self, t, dt, ej_masses = {}, sn_masses = {},
                            snII_counter = -9999, snIa_counter = -9999,
                            special_accumulator={}):
        pass

    def _assign_properties(self):
        """
        Compute stellar properties given M and Z
        """
        pass

#
# Star_types: 'star', 'new_WD'
#

cdef class Star(StarParticle):

    cdef double Mdot_ej

    def __init__(self, star_type = 'star', *args, **kwargs):

        StarParticle.__init__(self, *args, **kwargs)

        self.properties['type'] = star_type

        self._assign_properties()

        if 'abundances' in kwargs:
            self.write_abundance(kwargs['abundances'])

        return


    def evolve(self, double t, double dt, ej_masses = {}, sn_masses = {},
                            int snII_counter = -9999, int snIa_counter = -9999,
                            special_accumulator={}):
        """
        Evolve
        """

        self.age = t - self.tform

        #
        # check and update Mdot_ej from stellar winds
        #
        self.Mdot_ej = 0.0
        self.stellar_wind_parameters(self.age, dt)
        self.Mdot_ej = self.properties['Mdot_wind']

        cdef double SN_mass_loss = 0.0

        if (self.age + dt > self.properties['lifetime'] / (config.units.time)):

            if 'new' in self.properties['type']:
                #
                # star changed types in previous timestep, update to "old"
                # and do nothing else
                #
                self.properties['type'] = self.properties['type'].replace('new_','')
                self._clear_properties()

                SN_mass_loss = 0.0

            elif self.properties['type'] == 'star':
                if self.M_o > config.stars.SNII_mass_threshold  and\
                   self.M_o < config.stars.direct_collapse_mass_threshold:
                    #
                    # Core collapse supernova - change type and compute yields
                    #
                    self.set_SNII_properties()
                    self.properties['type'] = 'new_remnant'
                    SN_mass_loss = self.sn_ejecta_masses['m_tot']
                    snII_counter += 1
                elif self.M_o < config.stars.SNII_mass_threshold:
                    #
                    # Otherwise, form a white dwarf when dead and label as
                    # candidate for future SNIa
                    #
                    self.properties['type']               = 'new_WD'

                    if self.M_o > config.stars.SNIa_candidate_mass_bounds[0]\
                           and self.M_o < config.stars.SNIa_candidate_mass_bounds[1]:

                        self.properties['SNIa_candidate'] = True
                        self.properties['WD_lifetime']    = phys.WD_lifetime(t,
                                                                             self.tform,
                                                                             self.properties['lifetime']/config.units.time,
                                                                             config.stars.DTD_slope,
                                                                             config.stars.NSNIa,
                                                                             config.zone.current_redshift) * config.units.time

                    else:
                        self.properties['SNIa_candidate'] = False

                else:
                    #
                    # direct collapse to black hole - no supernova
                    #
                    self.properties['type'] = 'new_direct_collapse'

            # if this is a WD, need to check and see if it will explode
            if self.properties['type'] == 'WD':

                if self.properties['SNIa_candidate']:

                    if self.properties['WD_lifetime'] + self.tform <= t*config.units.time:

                        # go Type Ia supernova
                        self.properties['type'] = 'new_SNIa_remnant'
                        self.set_SNIa_properties()
                        SN_mass_loss = self.sn_ejecta_masses['m_tot']
                        snIa_counter += 1

        #
        # Compute total mass lost through supernova and wind
        #
        cdef double M_loss = self.Mdot_ej * (config.units.time) * dt + SN_mass_loss

        self.M = self.M - M_loss

        if self.M < 0.0 and not 'SNIa' in self.properties['type']:
            _my_print("ERROR IN STAR: Negative stellar mass in particle type " + self.properties['type'])
            _my_print("birth mass, mass, mdot_ej, mdot_ej*dt, sn_mass_loss, M_loss, self.age")
            _my_print("%3.3E %3.3E %3.3E %3.3E %3.3E %3.3E %3.3E"%(self.M_o, self.M, self.Mdot_ej, self.Mdot_ej*dt, SN_mass_loss, M_loss, self.age))
            _my_print(self.properties)
            _my_print("time, dt")
            _my_print("%3.3E %3.3E"%(t, dt))
            raise RuntimeError
        elif self.M < 0.0 and 'SNIa' in self.properties['type']:
            self.M = 0.0

        if self.properties['type'] == 'new_WD':
            self.M = phys.white_dwarf_mass(self.M_o)

        #
        # add in ejected mass for
        #   1) winds
        #   2) SN explosion
        #
        if self.properties['type'] == 'star' or\
           self.properties['type'] == 'new_WD':

            for key in ej_masses.keys():
                ej_masses[key] += self.wind_ejecta_abundances[key] * self.Mdot_ej

            if self.M_o > config.zone.track_massive_star_ejecta_mass:
                special_accumulator['m_massive'] += self.wind_ejecta_abundances['m_metal'] * self.Mdot_ej

        elif self.properties['type'] == 'new_remnant':
            # sn may have both winds and SN ejecta if explosion
            # happens between timesteps (almost always)
            for key in ej_masses.keys():
                ej_masses[key] += self.wind_ejecta_abundances[key] * self.Mdot_ej
                sn_masses[key] += self.sn_ejecta_masses[key]

            if self.M_o > config.zone.track_massive_star_ejecta_mass:
                special_accumulator['m_massive'] += self.wind_ejecta_abundances['m_metal']*self.Mdot_ej +\
                                                    self.sn_ejecta_masses['m_metal']

        elif self.properties['type'] == 'new_SNIa_remnant':

            for key in sn_masses.keys():
                sn_masses[key] += self.sn_ejecta_masses[key]

            if self.M_o > config.zone.track_massive_star_ejecta_mass:
                special_accumulator['m_massive'] += self.sn_ejecta_masses['m_metal']


        return None

    def set_SNIa_properties(self, check_mass = False):
        """
        If a SNIa candidate, sets SNIa ejecta masses in array. check_mass
        call is in place for versatility. One can set SNIa properties BEFORE turning
        particle into remnant using check_mass = False (default), or if using this
        code for post-processing, set check_mass = True to ensure that SNIa ejecta
        masses are non-zero ONLY if a candidate has actually exploded.
        """
        # need to rename and combine this and functions below

        if not config.stars.use_snIa:
            return

        if ((self.M_o > config.stars.SNIa_candidate_mass_bounds[0] and\
           self.M_o < config.stars.SNIa_candidate_mass_bounds[1] and (not check_mass)) or\
          (check_mass and (self.M_o > config.stars.SNIa_candidate_mass_bounds[0] and\
           self.M_o < config.stars.SNIa_candidate_mass_bounds[1] and self.M == 0.0))) :

            if len(self.wind_ejecta_abundances.keys()) > 0:
                yields = phys.SNIa_yields(self.wind_ejecta_abundances.keys())

                i = 0
                for e in self.sn_ejecta_masses.keys():
                    self.sn_ejecta_masses[e] = yields[i]
                    i = i + 1
                #print "------------------------", self.sn_ejecta_masses

        else:
            return NotImplementedError

        return yields

    def set_SNII_properties(self):

        if not config.stars.use_snII:
            return

        if len(self.wind_ejecta_abundances.keys()) > 0:

            if self.M_o < config.stars.direct_collapse_mass_threshold and\
               self.M_o > config.stars.SNII_mass_threshold :

                if self.M_o < config.data.yields_mass_limits[1]:
                    yields =  SN_YIELD_TABLE.interpolate([self.M_o, self.Z],
                                                          self.wind_ejecta_abundances.keys())
                elif config.stars.extrapolate_snII_yields:
                    yields = np.asarray(SN_YIELD_TABLE.interpolate([config.data.yields_mass_limits[1] * _interpolation_hack, self.Z],
                                                          self.wind_ejecta_abundances.keys()))
                    yields = yields * self.M_o / (config.data.yields_mass_limits[1] * _interpolation_hack)

            else:
                # direct collapse supernova - no SN mass injection
                yields = np.zeros(len(self.sn_ejecta_masses.keys()))

            i = 0
            for e in self.sn_ejecta_masses.keys():
                self.sn_ejecta_masses[e] = yields[i]
                i = i + 1

        else:
            raise NotImplementedError



#    @property
    cdef double mechanical_luminosity(self):
        return 0.5 * self.properties['Mdot_wind'] * const.Msun * self.properties['v_wind']**2

#    @property
    cdef double total_wind_thermal_energy(self):
        return 1.5 * const.k_boltz * (self.wind_ejecta_masses()['m_tot']*const.Msun)*\
                 self.properties['Teff'] / (0.65*const.m_p)

#    @property
    cdef double total_wind_kinetic_energy(self):
        return 0.5 * self.wind_ejecta_masses()['m_tot']*const.Msun * self.properties['v_wind']**2


    cdef double surface_gravity(self):
        return const.G * self.M * const.Msun / self.properties['R']**2

    cdef double surface_area(self):
        return 4.0 * np.pi * self.properties['R']**2

    cdef void _clear_properties(self):
        """
        zeroes certain properties after star dies
        """
        zero_properties = ['E0', 'E1', 'L_FUV', 'L_LW', 'Q0', 'Q1',
                           'luminosity', 'v_wind', 'Mdot_wind']

        self.Mdot_ej = 0.0
        for p in zero_properties:
            self.properties[p] = 0.0

        return

    cdef double ionizing_photons(self, photon_type):

        if not 'q' in photon_type:
            if photon_type == 'HI':
                photon_type = 'q0'
            elif photon_type == 'HeI':
                photon_type = 'q1'

        if self.properties.haskey(photon_type):
            return self.properties[photon_type] * self.surface_area()
        else:
            return 0.0

    cdef double fuv_luminosity(self):
        if self.properties.haskey('FUV_flux'):
            return self.properties['FUV_flux'] * self.surface_area()
        else:
            return 0.0

    cdef double LW_luminosity(self):
        if self.properties.haskey('LW_flux'):
            return self.properties['LW_flux'] * self.surface_area()
        else:
            return 0.0

    cdef void stellar_wind_parameters(self, age, dt):

        if not self.properties['type'] == 'star' or not config.stars.use_stellar_winds:
            return

        #if (self.M_o <= config.data.yields_mass_limits[1]) and\
        #   (self.M_o >= config.data.yields_mass_limits[0]):
        if True:
            # need to compute wind velocities for all stars
            # check if star's wind is ON
            do_wind = True

            if (self.M_o < config.stars.AGB_wind_phase_mass_threshold) and config.stars.use_AGB_wind_phase:
                if self.age + dt < self.properties['age_agb'] / config.units.time:
                    do_wind = False
                    wind_lifetime = 0.0
                else:
                    wind_lifetime = (self.properties['lifetime'] - self.properties['age_agb'])

            else: # else have star wind on for entire lifetime
                wind_lifetime = self.properties['lifetime']


            if wind_lifetime < dt * config.units.time:
                wind_lifetime = dt * config.units.time


            if do_wind and self.age * config.units.time < self.properties['lifetime']:
                Mdot   = self.properties['M_wind_total'] / wind_lifetime
            else:
                Mdot   = 0.0

            #
            # if difference b/t birth mass and mass after wind timestep is more than
            # total amount, set wind ejected to just the difference
            # this can happen when wind phase is < dt and lines up between timesteps
            #
            final_mass = self.M - Mdot * dt * config.units.time
            correct_final_mass = self.M_o - self.properties['M_wind_total']


            if final_mass < correct_final_mass:
                Mdot = (self.M - correct_final_mass) / wind_lifetime

        else:

            Mdot  = phys.s99_wind_mdot(self.properties['luminosity'], self.M_o,
                                      self.properties['Teff'], self.Z)


        if (self.M > config.stars.AGB_wind_phase_mass_threshold):
            vwind = phys.s99_wind_velocity( self.properties['luminosity'], self.M_o,
                                            self.properties['Teff'], self.Z)
        else:
            vwind = config.stars.AGB_wind_velocity * 1.0E5 # km/s -> cm/s


        self.properties['Mdot_wind'] = Mdot
        self.properties['v_wind']    = vwind
        return

    def compute_stellar_wind_yields(self):
        """ compute_stellar_wind_yields

        Computes yields from stellar winds for all considered species using
        either the NuGrid table (1 < M < 25) or the PARSEC table ( M > 25). If
        the parsec table is NOT used (config.stars.use_massive_star_yields == FALSE),
        then extrapolation from the NuGrid table is used (which is very wrong).

        Returns:
            numpy array of total yields over lifetime for each element, sorted in
            atomic number order
        """


        if( self.M_o < config.data.yields_mass_limits[1] ):

            yields = np.asarray(WIND_YIELD_TABLE.interpolate([self.M_o, self.Z],
                                                              self.wind_ejecta_abundances.keys()))
        elif (config.stars.use_massive_star_yields):
            # use yields from PARSEC massive star yields
            yields = np.asarray(MASSIVE_STAR_YIELD_TABLE.interpolate([self.M_o, self.Z],
                                                                     self.wind_ejecta_abundances.keys()))

        else:
            #
            # For stars off of the grid, scale most massive star
            # to current mass.
            #
            yields = np.asarray(WIND_YIELD_TABLE.interpolate([config.data.yields_mass_limits[1]*_interpolation_hack, self.Z], self.wind_ejecta_abundances.keys()))
            yields = yields * self.M_o / (config.data.yields_mass_limits[1] * _interpolation_hack)


        return np.asarray(yields)


    cdef void _assign_properties(self):

        # list of properties assigned in this function (remember to update!!)
        p_list = ['luminosity', 'radius',
                  'lifetime'  , 'age_agb', 'L_FUV', 'L_LW', 'agb_phase_length',
                  'Q0', 'E0', 'Q1', 'E1', 'Mdot_wind', 'v_wind', 'M_wind_total']

        if self.properties['type'] == 'unresolved_star':
            self.Mdot_ej           = 0.0
            for p in p_list:
                self.properties[p] = 0.0
            return

        L, T, R, lifetime, age_agb = SE_TABLE.interpolate([self.M_o,self.Z], ['L','Teff','R','lifetime','age_agb'])
        self.properties['luminosity']  = L * const.Lsun
        self.properties['Teff']        = T
        self.properties['R']           = R
        self.properties['lifetime']    = lifetime
        self.properties['age_agb']     = age_agb
        self.properties['agb_phase_length']  = lifetime - age_agb


        Q0, Q1, FUV, LW = RAD_TABLE.interpolate([self.properties['Teff'],
                                             self.surface_gravity(),
                                             self.Z], ['q0','q1','FUV_flux', 'LW_flux'])

        E0  = rad.average_energy(const.E_HI/ const.eV_erg, self.properties['Teff'])
        E1  = rad.average_energy(const.E_HeI/const.eV_erg, self.properties['Teff'])


        use_blackbody = False
        for a in [Q0, Q1, FUV, LW]:
            if a == 'offgrid':
                use_blackbody = True
                break;

        if use_blackbody:
            FUV = rad.fuv_flux_blackbody(self.properties['Teff'])
            LW  = rad.LW_flux_blackbody(self.properties['Teff'])
            Q0  = rad.compute_blackbody_q0(self.properties['Teff'])
            Q1  = rad.compute_blackbody_q1(self.properties['Teff'])

            if config.stars.normalize_black_body_to_OSTAR:
                if self.M_o < config.stars.black_body_correction_mass:
                    corr_ind = 0
                else:
                    corr_ind = 1

                Q0  *= config.stars.black_body_q0_factors[corr_ind]
                Q1  *= config.stars.black_body_q1_factors[corr_ind]
                FUV *= config.stars.black_body_FUV_factors[corr_ind]
                LW  *= config.stars.black_body_LW_factors[corr_ind]


        self.properties['Q0']    = Q0 * self.surface_area()
        self.properties['E0']    = E0
        self.properties['Q1']    = Q1 * self.surface_area()
        self.properties['E1']    = E1
        self.properties['L_FUV'] = FUV * self.surface_area()
        self.properties['L_LW']  = LW  * self.surface_area()

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

        self.set_SNII_properties()

        self.properties['Mdot_wind'] = 0.0
        self.properties['v_wind']    = 0.0
        return

    cdef dict wind_ejecta_masses(self):

        mass = {} # OrderedDict()

        for k in self.wind_ejecta_abundances.keys():
            mass[k] = self.wind_ejecta_abundances[k] * self.properties['M_wind_total']

        return mass

    cdef dict return_sn_ejecta_masses(self):

        self.set_SNII_properties()

        return self.sn_ejecta_masses

    cdef void write_abundance(self, abundances):
        """
        Write star particle abundances to file. In principle, if homogenous one-zone
        model is used this is not necessary as abundances will just be identical to gas
        in zone and can be obtained by matching star particle formation time to
        the abundance outputs in the summary files. However, this is useful when
        summary file is not written to every time step OR in the future when
        multi-zones are in use.

        In actuality this just adds the star to a buffer to limit IO
        """

        if not (config.io._abundance_buffer is None):

            config.io._abundance_buffer.flush() # flush if needed

            i = config.io._abundance_buffer.count
            config.io._abundance_buffer.buffer[i][0] = self.tform
            config.io._abundance_buffer.buffer[i][1] = self.id
            config.io._abundance_buffer.buffer[i][2] = self.M
            config.io._abundance_buffer.buffer[i][3] = self.Z
            config.io._abundance_buffer.buffer[i][4] = self.properties['lifetime']

            for ei,e in enumerate(config.zone.species_to_track):
                config.io._abundance_buffer.buffer[i][ei+5] = abundances[e]

            config.io._abundance_buffer.count += 1

        return

class StarList:
    """
    List of star objects with useful functions to handle operating
    on many stars at once.
    """

    def __init__(self, stars = None):


        if stars is None:
            if config.zone.maximum_stars != None and config.zone.optimize:
                self._stars           = [None] * config.zone.maximum_stars
                self._stars_optimized = True
            else:
                self._stars           = []
                self._stars_optimized = False


            self._are_there_new_stars = False
            self._N_stars             = 0
        else:
            self._stars = stars
            self._N_stars = len(self._stars)
            self._are_there_new_stars = True

        return

    def evolve(self, t, dt, *args, **kwargs):

        #map( lambda x : x.evolve(t, dt, *args, **kwargs), self.stars_iterable)
#        Was debating trying to multi-thread here but that might not do anything
#        if config.multiprocess:
#
#        else:
        for x in self.stars_iterable:
            x.evolve(t,dt,*args,**kwargs)

        return

    @property
    def stars(self):
        if self._stars_optimized:
            return self._stars[:self.N_stars]
        else:
            return self._stars

    @property
    def stars_iterable(self):
        return itertools.islice(self._stars, None, self.N_stars)

    @property
    def N_stars(self):

        #if self._are_there_new_stars:

        #    if self._stars_optimized :
        #        self._N_stars = len(self._stars) - self._stars.count(None)
        #    else:
        #        self._N_stars = len(self._stars)

        return self._N_stars


    def _values_outdated(self):

        return self._internal_time < config.global_values.time

    def add_new_star(self, new_star):
        """
        Adds a new star to the list. Either appends the star to the list
        or fills it in the list of memory optimization is on and a maximum
        number of stars is provided
        """

        config.global_values.profiler.start_timer('add_new_star', True)

        if self._stars_optimized:
            #
            # add to last element in list
            #
            self._stars[ self.N_stars ] = new_star
        else:

            self._append(new_star)

        self._N_stars += 1

        config.global_values.profiler.end_timer('add_new_star')
        return



    def _append(self, new_star):
        # apparently a possible bug in appending objects to list with gc
        gc.disable()
        self._stars.append(new_star)
        gc.enable()
        return

    def property_asarray(self, name, star_type = 'all', subset_condition = None):

        config.global_values.profiler.start_timer('property_asarray', True)

        if self.N_stars == 0:
            config.global_values.profiler.end_timer('property_asarray')
            return np.zeros(1)

        if not star_type == 'all':
            _star_subset = self.get_subset( lambda x : x.properties['type'] != star_type )
        else:
            _star_subset = self.stars_iterable

        if not subset_condition == None:
            for key in subset_condition.keys():
                _star_subset = self._get_subset( _star_subset, subset_condition[key])

        if name == 'mass' or name == 'Mass' or name == 'M':
            array = np.asarray( [x.M for x in _star_subset])
        elif name == 'initial_mass' or name == 'M_o' or name == 'birth_mass':
            array = np.asarray( [x.M_o for x in _star_subset])
        elif name == 'Z' or name == 'metallicity' or name == 'Metallicity':
            array = np.asarray( [x.Z for x in _star_subset])
        elif name == 'Mdot_ej':
            array = np.asarray( [x.Mdot_ej for x in _star_subset] )
        elif name == 'mechanical_luminosity':
            array = np.asarray( [x.mechanical_luminosity for x in _star_subset])
        elif name == 'id':
            array = np.asarray( [x.id for x in _star_subset])
        elif name == 'age':
            array = np.asarray( [x.age for x in _star_subset])
        else:
            try:
                array = np.asarray( [x.properties[name] for x in _star_subset] )
            except KeyError:
                _my_print( name + " star property or value not understood for " + star_type + " stars")
                raise KeyError

        #
        # as can happen if there are no stars in subset
        #
        config.global_values.profiler.end_timer('property_asarray')

        if len(array) == 0:
            if name == 'type':
                return [None]
            else:
                return np.zeros(1)
        else:
            return array




    def property_names(self, mode='unique', star_type = 'all'):
        if self.N_stars == 0:
            return None

        if not star_type == 'all':
            _star_subset = self.get_subset( lambda x : x.properties['type'] != star_type )
        else:
            _star_subset = self.stars_iterable


        all_keys = [x.properties.keys() for x in _star_subset]
        unique_keys = np.unique([item for sublist in all_keys for item in sublist])

        if mode == 'unique':   # get all unique keys from all stars in set
            return unique_keys
        elif mode == 'shared': # get only properties shared among all stars in set
            i = 0
            return np.unique([element for element in unique_keys if element in all_keys[i] for i in np.arange(0,len(all_keys))])

    def get_subset(self, expr):
        """
        Get subset of stars that have a FALSE value for the desired expression.
        Previously did list method, but now attempting to handle this with an iterable
        to (hopefully) substantially improve things. For example, to get stars
        that are of type star:

            >> obj.get_subset( lambda x : x.properties['type'] != 'star' )

        Or to get stars between (and including) 10 and 20 solar masses ( M over [10,20]):

            >> obj.get_subset( lambda x : (x.M < 10.0) + (x.M > 20.0) )

        Returns iterable
        """

        config.global_values.profiler.start_timer('get_subset',True)
        res = [ x for x in self.stars_iterable if not expr(x) ]
        config.global_values.profiler.end_timer('get_subset')
        return res

#        return itertools.ifilterfalse( expr,  self.stars)

    def _get_subset(self, iterable, expr):
        return [x for x in iterable if expr(x)]

    def species_asarray(self, name, star_type = 'all'):
        """
        Return either the ejecta rates or chemical tags for stars as array
        """

        if not star_type == 'all':
            _star_subset = self.get_subset( lambda x : x.properties['type'] != star_type )
        else:
            _star_subset = self.stars_iterable

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

            func = lambda x, y : x.sn_ejecta_masses[y]

        else:
            func = lambda x, y : x.abundances[y]

        return_list = np.asarray([ func(x, name) for x in _star_subset])

        if len(return_list) == 0:
            return np.zeros(1)
        else:
            return return_list

    def Z(self):
        """
        List comprehension to return all metallicities as numpy array
        """
        return np.asarray([x.Z for x in self.stars_iterable])

    def M(self):
        """
        List comprehension to return all masses as np array
        """
        return np.asarray([x.M for x in self.stars_iterable])

    def M_o(self):
        """
        Return all initial masses of stars as np array
        """

        return np.asarray([x.M_o for x in self.stars_iterable])


def _my_print(string):
    print('[Star]: ' + string)
