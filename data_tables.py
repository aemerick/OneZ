"""
Will contain functions and classes for loading in needed
data and data tables
"""

from collections import OrderedDict

#
# need to code this up as a global set in setup.py
# 
install_dir = '/home/emerick/code/onezone/'

class DataTable:

    def __init__(self, name):
        self.name = name

        # all inherited classes need these:
        self.x           = OrderedDict()
        self.y           = OrderedDict()
        self.nbins       = OrderedDict()
        self.ndim        = 0
        self.dim_names   = []
        self._array_size = 0

    def read_data(self):
        pass

    def interpolate(self, vals, yname, *args, **kwargs):
        """
        Vals is either dict or an ordered list where order
        matches the order of the dimensions in
        """

        if isinstance(vals, dict):  
            # make sure all keys exist
            not_exist = [not vals.has_key(x) for x in self.dim_names]
            if any(not_exist):
                print "Need to supply all values ", self.dim_names       
                print "only gave", self.vals.keys()
                raise KeyError

            vals_list = [ vals[x] for x in self.dim_names]
        else:
            vals_list = vals


        return self._interpolate(vals_list, self.x.values(), self.y[yname],
                                            *args, **kwargs)


    @classmethod
    def _interpolate(cls, vals, val_arrays, y, silence = False,
                           flag = "offgrid", special_errval = None, special_flag = 'errval'):
        """
        Interpolate along an n+1 axis (y) that is a function of
        n variables, where n should be the length of vals,
        the coordinates to interpolate, and val_arrays is a
        list of n axis points to interpolate along. Returns
        the interpolated value. Some of the underlying machinery
        is generalized to arbitrary interpolation degree, but
        not all... supports up to trillinear interpolation (n <=3).
        """

        # obtain interpoaltion coefficients and index for nearest
        # grid points
        c, id = _interpolation_coefficients(vals, val_arrays, 
                                                  silence, flag, special_errval)

        n       = len(vals) # number of dimensions

        if len(c) == 1 and len(id) == 1 and n != 1:
            # likely got a bad point, check and return flags
            if silence and c == flag and id == flag:
                return flag
            else: # something is wrong if this happens
                raise RuntimeError 
        else:
            print "interpolation coefficients don't match dimensions - something broke"
            raise RuntimeError # something is wrong
    
        # check if user supplied special errval where grid may exist but 
        # the values are erroneous... Let user know with their provided flag       
        if special_errval != None: 
            if self._check_for_special_errval(n, id, y, special_errval):
                return special_flag

        # otherwise, actually perform the interpolation

        if n == 1: # linear interpolation

            yval = (1.0 - c[0]) * y[id[0]] + (c[0]) * y[id[0]]

        elif n == 2: # bilinear interpolation
            i,j = id

            yval = (1.0 - c[0]) * (1.0 - c[1]) * y[i  ][j  ] +\
                   (1.0 - c[0]) * (      c[1]) * y[i  ][j+1] +\
                   (      c[0]) * (      c[1]) * y[i+1][j+1] +\
                   (      c[0]) * (1.0 - c[1]) * y[i+1][j  ]

        elif n == 3: # trilinear interpolation
            i,j,k = id

            yval = (1.0 - c[0])*(1.0 - c[1])*(1.0 - c[2]) * y[i  ][j  ][k  ] +\
                   (1.0 - c[0])*(      c[1])*(1.0 - c[2]) * y[i  ][j+1][k  ] +\
                   (      c[0])*(      c[1])*(1.0 - c[2]) * y[i+1][j+1][k  ] +\
                   (      c[0])*(1.0 - c[1])*(1.0 - c[2]) * y[i+1][j  ][k  ] +\
                   (1.0 - c[0])*(1.0 - c[1])*(      c[2]) * y[i  ][j  ][k+1] +\
                   (1.0 - c[0])*(      c[1])*(      c[2]) * y[i  ][j+1][k+1] +\
                   (      c[0])*(      c[1])*(      c[2]) * y[i+1][j+1][k+1] +\
                   (      c[0])*(1.0 - c[1])*(      c[2]) * y[i+1][j  ][k+1]

        elif n > 3: # does not support degree above 3
            print "We do not support n > 3 dimensional interpolation"
            raise RuntimeError

        return yval

    @classmethod
    def _check_for_special_errval(cls, n, id, y, errval):
        """
        This is meant to cat "special" grid values that "exist" on the grid but
        may have erroneous interpalted values. Catch these and return to user
        this only occurs if user desires
        """


        if n == 1:
            if y[id[0]] == special_errval or y[id[0]+1] == special_errval:
                return True

        elif n == 2:

            if y[i  ][j  ] == special_errval or y[i  ][j+1] == special_errval or\
               y[i+1][j+1] == special_errval or y[i+1][j  ] == special_errval:
                 
                return True

        elif n == 3:
            for k in [0,1]:
                for j in [0,1]:
                    for i in [0,1]:
                        if y[id[0] + i][id[1] + j][id[2] + k] == special_errval:
                            return True

    @classmethod
    def _interpolation_coefficients(cls, vals, val_arrays,
                                         silence=False, flag = "offgrid", special_errval=None):
        """
        vals is list of parameters.
        val_arrays is list of arrays to interpolate over
        """

        n = len(vals)

        coeff   = np.zeros(n)
        indeces = np.zeros(n)

        # perform a linear interpolation in each dimension
        for i in np.arange(n):
            coeff[i], indeces[i] = \
               self._linear_interpolation_coefficients(vals[i], val_arrays[i])

        if silence:
            if any(coeff) == flag:
                return flag, flag

        return coeff, indeces

    @classmethod
    def _linear_interpolation_coefficients(cls, x, xarray, 
                                           silence=False, flag="offgrid", special_errval=None):

        # check bounds in this dimension
        if (x < xarray[0] or x > xarray[-1]):

            # In some cases data off grid is O.K. In this case don't complain,
            # instead return the flag values so calling function can do somthing
            # else for data off grid
            if silence:     
                return flag, flag

            print x, xarray[0], xarray[-1]
            raise ValueError

        i = np.abs( x - xarray).argmin()
        if ( x < xarray[i]):
          i = i - 1

        t = (x - xarray[i]) / (xarray[i+1] - xarray[i])

        return t, i        

    @property
    def y_names(self):
        return self.y.keys()

    @property
    def y_values(self):
        return self.y.values()




# in rad data class, do 
# try:
#      interpolation coeff
# except ValueError:
#      off grid, use black body

class StellarEvolutionData(DataTable):

    def __init__(self, manual_table = False):
        DataTable.__init__("Stellar Evolution Data Table")

        self.ndim = 2
        self.dim_names = ['mass','metallicity']

        if not manual_table:
            self.read_data()

    def read_data(self, data_dir = None):

        if data_dir = None:
            self.data_dir = install_dir + 'Data/'

        # 
        self.x['mass'] = np.array( [  1.0,  2.0,  3.0,  4.0,   5.0,
                                      6.0,  7.0,  8.0, 10.0,  12.0,
                                     14.0, 16.0, 18.0, 20.0,  24.0,
                                     28.0, 30.0, 35.0, 40.0,  50.0,
                                     55.0, 60.0, 70.0, 90.0, 100.0, 
                                    120.0])

        self.x['metallicity'] = np.array( [0.0001, 0.0002, 0.0005, 0.001, 
                                           0.002, 0.004, 0.006, 0.008,
                                           0.01, 0.014, 0.017] )

        # hard code this for now but need to generalize
        self.nbins['mass'] = np.size(self.x['mass'])
        self.nbins['metallicity'] = np.size(self.x['metallicity'])

        self._array_size = int( np.prod(self.nbins.values()))


        # now read in each of the data sets
        for name in ['L', 'Teff', 'R', 'lifetime', 'agb_age']:

            self.y[name] = np.zeros(self._array_size)
            self.y[name] = self.y.reshape( self.nbins.keys() )

        # now read in the data
        i = 0; j = 0
        data = np.genfromtxt(self.data_dir + 'parsec_data.in')
        for line in data:
          
            for counter in [0,1,2]: # L T and R are logged - 1st 2 cols are M, Z
                self.y[ self.y_names[counter]  ][i][j] = 10.0**(line[counter+2])

            for counter in [3,4]: # lifetime and agb are not - skip first two cols
                self.y[ self.y_names[counter]  ][i][j] = line[counter+2]

            j = j + 1
            if j >= (self.nbins.values())[1]:
                j = 0
                i = i + 1

        return None     


class RadiationData(DataTable):

    def __init__(self, manual_table = False):
        DataTable.__init__("Radiation data table")

        self.ndim      = 3
        self.dim_names = ['temperature','surface_gravity','metallicity']


        if not manual_table:
            self.read_data()

    def read_data(self, data_dir = None):
        # if data file is not provided, use internal data

        if data_dir == None:
            self.data_dir = install_dir + 'Data/'

        # hard code this for now, but need to generalize
        self.nbins['temperature']     = 12
        self.nbins['surface_gravity'] = 8
        self.nbins['metalllicity']    = 10
        self._array_size              = int( np.prod(self.nbins.values()))

        # set values for each dimension
        self.x['temperature']     = np.arange(27500.0, 57500.0, 2500.0)
        self.x['surface_gravity'] = 10.0**(np.arange(3.0,5.0,0.25))
        self.x['metallicity']     = np.array([0.0,0.001,0.01,1/50.0,1/30.0,0.1,0.2,0.5,1.0,2.0])

        # now, read in each of the data sets
        self.y['q0'] = None; self.y['q1'] = None; self.y['FUV_flux'] = None
        self._data_file_names = {'q0' : 'q0_rates.in' , 'q1' : 'q1_rates.in' , 
                                 'FUV_flux' : 'FUV_rates.in'}

        # make the data arrays and shape to the correct dimensions and size
        for yi in self.y:
            self.y[yi] = np.zeros(self._array_size)
            self.y[yi] = self.y[yi].reshape( (self.nbins[self.dim_names[0]],
                                              self.nbins[self.dim_names[1]],
                                              self.nbins[self.dim_names[2]))

        #
        # The below is very gross and not generalizeable as is. This should
        # be fixed for better user experience - 5/2016
        #

        # now load from each file - 
        # q0 and q1 files atm have reverse ordered metalliciites atm
        # fuv flux file is in value order - this needs to be changes 5/2016
        for name in self.dim_names:
            data = np.genfromtxt(self.data_dir + self._data_file_names[name],
                                 usecols=(2,3,4,5,6,7,8,9,10,11))

            i = 0 ; j = 0; k = 0
            for line in data:
                for k in np.arange(np.size(line)):

                    # need to get rid of this if statement by unifying file format
                    if name == 'FUV_flux':
                        self.y[name][i][j][k] = line[k]
                    else:
                        self.y[name][i][j][np.size(line) - k - 1] = line[k]

                j = j + 1
                if ( j >= self.nbins[self.dim_names[1]] ):
                    j = 0
                    i = i + 1

        # un - log the q values
        for name in self.dim_names:
            self.y[name] = 10.0**(self.y[name])

        # flag FUV values that are off of the grid  
        self.y['FUV_flux'][ self.y['FUV_flux'] < 0.0] = -1

        return None


    def interpolate(self, vals, yname, silence = None):

        if silence == None: # default behavior is ignore off grid values 
            silence = True  # this usualy means radiation should be computed some other way

        if yname == 'q0' or yname == 'q1':

            return DataTable.interpolate(self, vals, yname, silence = silence,
                                        flag = 'offgrid', special_flag = 'offgrid',
                                        special_errval = 0.0)
        elif yname == 'FUV_flux':

            return DataTable.interpoloate(self, vals, yname, silence = silence,
                                        flag = 'offgrid', special_flag = 'offgrid',
                                        special_errval = -1.0)


    
# # #  do the below in the star... probably shouldnt be handled by the table class
#        if not interp_value == 'offgrid':
 #           return interp_value
  ##      elif allow_black_body == False:
    #        print "Radiation point off of interpolation grid and black body fix is turned off"
     #       print "either fix data table, interpolation point, or allow black body approx"
      #      raise RuntimeError

        # if we are here then use the black body approximation

   #     if yname == 'q0' or yname == 'q1':

        


    #    elif yname == 'FUV_flux':

        

