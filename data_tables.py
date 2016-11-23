
__author__ = "aemerick <emerick@astro.columbia.edu>"

# --- external ---
from collections import OrderedDict
import numpy as np

# --- internal ---


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

    def interpolate(self, vals, ynames, *args, **kwargs):
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

        single_output = False
        if isinstance(ynames, basestring):
            ynames = [ynames]
            single_output = True

        y_list = [ self.y[yname] for yname in ynames ]

        return_list =  self._interpolate(vals_list, self.x.values(), y_list,
                                                    *args, **kwargs)
        if single_output:
            if isinstance(return_list, basestring):
                return return_list
            else:
                return return_list[0]
        else:
            return return_list

    @classmethod
    def _interpolate(cls, vals, val_arrays, y_list, silence = False,
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
        c, id = cls._interpolation_coefficients(vals, val_arrays, 
                                                  silence, flag, special_errval)

        n       = len(vals) # number of dimensions

        if silence:
            # likely got a bad point, check and return flags
            if  c == flag or id == flag:
                return flag

        elif len(c) != n:
            print "interpolation coefficients don't match dimensions - something broke"
            raise RuntimeError # something is wrong

        return_list = [None] * len(y_list)    
        count = 0

        for y in y_list:
            # check if user supplied special errval where grid may exist but 
            # the values are erroneous... Let user know with their provided flag       
            if special_errval != None: 
                if cls._check_for_special_errval(n, id, y, special_errval):
                    return_list[count] =  special_flag
                    count = count + 1
                    break 
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

            return_list[count] = yval
            count = count + 1

        return return_list

    @classmethod
    def _check_for_special_errval(cls, n, id, y, errval):
        """
        This is meant to cat "special" grid values that "exist" on the grid but
        may have erroneous interpalted values. Catch these and return to user
        this only occurs if user desires
        """


        if n == 1:
            if y[id[0]] == special_errval or y[id[0]+1] == errval:
                return True

        elif n == 2:

            if y[i  ][j  ] == special_errval or y[i  ][j+1] == errval or\
               y[i+1][j+1] == special_errval or y[i+1][j  ] == errval:
                 
                return True

        elif n == 3:
            for k in [0,1]:
                for j in [0,1]:
                    for i in [0,1]:
                        if y[id[0] + i][id[1] + j][id[2] + k] == errval:
                            return True

    @classmethod
    def _interpolation_coefficients(cls, vals, val_arrays,
                                         silence=False, flag = "offgrid", special_errval=None):
        """
        vals is list of parameters.
        val_arrays is list of arrays to interpolate over
        """

        n = len(vals)

        coeff   = [None]*n
        indeces = [None]*n

        # perform a linear interpolation in each dimension
        for i in np.arange(n):
            coeff[i], indeces[i] = \
               cls._linear_interpolation_coefficients(vals[i], val_arrays[i],
                                                      silence, flag, special_errval)

        if silence:
            for i in np.arange(n):
                if coeff[i] == flag:
                    return flag, flag

        

        return np.array(coeff), np.array(indeces)

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

            print "value", x, "off of grid with bounds",  xarray[0], xarray[-1]
            raise ValueError

        i = np.abs( x - xarray).argmin()
        if ( x < xarray[i]):
          i = i - 1

        t = (x - xarray[i]) / (xarray[i+1] - xarray[i])

        return t, i        

    def y_names(self):
        return self.y.keys()

    def y_values(self):
        return self.y.values()




# in rad data class, do 
# try:
#      interpolation coeff
# except ValueError:
#      off grid, use black body

class StellarEvolutionData(DataTable):

    def __init__(self, manual_table = False):
        DataTable.__init__(self, "Stellar Evolution Data Table")

        self.ndim = 2
        self.dim_names = ['mass','metallicity']

        if not manual_table:
            self.read_data()

    def read_data(self, data_dir = None):

        if data_dir == None:
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
        for name in ['L', 'Teff', 'R', 'lifetime', 'age_agb']:

            self.y[name] = np.zeros(self._array_size)
            
            self.y[name] = (self.y[name]).reshape( tuple(self.nbins.values()) )

        # now read in the data
        i = 0; j = 0
        data = np.genfromtxt(self.data_dir + 'parsec_data.in')
        for line in data:
          
            for counter in [0,1,2]: # L T and R are logged - 1st 2 cols are M, Z
                self.y[ self.y_names()[counter]  ][i][j] = 10.0**(line[counter+2])

            for counter in [3,4]: # lifetime and agb are not - skip first two cols
                self.y[ self.y_names()[counter]  ][i][j] = line[counter+2]

            j = j + 1
            if j >= (self.nbins.values())[1]:
                j = 0
                i = i + 1

        return None     


class RadiationData(DataTable):

    def __init__(self, manual_table = False):
        DataTable.__init__(self, "Radiation data table")

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
        self.nbins['metallicity']    = 10
        self._array_size              = int( np.prod(self.nbins.values()))

        # set values for each dimension
        self.x['temperature']     = np.arange(27500.0, 57500.0, 2500.0)
        self.x['surface_gravity'] = 10.0**(np.arange(3.0,5.0,0.25))
        self.x['metallicity']     = np.array([0.0,0.001,0.01,1/50.0,1/30.0,0.1,0.2,0.5,1.0,2.0])

        # now, read in each of the data sets
        self.y['q0'] = None; self.y['q1'] = None; self.y['FUV_flux'] = None; self.y['LW_flux'] = None
        self._data_file_names = {'q0' : 'q0_rates.in' , 'q1' : 'q1_rates.in' , 
                                 'FUV_flux' : 'FUV_rates.in', 'LW_flux' : 'LW_rates.in'}

        # make the data arrays and shape to the correct dimensions and size
        for yi in self.y:
            self.y[yi] = np.zeros(self._array_size)
            self.y[yi] = self.y[yi].reshape( (self.nbins[self.dim_names[0]],
                                              self.nbins[self.dim_names[1]],
                                              self.nbins[self.dim_names[2]] ))

        #
        # The below is very gross and not generalizeable as is. This should
        # be fixed for better user experience - 5/2016
        #

        # now load from each file - 
        # q0 and q1 files atm have reverse ordered metalliciites atm
        # fuv flux file is in value order - this needs to be changes 5/2016
        for name in self.y.iterkeys():
            data = np.genfromtxt(self.data_dir + self._data_file_names[name],
                                 usecols=(2,3,4,5,6,7,8,9,10,11))

            i = 0 ; j = 0; k = 0
            for line in data:
                for k in np.arange(np.size(line)):

                    # need to get rid of this if statement by unifying file format
                    if name == 'FUV_flux' or name == 'LW_flux':
                        self.y[name][i][j][k] = line[k]
                    else:
                        self.y[name][i][j][np.size(line) - k - 1] = line[k]

                j = j + 1
                if ( j >= self.nbins[self.dim_names[1]] ):
                    j = 0
                    i = i + 1

        # un - log the q values
        for name in ['q0','q1']:
            self.y[name] = 10.0**(self.y[name])

        # flag FUV values that are off of the grid  
        self.y['FUV_flux'][ self.y['FUV_flux'] <= 0.0] = -1
        self.y['LW_flux' ][ self.y['LW_flux']  <= 0.0] = -1

        return None


    def interpolate(self, vals, ynames, silence = None):

        if silence == None: # default behavior is ignore off grid values 
            silence = True  # this usualy means radiation should be computed some other way

        single_output = False
        if isinstance(ynames, basestring):
            ynames = [ynames]
            single_output = True

        count = 0
        return_list = [None] * len(ynames)
        for yname in ynames:

            if yname == 'q0' or yname == 'q1':

                return_list[count] =  DataTable.interpolate(self, vals, yname, silence = silence,
                                            flag = 'offgrid', special_flag = 'offgrid',
                                            special_errval = 0.0)

            elif yname == 'FUV_flux' or yname == 'LW_flux':
    
                return_list[count] =  DataTable.interpolate(self, vals, yname, silence = silence,
                                            flag = 'offgrid', special_flag = 'offgrid',
                                            special_errval = -1.0)
            count = count + 1

        if single_output:
            if isinstance(return_list, basestring):
                return return_list
            else:
                return return_list[0]
        else:
            return return_list

    
class StellarYieldsTable(DataTable):

    def __init__(self, yield_type, name = None, manual_table = False):
        """ StellarYieldsDataTable subclass of DataTable

        Given an table type, reads in and constructs data table for 
        stellar yields. Assigns a name to that table as a descriptor
        (name is not used for anything else).

        Args:
           yield_type (str) : Type of table to load. Options are:
               'SNII', 'wind', and 'massive_star'. Former two are NuGrid
               yields for SNII and winds over 1 < M < 25, latter are wind
               only yields for stars 8 < M < 350 from Slemer et. al. 2017.
               Latter is meant to be used only for stars at M > 25.
           name (str, optional) : Name for table, default None. When none,
               a pre-determined descriptive name is used.

        """
        
        if yield_type == 'SNII' and name == None:
            name = 'SNII Stellar Yields Table'
        elif yield_type == 'wind' and name == None:
            name = 'Stellar Winds Stellar Yields Table'
        elif yield_type == 'massive_star' and name == None:
            name = 'Massive Star Stellar Winds Yields Table'
        elif yield_type == None:
            print "Error must set a yiled type as either SNII, wind, or massive_star"
            raise RuntimeError

        
        DataTable.__init__(self, name)

        self.ndim = 2
        self.dim_names = ['mass', 'metallicity']

        self.yield_type = yield_type

        if not manual_table:
            self.read_data()


    def read_data(self, yield_type = None , data_dir = None):
        """ read_data

            Read in data. This is done automatically in initilization, but can
            be repeated in case table needs to be reloaded for whatever reason.
        """

        if yield_type is None:
            yield_type = self.yield_type

        if data_dir == None:
            self.data_dir = install_dir + 'Data/'

        max_col = 87
        if yield_type == 'SNII':
            filename = 'stellar_yields.in'
        elif yield_type == 'wind':
            filename = 'stellar_yields_winds.in'
        elif yield_type == 'massive_star':
            filename = 'stellar_yields_massive_star.in'
            max_col  = 33
       
        tmp_data = np.genfromtxt(self.data_dir + filename, usecols = (0,1), names=True)

        self.x['mass']        = np.unique( tmp_data['M'] )
        self.x['metallicity'] = np.unique( tmp_data['Z'] )

        self.nbins['mass']        = np.size(self.x['mass'])
        self.nbins['metallicity'] = np.size(self.x['metallicity'])

        self._array_size = int( np.prod(self.nbins.values()))

        data   = np.genfromtxt(self.data_dir + filename, usecols=np.arange(2,max_col), names=True)

        for element in list(data.dtype.names):
            self.y[element]   = data[element].reshape(tuple(self.nbins.values()))


        return None


    def interpolate(self, vals, ynames, silence = True):
        """ interpolate
    
        Wrapper around base class interpolation routine to handle 
        edge cases better in this specific instance.
        """
        output = DataTable.interpolate(self, vals, ynames, silence)

        out_of_bounds = False
        if not isinstance(output, basestring):
            for a in output:
                if a == 'offgrid':
                    out_of_bounds = True
                    break
        else:
            if output == 'offgrid':
                out_of_bounds = True
          
        if out_of_bounds:
            output = np.zeros(len(ynames))

        return output
