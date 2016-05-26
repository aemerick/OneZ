"""
Will contain functions and classes for loading in needed
data and data tables
"""


class DataTable:

    def __init__(self, name):
        self.name = name

    def read_data(self):
        pass

    def _interpolate(self, y, vals, val_arrays):

        c, id = _interpolation_coefficients(vals, val_arrays)

        n       = len(vals) # number of dimensions

        if n == 1:
            yval = (1.0 - c[0]) * y[id[0]] + (c[0]) * y[id[0]]

        elif n == 2:
            i,j = id

            yval = (1.0 - c[0]) * (1.0 - c[1]) * y[i  ][j  ] +\
                   (1.0 - c[0]) * (      c[1]) * y[i  ][j+1] +\
                   (      c[0]) * (      c[1]) * y[i+1][j+1] +\
                   (      c[0]) * (1.0 - c[1]) * y[i+1][j  ]

        elif n == 3:
            i,j,k = id

            yval = (1.0 - c[0])*(1.0 - c[1])*(1.0 - c[2]) * y[i  ][j  ][k  ] +\
                   (1.0 - c[0])*(      c[1])*(1.0 - c[2]) * y[i  ][j+1][k  ] +\
                   (      c[0])*(      c[1])*(1.0 - c[2]) * y[i+1][j+1][k  ] +\
                   (      c[0])*(1.0 - c[1])*(1.0 - c[2]) * y[i+1][j  ][k  ] +\
                   (1.0 - c[0])*(1.0 - c[1])*(      c[2]) * y[i  ][j  ][k+1] +\
                   (1.0 - c[0])*(      c[1])*(      c[2]) * y[i  ][j+1][k+1] +\
                   (      c[0])*(      c[1])*(      c[2]) * y[i+1][j+1][k+1] +\
                   (      c[0])*(1.0 - c[1])*(      c[2]) * y[i+1][j  ][k+1]

        elif n > 3:
            print "We do not support n > 3 dimensional interpolation"
            raise RuntimeError

        return yval


    def _interpolation_coefficients(self, vals, val_arrays):
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

        return coeff, indeces

    def _linear_interpolation_coefficients(self, x, xarray):

        # check bounds in this dimension
        if (x < xarray[0] or x > xarray[-1]):
            print x, xarray[0], xarray[-1]

            raise ValueError

        i = np.abs( x - xarray).argmin()
        if ( x < xarray[i]):
          i = i - 1

        t = (x - xarray[i]) / (xarray[i+1] - xarray[i])

        return t, i        


# in rad data class, do 
# try:
#      interpolation coeff
# except ValueError:
#      off grid, use black body


class RadData(DataTable):

    def __init__(self):
        DataTable.__init__("Radiation data table")

        self.

