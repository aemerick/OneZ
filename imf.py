# -- external --
import numpy as np

class IMF(object):

    def __init__(self, M_min = 1.0 , M_max = 100.0):

        self._M_min = M_min
        self._M_max = M_max # min and max M of IMF

    def integrate(self, m_a, m_b):
        raise NotImplementedError

    def sample(self, N = None, M = None, npoints = 1000,
                     regenerate_table = False):
        """
        Sample the given imf using a cumulative distribution
        sampled at npoints. Sample until N stars are obtained
        or M mass is reached. 
        """

        # if needed, regenerated tabulated imf
        if regenerate_table or (not hasattr(self, '_tabulated_imf')):
            self._tabulate_imf(npoints)

        def _find_bin(x, array):
            i = np.abs(array - x).argmin()

            if x < array[i] and i > 0:
                i = i -1
            if x < array[i] and i > 0:
                i = i -1
            if x < array[i] and i > 0:
                print(i+2, array[i+2], x)
                print(i, array[i], x)
                print("Failure finding bin")

            return i

        if (not (N is None)) and M is None:
            random_numbers = np.random.rand(N)
            stars = np.zeros(N)

            for i in np.arange(N):
                bin = _find_bin(random_numbers[i], self._tabulated_imf)
                stars[i] = 10.0**(bin * self._tabulated_dm + self._m_o)

        elif (not (M is None)) and N is None:

           stars = np.zeros(int(M / self._M_min))
           i = -1
           total_mass = 0.0

           while total_mass <= M:
               i = i + 1
               rnum = np.random.rand()

               bin = _find_bin(rnum, self._tabulated_imf)
               stars[i] = 10.0**(bin * self._tabulated_dm + self._m_o)
               total_mass = np.sum(stars)

           # remove trailing stars
           stars = stars[0:i]

           # check if we need to remove last star
           if np.size( stars ) > 1:
               if np.abs( (total_mass) - M) < np.abs(np.sum(stars[:i]-M)):
                   stars = stars[0 : i -1]

        return stars


    def _tabulate_imf(self, npoints):

        # generate binned imf
        dm = np.log10(self._M_max / self._M_min) / (1.0*(npoints-1))

        m_o = np.log10(self._M_min)

        m   = 10.0**(m_o + np.arange(0,npoints) * dm)

        # 
        IMF_vals = np.cumsum( self.imf(m) )
        IMF_vals = IMF_vals / ( IMF_vals[-1] * 1.0)

        self._tabulated_dm  = dm
        self._tabulated_m   = m
        self._tabulated_imf = IMF_vals
        self._m_o           = m_o

        return
        
    def imf(self, M):
        """
        Evaluate the imf at a given mass. This is 
        where the functional form of the IMF 
        is housed
        """
        raise NotImplementedError

    def _retabulate_imf(self):
        """
        Function to retabulate IMF if it exists already and 
        if a property is changed.
        """

        # haven't tabulated yet, handle this at first call to IMF
        if not hasattr(self, '_tabulated_imf'):
            return

        npoints = np.size(self._tabulated_imf)
        self._tabulate_imf(npoints)

        return

    @property
    def M_min(self):
        return self._M_min

    @M_min.setter
    def M_min(self, value):
        self._M_min = value
        self._retabulate_imf()
        return

    @property
    def M_max(self):
        return self._M_max

    @M_max.setter
    def M_max(self, value):
        self._M_max = value
        self._retabulate_imf()
        return

#    @classmethod
#    def _M_min(cls):
#        return self.M_min

#    @classmethod
#    def _M_max(cls):
#        return self.M_max

class kroupa(IMF):

    def __init__(self, alpha = [0.3,1.3,2.3],
                       M_min = 0.1, M_max =120.0 ):


        self.alpha = alpha
        self._M_min = M_min
        self._M_max = M_max

    def imf(self, M):
        """
        Evaluates the imf at a point m
        """
        
        print(M)
        M, scalar_input = _check_scalar_input(M)

        print(M)
        low_mass = (M <= 0.08)
        mid_mass = (M  > 0.08)*(M <= 0.5)
        salpeter = (M  > 0.5 )

        dNdM = np.zeros(np.shape(M))

        
        dNdM[low_mass] = M[low_mass]**(- self.alpha[0])
        dNdM[mid_mass] = M[mid_mass]**(- self.alpha[1])
        dNdM[salpeter] = M[salpeter]**(- self.alpha[2])

        return _check_scalar_output(dNdM, scalar_input)

 

class salpeter(IMF):

    def __init__(self, alpha = 1.35, M_min = 1.0, M_max = 120.0):

        self._M_min = M_min
        self._M_max = M_max
        self.alpha = alpha

    def imf(self, M):
        return M**(-self.alpha)

def _check_scalar_input(x):

    x = np.asarray(x)
    scalar_input = False
    if x.ndim == 0:
        x = x[None]
        scalar_input = True

    return x, scalar_input


def _check_scalar_output(x, scalar_input):

    if scalar_input:
        return np.squeeze(x)
    else:
        return x
