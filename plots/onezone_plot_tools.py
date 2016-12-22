import numpy as np
from onezone import star

def star_sample(nstar, m, z, logM = False, logZ = True,
                abundances = {'m_tot' : 1.0}):
    """
    Returns sample of stars over a mass and metallicity range. 
    nstar is tuple with number of stars in each dimension

    """

    if np.size(m) == 2:
        if logM:
            m = np.logspace( m[0], m[1], nstar[0])
        else:
            m = np.linspace( m[0], m[1], nstar[0])

    elif not np.iterable(m):
        m = [m]

    if np.size(z) == 2:
        if logZ:
            z = np.logspace( z[0], z[1], nstar[1])
        else: 
            z = np.linspace( z[0], z[1], nstar[1])    

    elif np.size(z) > 2 and logZ:
        z = [10.0**(x) for x in z]

    elif not np.iterable(z):
        z = [z]

    total_count = np.size(m) * np.size(z)


    s = [None] * total_count

    for i in np.arange(np.size(z)):
        for j in np.arange(np.size(m)):
            s[i*np.size(m) + j] = star.Star(M = m[j], Z = z[i], id= j + i * np.size(m),
                                            abundances=abundances)

    s = np.reshape(s, (np.size(z), np.size(m)))


    return s, m, z



def _generate_label_dictionary():

    x = {'L_LW'       : r'L$_{\rm LW}$ (erg s$^{-1}$)',
         'L_FUV'      : r'L$_{\rm FUV}$ (erg s$^{-1}$)',
         'Q0'         : r"Q$_{\rm 0}$ (s$^{-1}$)",
         'Q1'         : r"Q$_{\rm 1}$ (s$^{-1}$)",
         'luminosity' : r"L$_{\rm bol}$ (erg s$^{-1}$)",
         'v_wind'     : r"v$_{\rm wind}$ (km s$^{-1}$)",
         'E0'         : r"<E$_{\rm HI}$> (eV)",
         'E1'         : r"<E$_{\rm HeI}$> (eV)"}


    return x
