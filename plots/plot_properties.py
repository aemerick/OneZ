def plot_properties(abundances, yield_mode, fractional=False):


    s, m, z  = ptools.star_sample( (200, 4), [1.0, 100.0],
                                   [-4, -3, -2, np.log10(0.017)],
                                   abundances=abundances)

    s  = np.r
