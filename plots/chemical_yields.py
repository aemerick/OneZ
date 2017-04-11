import numpy as np
import matplotlib.pyplot as plt

import onezone_plot_tools as ptools
from collections import OrderedDict

from onezone import config as config

def plot_yields(abundances, yield_mode, fractional=False, IMF_weighted=False):


    s, m, z  = ptools.star_sample( (200, 4), [1.0, 100.0],
                                   [-4, -3, -2, np.log10(0.017)],
                                   abundances=abundances)


    yields = OrderedDict()
    s    = np.reshape( s, (np.size(z)*np.size(m),))

    if yield_mode == 'wind':
        for a in abundances:
            yields[a] = [item.wind_ejecta_masses()[a] for item in s]
            yields[a] = np.reshape(yields[a], (np.size(z),np.size(m)))   
    elif yield_mode == 'SN':

        for a in abundances:
            yields[a] = [item.return_sn_ejecta_masses()[a] for item in s]
            yields[a] = np.reshape(yields[a], (np.size(z),np.size(m)))


    
    for nz in np.arange(np.size(z)):
        fig, ax = plt.subplots(1)

        ls_count = 0 ; color_count = 0
        colors = ['black', 'blue', 'green', 'red', 'purple', 'orange']
        ls     = ['-','--','-.',':','-']

        for a in abundances:
            norm = 1.0

            if fractional:
                norm = yields['m_tot'][nz]

            if IMF_weighted:
                weights = config.zone.imf.imf(m)
                yields[a][nz] = yields[a][nz] * weights / np.sum(weights)

            ax.plot(m, yields[a][nz] / norm, ls = ls[ls_count], color = colors[color_count],
                       label = a, lw = 3)

            color_count = color_count + 1
            if color_count >= len(colors):
                color_count = 0
                ls_count = ls_count + 1

        ax.set_xlabel(r'Stellar Mass (M$_{\odot}$)')

        if fractional:
            ylabel = r'Fractional Ejected Mass (M$_{\odot}$)'
        elif IMF_weighted:
            ylabel = r'IMF Weighted Ejected Mass (M$_{\odot}$)'
        else:
            ylabel = r'Ejected Mass (M$_{\odot}$)'
 
        ax.set_ylabel(ylabel)
        ax.semilogy()
        ax.set_ylim(1.0E-9, 1.0E2)
        ax.legend(loc='lower right', ncol=4)
        plt.tight_layout()
        fig.set_size_inches(8,8)


        fractional_out = ''
        if fractional:
            fractional_out = 'fractional_'
        elif IMF_weighted:
            fractional_out = 'IMF_weighted_'

        if yield_mode == 'wind':
            yield_out = 'wind'
        elif yield_mode == 'SN':
            yield_out = 'SN'

        fig.savefig('./yields/' + fractional_out + yield_out + '_yields_z=%.4f.png'%(z[nz]))


        plt.close()


abundances = OrderedDict()
abundances['m_tot'] = 1.0
abundances['m_metal'] = 0.0

for a in ['H','He','C','N','O','Mg','Si','Mn','Fe','Ni','Y','Ba','Eu']:
    abundances[a] = 0.0

plot_yields(abundances, yield_mode = 'wind')
plot_yields(abundances, yield_mode = 'wind', IMF_weighted=True)

plot_yields(abundances, yield_mode = 'SN')
plot_yields(abundances, yield_mode = 'SN', IMF_weighted = True)




plot_yields(abundances, yield_mode = 'wind')
plot_yields(abundances, yield_mode = 'wind', fractional = False)

plot_yields(abundances, yield_mode = 'SN')
plot_yields(abundances, yield_mode = 'SN', fractional = False)


