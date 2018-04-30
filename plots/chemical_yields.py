import numpy as np
from galaxy_analysis.plot.plot_styles import *

import matplotlib.pyplot as plt
import onezone_plot_tools as ptools

from collections import OrderedDict

from onezone import config as config

def plot_yields(abundances, yield_mode, fractional=False, IMF_weighted=False):


    s, m, z  = ptools.star_sample( (1000), [1.0, 100.0],
                                   4.3E-4, #, -4, -3, -2, np.log10(0.017)],
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



    labels = {}
    for a in abundances:
        labels[a] = a
    if 'm_tot' in labels.keys():
        labels['m_tot'] = 'Total Mass'
    if 'm_metal' in labels.keys():
        labels['m_metal'] = 'Total Metals'
    
    for nz in np.arange(np.size(z)):
        fig, ax = plt.subplots(1)

        ls_count = 0 ; color_count = 0
        colors = ['black','C0','C1','C2','C3','C4','C5','C6','C7','C8','C9']
        ls     = ['-','--','-.',':','-']

        for a in abundances:
            norm = 1.0

            if fractional:
                norm = yields['m_tot'][nz]

            if IMF_weighted:
                weights = config.zone.imf.imf(m)
                yields[a][nz] = yields[a][nz] * weights / np.sum(weights)

            _ls = ls[ls_count]
            if a in ['Sr','Y','Ba']:
                _ls = '--'
           

            ax.plot(m, yields[a][nz] / norm, ls = _ls, color = colors[color_count],
                       label = labels[a], lw = line_width)

            color_count = color_count + 1
            if color_count >= len(colors):
                color_count = 0
                ls_count = ls_count + 1

        ax.set_xlabel(r'Stellar Mass (M$_{\odot}$)')
        ax.set_xlim(0.0,25.0)

        if fractional:
            ylabel = r'Fractional Ejected Mass (M$_{\odot}$)'
        elif IMF_weighted:
            ylabel = r'IMF Weighted Ejected Mass (M$_{\odot}$)'
        else:
            ylabel = r'Ejected Mass (M$_{\odot}$)'
 
        ax.set_ylabel(ylabel)
        ax.semilogy()
        ax.set_ylim(1.0E-15, 1.0)
        ax.legend(loc='upper right', ncol=3)
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

        ax.set_ylim(ax.get_ylim())
        ax.plot([8.0,8.0],ax.get_ylim(), lw = line_width, ls = '--', color = 'black')
        ax.plot([25.0,25.0],ax.get_ylim(), lw = line_width, ls = '-', color = 'black')
        plt.tight_layout()
        fig.savefig('./yields/' + fractional_out + yield_out + '_yields_z=%.4f.png'%(z[nz]))


        plt.close()


if __name__ == "__main__":

    abundances = OrderedDict()
    abundances['m_tot']   = 1.0
    abundances['m_metal'] = 0.0

#    for a in ['C','N','O','Na''Mg','Si','S','Ca','Mn','Fe','Ni','As','Sr','Y','Ba']:
    for a in ['C','N','O','Mg','S','Sr','Y','Ba']:
        abundances[a] = 0.0

    plot_yields(abundances, yield_mode = 'wind')
    plot_yields(abundances, yield_mode = 'wind', fractional = True)

    plot_yields(abundances, yield_mode = 'SN')
    plot_yields(abundances, yield_mode = 'SN', fractional = True)

    plot_yields(abundances, yield_mode = 'wind')
    plot_yields(abundances, yield_mode = 'wind', IMF_weighted=True)

    plot_yields(abundances, yield_mode = 'SN')
    plot_yields(abundances, yield_mode = 'SN', IMF_weighted = True)

    plot_yields(abundances, yield_mode = 'wind')
    plot_yields(abundances, yield_mode = 'wind', fractional = False)

    plot_yields(abundances, yield_mode = 'SN')
    plot_yields(abundances, yield_mode = 'SN', fractional = False)


