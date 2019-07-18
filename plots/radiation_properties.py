
from galaxy_analysis.plot.plot_styles import *

import numpy as np
#from matplotlib import rc

#fsize = 17
#rc('text', usetex=False)
#rc('font', size=fsize)#, ftype=42)
#line_width = 3
#point_size = 30

rc('font', size = 24)

import matplotlib.pyplot as plt

import onezone_plot_tools as ptools
from onezone import config as config

def plot_property(name, IMF_weighted = False, xlim = None, ylim = None):

    s, m, z  = ptools.star_sample( (1000, 4), [1.0, 100.0], [-4, -3, -2, np.log10(0.017)])

    s    = np.reshape( s, (np.size(z)*np.size(m),))
    prop = [item.properties[name] for item in s]

    #prop = [item.properties[name] for sublist in s for item in sublist]
    prop = np.reshape(prop, (np.size(z),np.size(m)))
    
    

    label_dict = ptools._generate_label_dictionary()

    fig, ax = plt.subplots(1)
    fig.set_size_inches(8,8)

    ls = ['-','--','-.',':','-']


    if name == 'E0' or name == 'E1':
        prop = np.array(prop) * 6.242e+11

    if name == 'agb_phase_length':


        for i in np.arange(np.size(z)):
            mtemp = m[prop[i] > 0.0]
            ptemp = prop[i][prop[i] > 0.0]
            pmin = np.min(prop[prop>0.0])
            pmax = np.max(prop[prop>0.0])
            ax.plot(mtemp, ptemp, ls = ls[i], lw = 2.0, color = 'black', label='Z= %.1E'%(z[i]))
    else:
    

        for i in np.arange(np.size(z)):
            y = 1.0 * prop[i]
            if IMF_weighted:
                y = y * config.zone.imf.imf(m) / np.sum(config.zone.imf.imf(m)) 

            ax.plot(m, y, ls = ls[i], lw = 2.0, color = 'black', label='Z= %.1E'%(z[i]))

    ax.set_xlabel(r'Stellar Mass (M$_{\odot}$)')

    if name in label_dict:
        ax.set_ylabel(label_dict[name])
    else:
        ax.set_ylabel(name)

    
    if name == 'E0' or name == 'E1':   
        if (ylim is None): ax.set_ylim(np.min(prop), np.max(prop))
        if (xlim is None): ax.set_xlim(m[0], m[-1])
    elif name == 'agb_phase_length':
        if (ylim is None): ax.set_ylim(pmin, pmax)
        ax.semilogy()
        if (xlim is None): ax.set_xlim(1.0, 8.0)
    else:
        if (ylim is None): ax.set_ylim(np.min(prop), np.max(prop))
        if (xlim is None): ax.set_xlim(m[0], m[-1])

        ax.semilogy()

    if not (xlim is None):
        ax.set_xlim(xlim)
    if not (ylim is None):
        ax.set_ylim(ylim)

    ax.minorticks_on()
    ax.legend(loc='best')
    plt.tight_layout()
    fig.savefig(name + '.png')
    plt.close()

plot_property('agb_phase_length')
plot_property('L_FUV', xlim = [1,100], ylim = [1.0E35, 1.0E39])
plot_property('L_LW', xlim = [1,100], ylim = [1.0E35, 1.0E39])
plot_property('Q0', xlim = [1,100], ylim = [1.0E45, 1.0E50])
plot_property('Q1', xlim = [1,100], ylim = [1.0E45, 1.0E50])
plot_property('E0', xlim = [1,100], ylim = [14,34])
plot_property('E1', xlim = [1,100], ylim = [14,34])
plot_property('luminosity')
plot_property('lifetime')
plot_property('age_agb')
plot_property('R')
plot_property('Teff')

#plot_property('v_wind')

