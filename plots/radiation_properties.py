import numpy as np
import matplotlib.pyplot as plt

import onezone_plot_tools as ptools


def plot_property(name):

    s, m, z  = ptools.star_sample( (1000, 4), [1.0, 100.0], [-4, -3, -2, np.log10(0.017)])



    s    = np.reshape( s, (np.size(z)*np.size(m),))
    prop = [item.properties[name] for item in s]

    #prop = [item.properties[name] for sublist in s for item in sublist]
    prop = np.reshape(prop, (np.size(z),np.size(m)))
    
    

    label_dict = ptools._generate_label_dictionary()

    fig, ax = plt.subplots(1)

    ls = ['-','--','-.',':','-']


    if name == 'E0' or name == 'E1':
        prop = np.array(prop) * 6.242e+11

    if name == 'agb_phase_length':


        for i in np.arange(np.size(z)):
            mtemp = m[prop[i] > 0.0]
            ptemp = prop[i][prop[i] > 0.0]
            pmin = np.min(prop[prop>0.0])
            pmax = np.max(prop[prop>0.0])
            ax.plot(mtemp, ptemp, ls = ls[i], lw = 3, color = 'black', label='Z= %.1E'%(z[i]))
    else:
        for i in np.arange(np.size(z)):
            ax.plot(m, prop[i], ls = ls[i], lw = 3, color = 'black', label='Z= %.1E'%(z[i]))

    ax.set_xlabel(r'Stellar Mass (M$_{\odot}$)')

    if label_dict.has_key(name):
        ax.set_ylabel(label_dict[name])
    else:
        ax.set_ylabel(name)




    if name == 'E0' or name == 'E1':   
        
        ax.set_ylim(np.min(prop), np.max(prop))
        ax.set_xlim(m[0], m[-1])
    elif name == 'agb_phase_length':
        ax.set_ylim(pmin, pmax)
        ax.semilogy()
        ax.set_xlim(1.0, 8.0)
    else:
        ax.set_ylim(np.min(prop), np.max(prop))
        ax.set_xlim(m[0], m[-1])

        ax.semilogy()
    ax.legend(loc='best')
    plt.tight_layout()
    fig.set_size_inches(8,8)
    fig.savefig(name + '.png')
    plt.close()

plot_property('agb_phase_length')
plot_property('L_FUV')
plot_property('L_LW')
plot_property('Q0')
plot_property('Q1')
plot_property('E0')
plot_property('E1')
plot_property('luminosity')
plot_property('lifetime')
plot_property('age_agb')
plot_property('R')
plot_property('Teff')

#plot_property('v_wind')

