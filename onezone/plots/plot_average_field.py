from galaxy_analysis.plot.plot_styles import *
import numpy as np
import matplotlib.pyplot as plt
import sys, glob
from scipy.integrate import cumtrapz

def compute_property(fpath, property, integrate = False, normalize = None):

    # load file to get names
    files = np.sort(glob.glob(fpath + 'run????_summary_output.txt'))
    f0    = np.genfromtxt(files[0], names = True)
#    col   = np.where( np.array(f0.dtype.names) == property)[0][0]

    # resample output - don't assume output has same dt (and dt will change)
    sample_times = np.arange(0.0, 150.0, 1.0)

    nfiles   = np.size(files)
    npoints  = np.size(sample_times)
    all_data = np.zeros( (nfiles, npoints) )


    for i in np.arange(nfiles):

        temp = np.genfromtxt(files[i], names = True)

        select = temp['M_star_o'] > 0

        x    = temp['t'][select] - np.min(temp['t'][select])

        # renormalize time
#        x    = x - np.min(x[M_s>0])

        norm = 1.0
        if not (normalize is None):
            if normalize == 'M_o':
                norm = 1.0 / temp['M_star_o'][select]
            else:
                norm = normalize

        if 'L_' in property:
            y = temp['int_mass_' + property.replace('_','')] +\
                temp['high_mass_' + property.replace('_','')] +\
                temp['vhigh_mass_' + property.replace('_','')]
            y = y[select]
        else:
            y    = temp[property][select]

        newy = np.interp( sample_times, x, y)
        if np.size(norm) > 1:
            norm = np.interp(sample_times, x, norm)

        all_data[i] = newy
        x           = 1.0 * sample_times

        if integrate:
#            temp_data = np.zeros(npoints)
#            for j in np.arange(1,npoints):
#                temp_data[j] = np.trapz( [all_data[i][j-1],all_data[i][j]], [x[j-1], x[j]]) + temp_data[j-1]
            all_data[i][:-1] = cumtrapz( all_data[i], x = x)
            all_data[i][-1] = all_data[i][-2]

        all_data[i] = all_data[i] * norm

    stats = {}
    stats['median'] = np.median(all_data, axis = 0)
    stats['q1']     = np.percentile(all_data, 25, axis = 0)
    stats['q3']     = np.percentile(all_data, 75, axis = 0)
    stats['min']    = np.min(all_data, axis = 0)
    stats['max']    = np.max(all_data, axis = 0)
    stats['t'] = sample_times

    return stats

def plot_property(fpath, property, ylabel = None, integrate = False, normalize = None):

    stats = compute_property(fpath, property, ylabel, integrate, normalize)


    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)


    ax.plot(stats['t'], stats['median'], lw = 3, color = 'black', ls = '--')
    ax.plot(stats['t'], stats['q1'],     lw = 1, color = 'black', ls = '-')
    ax.plot(stats['t'], stats['q3'],     lw = 1, color = 'black', ls = '-')
    ax.fill_between(stats['t'], stats['q1'], stats['q3'], color = 'black', alpha = 0.5, interpolate = True)
    ax.fill_between(stats['t'], stats['min'], stats['max'], color = "black", alpha = 0.25, interpolate = True)

    ax.set_xlabel(r'Time (Myr)')

    if ylabel is None:
        ylabel = property

    ax.set_ylabel(ylabel)
    ax.semilogy()
    plt.tight_layout()
    plt.minorticks_on()
    fig.savefig(fpath + property + '_median_evolution.png')

    plt.close()
    return


if __name__ == "__main__":
    property = 'L_Q0'
    integrate = True
    ylabel   = 'HI Ionizing Luminosity'
    normalize = 3.162247E13 # / 6.2415E11   #  Myr -> s

    if len(sys.argv) > 1:
        property = string(sys.argv[1])
        integrate = False
        ylabel   = None
        normalize = None

    plot_property('./', property, integrate = integrate, ylabel = ylabel, normalize = normalize)
