import matplotlib.pyplot as plt
from astropy import units as u
from onezone.analysis import analysis_tools
import numpy as np

import os
import sys
#
#

ejection_dict = {}
for e in ['O','Mg','Fe']:
    ejection_dict[e] = 1.0 #(1.0 - 0.95)
for e in ['Ba','N']:
    ejection_dict[e] = (1.0 - 0.6)
for e in ['H','He']:
    ejection_dict[e] = 1.0

fake_yield={}
fake_yield['Eu'] = 8.0E-5 # * (1.0-0.95)

def plot_2D(data, yname, xname):

    xlabel = xname
    ylabel = yname

    if "over" in yname:
        ynum, ydenom = yname.split("_over_")

        if ynum in list(fake_yield.keys()):
            ynum_mass = fake_yield[ynum] * np.ones(np.size(data['t']))
        else:
            ynum_mass = data[ynum+'_mass']*ejection_dict[ynum]

        ydenom_mass = data[ydenom+'_mass']*ejection_dict[ydenom]

        y = analysis_tools.abundance_ratio_array(ynum, ynum_mass,
                                                 ydenom, ydenom_mass,
                                             input_type = "mass")
        ylabel = r"["+ynum+"/"+ydenom+"]"

        select = np.logical_not(np.isnan(y))
        ylim = (-4,4)
        if ydenom == 'H':
            ylim = (-6,0)

    else:
        y = data[yname]
        select = np.array([True]*np.size(y))
        ylim = (np.min(y),np.max(y))

    if "over" in xname:
        xnum,xdenom = xname.split("_over_")
        x = analysis_tools.abundance_ratio_array(xnum, data[xnum+'_mass']*ejection_dict[xnum],
                                                 xdenom, data[xdenom+'_mass']*ejection_dict[xdenom],
                                             input_type = "mass")
        xlabel = r"["+xnum+"/"+xdenom+"]"
        select = select*np.logical_not( np.isnan(x))
        xlim = (-3,3)
        if xdenom == 'H':
            xlim = (-6,0)
    else:
        x = data[xname]
        select = select* np.array([True]*np.size(x))
        xlim = (np.min(x),np.max(x))

    x = x[select]
    y = y[select]

    fig, ax = plt.subplots()
    fig.set_size_inches(6,6)
    ax.plot(x,y,lw=3,color='black',ls='-')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    plt.tight_layout()

    fig.savefig(yname+ "_vs_"+ xname + ".png")
    plt.close()
    return


def panel_plot(xname, yname, data, dim = [0,0], logscale = True,
                xlabel = '-', ylabel = '-', normalize = False, *args, **kwargs):

    if isinstance(xname, str) and isinstance(yname, str):
        xname  = [xname]
        yname  = [yname]
        xlabel = [xlabel]
        ylabel = [ylabel]
    elif isinstance(xname, str) and not isinstance(yname, str):
        nplots = len(yname)
        xname  = [xname] * nplots
        xlabel = [xlabel] * nplots
    elif not isinstance(xname, str) and not isinstance(yname, str):
        # assume both are lists of names
        if not ( len(xname) == len(yname) ):
            print("lengths of x and y axis names must match")
            raise TypeError

        nplots = len(yname)

    if isinstance(ylabel, str):
        ylabel = [ylabel] * nplots

    #
    # figure out the dimensions
    # 
    if dim == [0,0]:
        dim = np.ceil(np.sqrt(nplots))
        dim = [dim, dim]
    elif len(dim) == 0:
        dim = [dim, dim]
    
    fig, ax = plt.subplots(dim[0], dim[1])

    if not hasattr(ax, '__iter__'):
        ax = [ax]

    i = 0

    for a, x, y in zip(ax.flatten(), xname, yname):

        if '/' in x:
            xsplit = x.split('/')
            if normalize:
                xdata = analysis_tools.normalize_abundance_ratio( (xsplit[0], data[xsplit[0]] * u.Msun),
                                                   (xsplit[1], data[xsplit[1]] * u.Msun),
                                                   input_type = 'mass' )
            else:
                xdata = np.log10( data[xsplit[0]] / data[xsplit[1]] )

            if logscale:
                logscale = 'y'
        else:
            xdata = data[x]

        if '/' in y:
            ysplit = y.split('/')
            if normalize:
                ydata = analysis_tools.normalize_abundance_ratio( (ysplit[0], data[ysplit[0]] * u.Msun),
                                                   (ysplit[1], data[ysplit[1]] * u.Msun),  
                                                   input_type = 'mass' )

            else:
                ydata = np.log10( data[ysplit[0]] / data[ysplit[1]] )

            if logscale == 'y':
                logscale = False
        else:
            ydata = data[y]
        

        if '/' in y:
            a.plot( xdata, ydata, label = ylabel[i], **kwargs)
            a.set_ylabel(r'log X/' + ysplit[1])
            a.legend(loc='best', frameon=False)
        else:
            a.plot( xdata, ydata, **kwargs)
            a.set_ylabel(ylabel[i])

        a.set_xlabel(xlabel[i])

        if logscale == 'x':
            a.semilogx()
        elif logscale == 'y':
            a.semilogy()
        elif logscale:
            a.loglog()

        i = i + 1
                

    return fig, ax
        


if __name__ == "__main__":

    print(sys.argv)
    if sys.argv[1] == "2D":
        data = np.genfromtxt( str(sys.argv[2]), names = True)

        plot_2D( data, str(sys.argv[3]), str(sys.argv[4]))
