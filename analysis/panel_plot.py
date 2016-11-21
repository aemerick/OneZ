import matplotlib.pyplot as plt
from astropy import units as u
from onezone.analysis import analysis_tools
import numpy as np
#
#
def panel_plot(xname, yname, data, dim = [0,0], logscale = True,
                xlabel = '-', ylabel = '-', normalize = False, *args, **kwargs):

    if isinstance(xname, basestring) and isinstance(yname, basestring):
        xname  = [xname]
        yname  = [yname]
        xlabel = [xlabel]
        ylabel = [ylabel]
    elif isinstance(xname, basestring) and not isinstance(yname, basestring):
        nplots = len(yname)
        xname  = [xname] * nplots
        xlabel = [xlabel] * nplots
    elif not isinstance(xname, basestring) and not isinstance(yname, basestring):
        # assume both are lists of names
        if not ( len(xname) == len(yname) ):
            print "lengths of x and y axis names must match"
            raise TypeError

        nplots = len(yname)

    if isinstance(ylabel, basestring):
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
        
