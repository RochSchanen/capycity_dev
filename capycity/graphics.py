# file: graphics.py

# from package: "https://matplotlib.org/"
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import contourf
from matplotlib import pyplot as pp

# from package: "https://numpy.org/"
from numpy import full
from numpy import arange
from numpy import array
from numpy import copy
from numpy import empty
from numpy import square

from numpy import concatenate
from numpy import log

# from numpy import save
# from numpy import load

from capycity.config import DATATYPE, ZER, MAX

###################################################### POTENTIAL CONTOUR FILLED

def plotPotentialContourFill(solver, n = 30, *args, pdfdoc = None, **kwargs):
    """  """

    # some parameters
    plot_ticks, contour_ticks = 0.1, 0.1

    # create figure
    fg = pp.figure("Contour")
    
    # set A4 paper dimensions
    fg.set_size_inches(8.2677, 11.6929)
    
    # create title
    fg.suptitle(f"Potential Distribution")
    
    # create axis
    w, h, k = [*list(array([1, 1 / 1.4143])*0.7), 0.03]
    x, y = (1-w-k)/3, (1-h)/2
    ax = fg.add_axes([x, y, w, h]) 
    
    # x-coordinates
    d = solver.l/(solver.n-1)
    l = d*(solver.n-1)/2
    X = arange(-l, l+d, d)
    
    # y-coordinates
    Y = solver.D[1:-1, 1:-1] / MAX
    
    # make contour plot
    c = [pp.cm.inferno, pp.cm.terrain][0]
    mapargs = {"cmap":c, "vmin":0.0, "vmax":1.0}
    pl = ax.contourf(X, X, Y, n, *args, **mapargs, **kwargs)
    
    # adjust ticks
    ax.set_xticks(list(arange(-0.5, 0.5 + plot_ticks, plot_ticks)))
    ax.set_yticks(list(arange(-0.5, 0.5 + plot_ticks, plot_ticks)))
    
    # add decors
    for g in solver.G:
        ax.add_artist(g.contour_decor())
    
    # add new axis and color bar 
    ax = fg.add_axes([x + w + x, y, k, h])
    cb = fg.colorbar(pl, cax = ax)
    cb.set_ticks(list(arange(0, 1 + contour_ticks, contour_ticks)))
    
    # export figure to pdf
    if pdfdoc: pdfdoc.savefig(fg)
    
    # done
    return

####################################################### POTENTIAL CROSS SECTION

def plotPotentialCrossSection(
    solver,             # solver with data set
    *args,              # arguments passed to "plot"
    name   = None,      # figure name to add the plot to
    pdfdoc = None,      # pdf document to add the figure to
    **kwargs):          # more arguments passed to "plot"
    """ explain """

    # some default parameters
    plot_ticks = 0.1

    # set generic name if undefined
    if name is None:
        # Further plot will be added to "plotCrossSection"
        name = "plotCrossSection"

    if not pp.fignum_exists(name):

        # create figure
        fg = pp.figure(name)

        # set A4 paper dimensions
        fg.set_size_inches(8.2677, 11.6929)

        # create title
        fg.suptitle(f"Potential Cross Section")

        # create axis
        w, h = array([1, 1 / 1.4143])*0.7
        x, y = (1-w)/2, (1-h)/2
        ax = fg.add_axes([x, y, w, h])

        # x-coordinates
        d = solver.l/(solver.n-1)
        l = d*(solver.n-1)/2
        X = arange(-l, l+d, d)

        # y-coordinates
        i = int(solver.n // 2 + 1)
        Y = solver.D[i, 1:-1] / MAX

        # make plot
        pl = ax.plot(X, Y, *args, **kwargs)

        # adjust ticks
        ax.set_xticks(list(arange(-0.5, 0.5 + plot_ticks, plot_ticks)))
        ax.set_yticks(list(arange(0.0, 1.0 + plot_ticks, plot_ticks)))

        # add decors
        for g in solver.G:
            ax.add_artist(g.section_decor())

    else:

        # select figure
        fg = pp.figure(name)

        # select axes
        ax = fg.get_axes()[0]

        # x-coordinates
        d = solver.l/(solver.n-1)
        l = d*(solver.n-1)/2
        X = arange(-l, l+d, d)

        # y-coordinates
        i = int(solver.n // 2 + 1)
        Y = solver.D[i, 1:-1] / MAX

        # make plot
        pl = ax.plot(X, Y, *args, **kwargs)

    # export figure to pdf
    if pdfdoc: pdfdoc.savefig(fg)

    # done
    return

####################################################### POTENTIAL THEORY

def plotTheory(
    r1, r2,             # inner and outer radius
    v1, v2,             # inner and outer potential
    *args,              # arguments passed to "plot"
    name   = None,      # figure name to add the plot to
    pdfdoc = None,      # pdf document to add the figure to
    **kwargs):          # more arguments passed to "plot"
    """ explain """

    # some default parameters
    plot_ticks = 0.1

    # set generic name if undefined
    if name is None:
        # Further plot will be added to "plotCrossSection"
        name = "plotCrossSection"

    if not pp.fignum_exists(name):

        # create figure
        fg = pp.figure(name)

        # set A4 paper dimensions
        fg.set_size_inches(8.2677, 11.6929)

        # create title
        fg.suptitle(f"Potential Cross Section")

        # create axis
        w, h = array([1, 1 / 1.4143])*0.7
        x, y = (1-w)/2, (1-h)/2
        ax = fg.add_axes([x, y, w, h])

        # x-coordinates
        XL = arange(-r2, -r1+(r2-r1)/50, (r2-r1)/50)
        XR = arange(+r1, +r2+(r2-r1)/50, (r2-r1)/50)

        # y-coordinates
        YL = log(-XL/r2)/log(r1/r2)
        YR = log(+XR/r2)/log(r1/r2)

        # make plot
        pl = ax.plot(XR, YR, *args, **kwargs)
        pl = ax.plot(XL, YL, *args, **kwargs)

        # adjust ticks
        ax.set_xticks(list(arange(-0.5, 0.5 + plot_ticks, plot_ticks)))
        ax.set_yticks(list(arange(0.0, 1.0 + plot_ticks, plot_ticks)))

    else:

        # select figure
        fg = pp.figure(name)

        # select axes
        ax = fg.get_axes()[0]

        # x-coordinates
        XL = arange(-r2, -r1+(r2-r1)/50, (r2-r1)/50)
        XR = arange(+r1, +r2+(r2-r1)/50, (r2-r1)/50)

        # y-coordinates
        YL = log(-XL/r2)/log(r1/r2)
        YR = log(+XR/r2)/log(r1/r2)

        # make plot
        pl = ax.plot(XR, YR, *args, **kwargs)
        pl = ax.plot(XL, YL, *args, **kwargs)

    # export figure to pdf
    if pdfdoc: pdfdoc.savefig(fg)

    # done
    return



def addplot(
    *args,              # arguments passed to "plot"
    name        = None, # figure name to add the plot to
    pdfdoc      = None, # pdf document to add the figure to
    **kwargs,           # more arguments passed to "plot"
    ):         

    # some default parameters
    plot_ticks = 0.1

    # set generic name if undefined
    if name is None:
        # Further plot will be added to "plotCrossSection"
        name = "plotCrossSection"

    if not pp.fignum_exists(name):

        # create figure
        fg = pp.figure(name)

        # set A4 paper dimensions
        fg.set_size_inches(8.2677, 11.6929)

        # create title
        fg.suptitle(f"Potential Cross Section")

        # create axis
        w, h = array([1, 1 / 1.4143])*0.7
        x, y = (1-w)/2, (1-h)/2
        ax = fg.add_axes([x, y, w, h])

        # x-coordinates
        # XL = arange(-r2, -r1+(r2-r1)/50, (r2-r1)/50)
        # XR = arange(+r1, +r2+(r2-r1)/50, (r2-r1)/50)
        d = solver.l/(solver.n-1)
        l = d*(solver.n-1)/2
        X = arange(-l, l+d, d)

        # y-coordinates
        YL = log(-X/r2)/log(r1/r2)
        YR = log(+X/r2)/log(r1/r2)
        i = int(solver.n // 2 + 1)
        Y = solver.D[i, 1:-1] / MAX

        # make plot
        pl = ax.plot(X, YR, *args, **kwargs)
        pl = ax.plot(X, YL, *args, **kwargs)
        pl = ax.plot(X, Y, *args, **kwargs)

        # adjust ticks
        ax.set_xticks(list(arange(-0.5, 0.5 + plot_ticks, plot_ticks)))
        ax.set_yticks(list(arange(0.0, 1.0 + plot_ticks, plot_ticks)))

    else:

        # select figure
        fg = pp.figure(name)

        # select axes
        ax = fg.get_axes()[0]

        # x-coordinates
        XL = arange(-r2, -r1+(r2-r1)/50, (r2-r1)/50)
        XR = arange(+r1, +r2+(r2-r1)/50, (r2-r1)/50)
        d = solver.l/(solver.n-1)
        l = d*(solver.n-1)/2
        X = arange(-l, l+d, d)

        # y-coordinates
        YL = log(-XL/r2)/log(r1/r2)
        YR = log(+XR/r2)/log(r1/r2)
        i = int(solver.n // 2 + 1)
        Y = solver.D[i, 1:-1] / MAX

        # make plot
        pl = ax.plot(XR, YR, *args, **kwargs)
        pl = ax.plot(XL, YL, *args, **kwargs)
        pl = ax.plot(X, Y, *args, **kwargs)

    # export figure to pdf
    if pdfdoc: pdfdoc.savefig(fg)

    # done
    return
