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

# from numpy import save
# from numpy import load

from numpy import uint32 as DATATYPE
ZER, MAX = array([0, -1], DATATYPE)

########################################################### POTENTIAL CONTOUR F

def plotPotentialContourFill(solver, n = 30, pdfdoc = None):
    """  """
    # some parameters
    plot_ticks, contour_ticks = 0.1, 0.1
    # create figure
    fg = pp.figure()
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
    args = {"cmap":c, "vmin":0.0, "vmax":1.0}
    pl = ax.contourf(X, X, Y, n, **args)
    # adjust ticks
    ax.set_xticks(list(arange(-0.5, 0.5 + plot_ticks, plot_ticks)))
    ax.set_yticks(list(arange(-0.5, 0.5 + plot_ticks, plot_ticks)))
    # add decors
    for g in solver.G: ax.add_artist(g.decor())
    # add new axis and color bar 
    ax = fg.add_axes([x + w + x, y, k, h])
    cb = fg.colorbar(pl, cax = ax)
    cb.set_ticks(list(arange(0, 1 + contour_ticks, contour_ticks)))
    # export figure to pdf
    if pdfdoc: pdfdoc.savefig(fg)
    # done
    return

####################################################### POTENTIAL CROSS SECTION

    def plotCrossSection(pdfdoc = None):
        """  """
        # some parameters
        # plot_ticks, contour_ticks = 0.1, 0.1

        # create figure
        fg = pp.figure()
        # set A4 paper dimensions
        fg.set_size_inches(8.2677, 11.6929)
        # create title
        fg.suptitle(f"Potential Cross Section")

        # create axis
        w, h, k = [*list(array([1, 1 / 1.4143])*0.7), 0.03]
        x, y = (1-w-k)/3, (1-h)/2
        ax = fg.add_axes([x, y, w, h]) 

        # x-coordinates
        d = self.l/(self.n-1)
        l = d*(self.n-1)/2
        X = arange(-l, l+d, d)

        # y-coordinates
        Y = self.D[1:-1, 1:-1] / MAX

        # make contour plot
        c = [pp.cm.inferno, pp.cm.terrain][0]
        args = {"cmap":c, "vmin":0.0, "vmax":1.0}
        pl = ax.contourf(X, X, Y, n, **args)

        # adjust ticks
        ax.set_xticks(list(arange(-0.5, 0.5 + plot_ticks, plot_ticks)))
        ax.set_yticks(list(arange(-0.5, 0.5 + plot_ticks, plot_ticks)))

        # add decors
        for g in self.G: ax.add_artist(g.decor())

        # add new axis and color bar 
        ax = fg.add_axes([x + w + x, y, k, h])
        cb = fg.colorbar(pl, cax = ax)
        cb.set_ticks(list(arange(0, 1 + contour_ticks, contour_ticks)))

        # export figure to pdf
        if pdfdoc: pdfdoc.savefig(fg)

        # done
        return

