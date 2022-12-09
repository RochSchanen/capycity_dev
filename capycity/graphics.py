# file: graphics.py

# from package: "https://matplotlib.org/"
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import contourf
from matplotlib.pyplot import sca
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

def plotcontour(solver, *args, name = "contour", **kwargs):

    # some parameters
    ptk, ctk = 0.1, 0.1

    # create figure
    fg = pp.figure(name)
    
    # set A4 paper dimensions
    fg.set_size_inches(8.2677, 11.6929)
    
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
    pl = ax.contourf(X, X, Y, *args, **mapargs, **kwargs)
    
    # adjust ticks
    ax.set_xticks(list(arange(-0.5, 0.5 + ptk, ptk)))
    ax.set_yticks(list(arange(-0.5, 0.5 + ptk, ptk)))
    
    # add decors
    for g in solver.G:
        ax.add_artist(g.contour2D_decor())
    
    # add new axis and color bar 
    ax = fg.add_axes([x + w + x, y, k, h])
    cb = fg.colorbar(pl, cax = ax)
    cb.set_ticks(list(arange(0, 1 + ctk, ctk)))

    # set contour as the current axis
    sca(fg.get_axes()[0])

    # done
    return

###################################################### POTENTIAL IMAGE MESH

def plotmesh(solver, *args, name = "mesh", **kwargs):

    # some parameters
    ptk, ctk = 0.1, 0.1

    # create figure
    fg = pp.figure(name)
    
    # set A4 paper dimensions
    fg.set_size_inches(8.2677, 11.6929)
    
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

    from matplotlib.colors import BoundaryNorm
    from matplotlib.ticker import MaxNLocator
    levels = MaxNLocator(nbins=15).tick_values(0.0, 1.0)
    norm = BoundaryNorm(levels, ncolors=c.N, clip=True)
    mapargs = {
        "cmap":c,
        # "vmin":0.0,
        # "vmax":1.0,
        "norm":norm,
        "shading":'auto',
        }

    pl = ax.pcolormesh(X, X, Y, *args, **mapargs, **kwargs)
    
    # adjust ticks
    ax.set_xticks(list(arange(-0.5, 0.5 + ptk, ptk)))
    ax.set_yticks(list(arange(-0.5, 0.5 + ptk, ptk)))
    
    # add decors
    for g in solver.G:
        ax.add_artist(g.contour2D_decor())
    
    # add new axis and color bar 
    ax = fg.add_axes([x + w + x, y, k, h])
    cb = fg.colorbar(pl, cax = ax)
    cb.set_ticks(list(arange(0, 1 + ctk, ctk)))

    # set contour as the current axis
    sca(fg.get_axes()[0])

    # done
    return

# FIGURE FUNCTIONS FOR EXPORT

def selectfigure(name):         
    if not pp.fignum_exists(name):
        # create figure
        fg = pp.figure(name)
        # set A4 paper dimensions
        fg.set_size_inches(8.2677, 11.6929)
        # create square axis
        w, h = array([1, 1 / 1.4143])*0.7
        x, y = (1-w)/2, (1-h)/2
        ax = fg.add_axes([x, y, w, h])
    else:
        # select figure
        fg = pp.figure(name)
        # select axes
        ax = fg.get_axes()[0]
    # done
    return fg, ax

# PDF DOCUMENT FUNCTION FOR EXPORT

_DOC = None

def opendocument(pathname):
    global _DOC
    _DOC = PdfPages(pathname)
    return _DOC

def exportfigure(name):
    fg, ax = selectfigure(name)
    _DOC.savefig(fg)
    return

def closedocument():
    _DOC.close()
    return
