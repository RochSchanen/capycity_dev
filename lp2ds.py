
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
from numpy import save
from numpy import load

# define data type
from numpy import uint32 as DATATYPE
""" the precision does not influence the error when storing the data.
the idea that the result precision is 1/256 for the uint8 data type
turns out to be naive. Even when the number of iteration has defintely
passed the test of stability for a laplace solution. The final result
can still be largely influenced by the integer precision. This can
easily be checked empirically. As it turns out, a size of 32 bits seems
to be the minimum required to provide a reliable solution. (This is
higher than I personaly expected. I thought uint16 would be largely
sufficient but it is not...). For development purposes, using 8 bits
is a good ideas when checking operations on arrays. just uncomment the
line next for this purpose """
from numpy import uint8 as DATATYPE

# define ZER and MAX values
ZER, MAX = array([0, -1], DATATYPE)
""" ZER is always represented by a chain of "zero" bits.
MAX is always represented by a chain of "one" bits.
The size of the chain depends on the integer size. """

class laplace2DSolver():

    """ describe the genral use of this class """

    p = None    # power of two
    n = None    # number of intervals
    l = None    # physical length

    G = []      # mask geometries
    T = []      # mask types

    def __init__(self, p = 3):
        """ For computational reason, the grid is square and its size
        is always a power of two. p is the power value and n the grid
        size. Both masks are the same size as the data array. An extra
        layer of zeroes is added at the edge of the square arrays for. """

        # debug
        print(f"Run Laplace 2D Solver")

        # set resolution
        self.p = p
        self.n = 2**p

        # set default data (zeros)
        self.D = full([self.n+2]*2, ZER, DATATYPE) # data
        """ The intial set of data is now set a zero. As it turns out
        empirically, a set of random numbers might converge slightly
        faster towards the solution. However, the factor by which the
        convergence is increased is on the order of one. This is not
        a major increse in convergence speed. """

        # set default masks (transparent)
        self.C = full([self.n+2]*2, MAX, DATATYPE) # clear mask
        self.S = full([self.n+2]*2, ZER, DATATYPE) # set mask

        # set array physical length [m]
        self.l = +1.0E+0

        # done
        return

    def applyMasks(self):
        """ The resetting of the boundary conditions is
        performed by using two simple bitwise logical
        operations: "OR" and "AND" """

        # apply clear mask
        self.D &= self.C
        """ The "C" mask force data point to the value "ZER".
        The technique used is to apply an bitwise logical "AND"
        to the data array using the currently defined "clear" mask.
        Mask regions with the MAX value leave the data unchanged,
        while mask regions with the ZER value will force the data
        to "ZER". A logical AND instead of a product is prefered
        for speed. """

        # apply the set mask
        self.D |= self.S 
        """ The "S" mask force data point to the value "MAX".
        The technique used is to apply an bitwise logical "OR"
        to the data array using the currently defined "set" mask.
        Mask regions with the ZER value leave the data unchanged,
        while mask regions with the MAX value will force the data
        to "MAX". A logical OR instead of a product is prefered
        for speed. """
        
        # done
        return

    def jacobiStep(self):
        """ operations a performed directly on arrays
        which should reduces the time of execution.
        four array copies are fast operations. The
        logical shifts to the right are also fast.
        Unfortunately, they have to be performed 
        before the sums to prevent from any sum
        overflow. The additional memory used for
        storing the extra arrays should not be an
        issue. maybe we should lower the MAX value
        by 1/4 to associate the shift operations
        into a single one while preventing overflow.
        To do: investigate on that lamda factor to
        overshoot the Jacobi step and maybe converge
        substantially faster. """
        
        # compute shifted arrays
        l = copy(self.D[+1:-1, 0:-2]) # left
        r = copy(self.D[+1:-1, 2:  ]) # right
        t = copy(self.D[0:-2, +1:-1]) # top
        b = copy(self.D[2:  , +1:-1]) # bottom

        # compute average (Jacobi step)
        self.D[1:-1, 1:-1] = l//4 + r//4 + t//4 + b//4

        # reset boundaries
        self.applyMasks()

        # done
        return

    def mergeMask(self, maskGeometry, maskType = "S"):

        # add new geometry
        self.G.append(maskGeometry)
        self.T.append(maskType)
        
        # merge according to type
        {   "C": self._mergeClrMask,
            "S": self._mergeSetMask,
        }[maskType](maskGeometry.mask(self.D))
        
        # done
        return

    # merge type C mask
    def _mergeClrMask(self, C):
        self.C &= C
        return

    # merge type S mask
    def _mergeSetMask(self, S):
        self.S |= S
        return

    def incrementResolution(self):
        
        """ increase the resolution by a factor two. """

        # set new resolution
        self.p += 1
        self.n = 2**self.p

        # save current data for the upgrade
        D = copy(self.D[1:-1, 1:-1])
        """ edges are ingnored. """

        # reserve memory for the new data set
        self.D = full([self.n+2]*2, ZER, DATATYPE)
        """ the new data set is initialised to ZER. This is used to
        clear the data set edges. The rest of the array is overidden
        during the "quadruplication" performed next. """

        # quadruplicate the data set
        self.D[1:-1:2, 1:-1:2] = D
        self.D[2:-1:2, 1:-1:2] = D
        self.D[1:-1:2, 2:-1:2] = D
        self.D[2:-1:2, 2:-1:2] = D
        """ the initial increase """

        # reset default masks (transparent)
        self.C = full([self.n+2]*2, MAX, DATATYPE) # clear mask
        self.S = full([self.n+2]*2, ZER, DATATYPE) # set mask

        # merge masks
        for g, t in zip(self.G, self.T):
            # merge according to type
            {   "C": self._mergeClrMask,
                "S": self._mergeSetMask,
            }[t](g.mask(self.D))

        # done
        return

    def plotPotentialDistribution(self, n = 30, pdfdoc = None):
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

#####################################################################
###                            SOLVER                            ####
#####################################################################

from capycity.geometry import Disk

SETUP, SOLVE, DISPLAY = 1, 1, 0

# instanciate solver
l = laplace2DSolver(3)

if SETUP:

    l.mergeMask(Disk(0.20), "S") # MAX
    l.mergeMask(Disk(0.40), "C") # ZER

if SOLVE:

    # compute
    print(l.D)
    print()

    l.jacobiStep()
    print(l.D)
    print()

    l.incrementResolution()
    print(l.D)
    print()

    for i in range(23):
        l.jacobiStep()
    print(l.D)
    print()

    # update data after solving
    # save(f"potential.npy", l.D)

if DISPLAY:

    l.D = load(f"potential.npy")

    # export graphs
    p = PdfPages("results.pdf")
    l.plotPotentialDistribution(n = 50, pdfdoc = p)
    p.close()
