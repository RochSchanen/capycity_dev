
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
to be the minimum requirement to provide a reliable solution. (This is
higher than I personaly expected. I thought uint16 would be largely
sufficient but it is not...) """

# define ZER and MAX values
ZER, MAX = array([0, -1], DATATYPE)
""" ZER is always represented by a chain of "zero" bits.
MAX is always represented by a chain of "one" bits.
The size of the chain depends on the integer size. """

class laplace2DSolver():

    """ """    

    l = None    # length
    p = None    # power of two
    n = None    # intervals

    def __init__(self, p = 3):
        """For computational reason, the grid is square and its size
        is always a power of two. p is the power value and n the grid
        size. Both masks are the same size as the data array. An extra
        layer of zeroes is added at the edge of the square arrays for
        ."""

        # debug
        print(f"Run Laplace 2D Solver")

        # set resolution
        self.p = p
        self.n = 2**p

        """ The intial set of data is now set a zero. As it turns out
        empirically, a set of random numbers might converge slightly
        faster towards the solution. However, the factor by which the
        convergence is increased is on the order of one. This is not
        a major increse in convergence speed. """
        # set default data (zeros)
        self.D = full([self.n+2]*2, ZER, DATATYPE) # data

        # set default masks (transparent)
        self.C = full([self.n+2]*2, MAX, DATATYPE) # mask clear 
        self.S = full([self.n+2]*2, ZER, DATATYPE) # mask set

        # set full size
        self.l = +1.0E+0  # grid length [m]

        # done
        return


    def applyMasks(self):
        """ the re-set of the boundary conditions is reduced to two
        elementary bitwise logical operations "OR" and "AND" that
        should be very fast to execute on our integer arrays """
        # apply the clear mask
        self.D &= self.C 
        # apply the set mask
        self.D |= self.S 
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
        Also, investigate on that lamda factor to
        overshoot the Jacobi step and maybe converge
        substantially faster. """
        
        # compute shifted array
        l = copy(self.D[+1:-1, 0:-2]) # left
        r = copy(self.D[+1:-1, 2:  ]) # right
        t = copy(self.D[0:-2, +1:-1]) # top
        b = copy(self.D[2:  , +1:-1]) # bottom

        # compute average (Jacobi step)
        self.D[1:-1, 1:-1] = l//4 + r//4 + t//4 + b//4

        # re-set boundaries
        self.applyMasks()

        # done
        return

    # todo: use a flag to determine the type of merging:
    # merge onto the "CLEAR" mask, or onto the "SET" mask
    def mergeCMask(self, C):
        """ """
        self.C &= C # clear mask
        return

    def mergeSMask(self, S):
        """ """
        self.S |= S # clear mask
        return

    def plotPotentialContourf(self, n = 30, pdfdoc = None):
        """ """
        # some parameters
        plot_ticks, contour_ticks = 0.1, 0.1
        # create figure
        fg = pp.figure()
        # set A4 paper dimensions
        fg.set_size_inches(8.2677, 11.6929)
        # create title
        fg.suptitle(f"Distribution")
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
        # if self.decors:
        #     for r in self.decors:
        #         ax.add_artist(self.getDiskDecor(r))
        # add new axis and color bar 
        ax = fg.add_axes([x + w + x, y, k, h])
        cb = fg.colorbar(pl, cax = ax)
        cb.set_ticks(list(arange(0, 1 + contour_ticks, contour_ticks)))
        # export figure to pdf
        if pdfdoc: pdfdoc.savefig(fg)
        # done
        return

from capycity.geometry import Disk

SOLVE, DISPLAY = 1, 1

# instanciate solver
l = laplace2DSolver(7)

if SOLVE:

    # create disk
    D1 = Disk(0.15)
    M = D1.getMask(l.D)
    l.mergeSMask(M)

    # create shell
    D2 = Disk(0.45)
    M = D2.getMask(l.D)
    l.mergeCMask(M)

    # compute
    for i in range(1000):
          l.jacobiStep()

    # update data after solving
    save(f"potential.npy", l.D)

if DISPLAY:

    l.D = load(f"potential.npy")

    # export graphs
    p = PdfPages("results.pdf")
    l.plotPotentialContourf(n = 20, pdfdoc = p)
    p.close()
