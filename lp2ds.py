
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
line next to use 8 bits range data type """
# from numpy import uint8 as DATATYPE

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

    def jacobiSteps(self, n):

        # n jacobiStep iterations
        for i in range(n):
            self.jacobiStep()

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
        """ increase the resolution by a factor two. One can notice that
        with large resolution, the jacobi iteration "propagates" the potential
        "wave" at a slower "velocity". At lower resolutions the convergence 
        to a solution is quicker. In order to benefit from that fact, an  """

        # update new resolution
        self.p += 1
        self.n = 2**self.p

        # save current data for the upgrade
        D = copy(self.D[1:-1, 1:-1])
        """ ignore edges. """

        # reserve memory for the new data set
        self.D = full([self.n+2]*2, ZER, DATATYPE)
        """ the new data set is initialised to zeros.
        This is clearing up the edges. The rest of the
        array is defined during the "quadruplication """

        # quadruplicate the data set
        self.D[1:-1:2, 1:-1:2] = D
        self.D[2:-1:2, 1:-1:2] = D
        self.D[1:-1:2, 2:-1:2] = D
        self.D[2:-1:2, 2:-1:2] = D

        # re-build masks
        self.C = full([self.n+2]*2, MAX, DATATYPE) # clear mask
        self.S = full([self.n+2]*2, ZER, DATATYPE) # set mask
        for g, t in zip(self.G, self.T):
            {   "C": self._mergeClrMask,
                "S": self._mergeSetMask,
            }[t](g.mask(self.D))

        # done
        return

#####################################################################
###                            SOLVER                            ####
#####################################################################

from capycity.geometry import Disk
from capycity.graphics import plotPotentialContourFill

p = PdfPages("results.pdf")

# instanciate solver
l = laplace2DSolver(6) # 128 X 128

l.mergeMask(Disk(0.10), "S") # CENTER 0.2 diameter
l.mergeMask(Disk(0.40), "C") # SHELL 0.8 diameter 

# apply boundaries and display
l.applyMasks()
plotPotentialContourFill(l, n = 30, pdfdoc = p)

# run jacobi steps and display
l.jacobiSteps(700)
plotPotentialContourFill(l, n = 30, pdfdoc = p)

p.close()
