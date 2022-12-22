# file: solvers.py
# author: Roch Schanen
# created: 

# from package: "https://numpy.org/"
from numpy import full, copy, shape, meshgrid
from numpy import load, save
from numpy import invert, gradient, arange

# from package: capycity (under heavy development)
from capycity.config import DATATYPE, ZER, MAX

##################################### Laplace equation solver in two dimensions

class laplace2DSolver():

    """ describe here the general use of this class """

    # p = None    # power of two
    # n = None    # number of intervals
    # l = None    # physical length
    # G = []      # mask geometries
    # T = []      # mask types

    def __init__(self, p = 3):
        """ For computational reason, the grid is square and its size
        is always a power of two. p is the power value and n the grid
        size. Both masks are the same size as the data array. An extra
        layer of zeroes is added at the edge of the square arrays for. """

        # local
        self.G, self.T = [], []
        self.dD = None, None

        # debug
        # print(f"Run Laplace 2D Solver")

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

        # compute mesh
        self._computemesh()

        # done
        return

    def _computemesh(self):
        d = self.l / self.n         # interval
        l = (self.l - d) / 2.0      # limit
        D = arange(-l, +l + d, d)   # domain
        self.mesh = meshgrid(D, D)  # mesh
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

        # done`
        return

    def mergeMask(self, maskGeometry, maskType = "1"):
        """ 1 is for SET mask,  0 is for CLEAR mask """

        # add new geometry
        self.G.append(maskGeometry)
        self.T.append(maskType)
        
        # merge according to type
        {   "0": self._mergeClrMask,
            "1": self._mergeSetMask,
        }[maskType](maskGeometry.mask(self))
        
        # done
        return

    # merge type C mask
    def _mergeClrMask(self, C):
        self.C &= invert(C)
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

        VERBOSE = False

        if VERBOSE: print("<-", end = "")

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

        if VERBOSE: print("quad-", end = "")

        # quadruplicate the data set
        self.D[1:-1:2, 1:-1:2] = D
        self.D[2:-1:2, 1:-1:2] = D
        self.D[1:-1:2, 2:-1:2] = D
        self.D[2:-1:2, 2:-1:2] = D

        if VERBOSE: print("mask-", end = "")

        # re-build masks
        self.C = full([self.n+2]*2, MAX, DATATYPE) # clear mask
        self.S = full([self.n+2]*2, ZER, DATATYPE) # set mask
        for g, t in zip(self.G, self.T):
            {   "0": self._mergeClrMask,
                "1": self._mergeSetMask,
            }[t](g.mask(self))

        if VERBOSE: print("mesh-", end = "")
        # re-compute mesh
        self._computemesh()

        if VERBOSE: print(">")

        # done
        return

    def save(self, fp):
        # save array without egdes
        save(fp, self.D[+1:-1, +1:-1])
        # done
        return

    def load(self, fp):
        # load the array
        D = load(fp)
        # retrieve size
        self.n, n = shape(D)
        # compute power of two by
        # counting binary shifts
        p = 0
        while n:
            n >>= 1
            p  += 1
        # record parameter
        self.p = p - 1
        # reserve memory
        self.D = full([self.n+2]*2, ZER, DATATYPE)
        """ the new data set is initialised to zeros.
        This automatically clears the edges. """
        # copy bulk
        self.D[+1:-1, +1:-1] = D
        # build masks
        self.C = full([self.n+2]*2, MAX, DATATYPE) # clear mask
        self.S = full([self.n+2]*2, ZER, DATATYPE) # set mask
        for g, t in zip(self.G, self.T):
            {   "0": self._mergeClrMask,
                "1": self._mergeSetMask,
            }[t](g.mask(self))
        # compute mesh
        self._computemesh()
        # done
        return

    def computegradient(self):
        # compute interval
        d = self.l / self.n
        # compute gradient without edges
        dDy, dDx = gradient(self.D[1:-1, 1:-1] / MAX)
        # fix units and record gradient without edges
        self.dD = dDx / d, dDy / d
        # done
        return
