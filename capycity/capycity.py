
# VERSION
version = 0.01
"""
    - remove the powers of two constrain.
    - keep the square shape

    improve:

        - use diff instead of gradient to prevent usage of float
        and keep the use of integers to accelerate calculations and
        prevent from using too much memory (see also multi grid).

        - decrease the maximum value by a factor 4. Sum the arrays
            first and then divide (double shift) the sum by fouronly
            once instead of four times as it is implemented now.

    next steps:

        - implement multi grid:
            keep the square shape. use concentric increase of the
            resolution with a list of concentric sizes: it is a
            re-implementation of the JacobiStepSeries where the size
            of the higher resolution map is redefined as a subset of
            the original map. The outer coarser map works where the
            potential gradients are small (towards the edges with a
            shield of large size).
        
        - allow rectangular shape
"""

# LIBRARIES

# from package: "https://numpy.org/"
from numpy import linspace, logspace
from numpy import ceil, floor
from numpy import log, log10, exp
from numpy import full
from numpy import copy
from numpy import diff
from numpy import array
from numpy import invert
from numpy import arange
from numpy import asarray
from numpy import gradient
from numpy import meshgrid
from numpy import sum as SUM
from numpy import multiply as MLT

# from package "https://python-pillow.org/"
from PIL import Image
from PIL import ImageDraw

# from package "https://scipy.org/"
from scipy.constants import epsilon_0 as EPS0

# from package: "https://matplotlib.org/"
from matplotlib.pyplot import sca
from matplotlib.pyplot import Circle
from matplotlib.pyplot import figure
from matplotlib.pyplot import Rectangle
from matplotlib.pyplot import fignum_exists
from matplotlib.pyplot import cm
from matplotlib.backends.backend_pdf import PdfPages

# DATA TYPE

""" The data type precision is fixed here. Empirically, the algorithm
converges well when using 32 bits unsigned integers. A lower 8 bits
unsigned integer type is not sufficient but can be used for
debugging. Select the required length of bits by un-commenting one of
following three lines: """

# from numpy import uint8  as _DT
# from numpy import uint16  as _DT
from numpy import uint32 as _DT

# EXTREMA (zero value and maximum value)
ZV, MV = array([0, -1], _DT)

# DEBUG DISPLAY
VERBOSE_MAP = False
VERBOSE = False

class SolverTwoDimensions():

    def __init__(self, nx = 8, ny = None, lx = 1.0, ly = None):

        if VERBOSE: print(f"instantiate 'SolverTwoDimensions'")

        # GRID
        
        """ the grid size is set for all the maps and masks.
        When the grid resolution is changed, maps and masks
        are upgraded to fit the new grid size. """
        
        # compute initial size
        if ny is None: ny = nx
        if VERBOSE: print(f"nx = {nx}, ny = {ny}")

        # LENGTH

        """ the length is set for a square shape only so far.
        therefore, just use square grids for the moment. it will
        be extended to two independent lengths later. The issue
        is that the grid size and length ratio must be kept
        identical """

        if ly is None: ly = ny*lx/nx
        if VERBOSE: print(f"lx = {lx}, ly = {ly}")

        # local
        self.ll = lx, ly    # grid length
        self.nn = nx, ny    # grid size

        self._computemesh()

        # MAPS AND MASKS

        """ Masks define the geometry of the conducting materials
        that compose the the capacitance to calculate. There is
        one potential field map to compute per conductor/mask. New
        maps are automatically created when new masks are added.
        Masks are build using pre-defined or user-defined methods.
        A single mask can be obtained by merging several mask."""
        
        self.M = {}

        # current step number
        self.k = 1

        # done
        return

    def addMap(self, name, mask = None):

        """ create a new map to hold a new component to the
        capacitance. A CLEAR map and a ANCHOR map are automatically
        added for the computation. The ANCHOR map should contain
        the shape of the solid. Where the maps is solid, the value is
        VM, otherwise, the value is ZV. A CLEAR map is build from all
        the other maps already available. """
        
        if VERBOSE: print(f"new map '{name.upper()}'")
        # instantiate new map class
        self.M[name] = _map(self.nn)
        # build CLEAR mask
        for k in self.M.keys():
            if k == name: continue
            # merge clear map 
            self.M[name].C &= invert(self.M[k].A)
        # no mask, done
        if mask is None: return
        # add mask
        self.addMask(name, mask)
        # done
        return

    def addMask(self, name, mask):
        # merge new mask to the "named map"
        self.M[name].setMask(mask, self.nn, self.ll)
        # update clear maps
        I = invert(self.M[name].A)
        for k in self.M.keys():
            if k == name: continue
            # merge clear map
            self.M[k].C &= I
        # done
        return

    def _computemesh(self):
        nx, ny = self.nn
        lx, ly = self.ll
        # intervals
        dx, dy = lx / nx, ly / ny
        # boundaries
        bx, by = (lx - dx) / 2.0, (ly - dy) / 2.0
        # domains
        Dx, Dy = linspace(-bx, +bx, nx), linspace(-by, +by, ny)
        # mesh
        self.mesh = meshgrid(Dx, Dy)
        # done
        return

    def jacobiSteps(self, n, C1, C2, SavePattern = None):
        """ using the same number of steps on every map is justified
        if the convergence rate depends mostly on the size/resolution
        of the grid and not so much on the geometry of the system.
        This should also allow to devise an alternating series of
        incrementing the resolution followed by a steps series that
        could be used for all cases: see the next method """
        K, C = [], []
        for i in range(n):
            # n Jacobi step iterations over C1 and C2
            self.M[C1].jacobiStep()
            self.M[C2].jacobiStep()
            # save if step match pattern element
            if (self.k+i) in SavePattern:
                K.append(self.k+i)
                C.append(self.computeCapacitance(C1, C2))
        # increment step counter
        self.k += n
        # return a single value
        if SavePattern is None:
            return self.computeCapacitance(C1, C2)
        # return saved data
        return K, C

    def jacobiStepSeries(self, S, C1, C2, n = 500):
        """ S is a list of step series numbers. After each step
        series, the resolution of the grid is incremented by a factor
        two. This is why the resolution of the grid is always a power
        of two. C1 and C2 are the names of the two conductors for
        which we compute the mutual capacitance. The other conductors
        affect the computation only by their presence (geometry) and
        the value of their potential is kept at zero (See theoretical
        considerations from the report). At the end of each series,
        the total number of steps and the value of the capacitance is
        computed and recorded in list K and R. This allows to monitor
        the convergence the series. """
        if VERBOSE:
            nx, ny = self.nn
            print(f"resolution: {nx}x{ny}")
        # compute saving pattern (logarithmic spacing)
        decade = ceil(log(sum(S)) / log(10))
        P = list(logspace(0, decade, n, dtype = 'int'))
        # make sure the last step is recorded
        P.append(sum(S))
        # declare lists
        K = []  # steps
        C = []  # capacitances
        # loop through series
        for n in S[:-1]:
            # compute steps series
            k, c = self.jacobiSteps(n, C1, C2, P)
            # append data points
            K += k; C += c
            # increment resolution
            self.incrementResolution()
            if VERBOSE:
                nx, ny = self.nn
                print(f"resolution: {nx}x{ny}")
            # add first point after resolution increment
            K += [self.k]
            C += [self.computeCapacitance(C1, C2)]
            self.k += 1
        # compute last steps series
        k, c = self.jacobiSteps(S[-1], C1, C2, P)
        # append last data points
        K += k; C += c
        # done 
        return K, C

    def computegradient(self, name):
        """ The electric field is directly computed
        as the gradient of the electric scalar potential """
        nx, ny = self.nn
        lx, ly = self.ll
        # compute intervals
        dx, dy = lx / nx, ly / ny
        # compute gradient without edges
        dDy, dDx = gradient(self.M[name].D[1:-1, 1:-1] / MV)
        # fix units and record gradient without edges
        self.M[name].dD = dDx / dx, dDy / dy
        # done
        return

    def computeCapacitance(self, name1, name2):
        nx, ny = self.nn
        lx, ly = self.ll
        # compute interval (intervals should be identical)
        dx, dy = lx / nx, ly / ny
        # compute C1 map gradient (electric field)
        self.computegradient(name1)
        # get component fields
        EX1, EY1 = self.M[name1].dD
        # compute C2 map gradient (electric field)
        self.computegradient(name2)
        # get component fields
        EX2, EY2 = self.M[name2].dD
        # compute 'field dot products'
        integrand = MLT(EX1, EX2) + MLT(EY1, EY2)
        # compute numeric integral
        integral = SUM(integrand)*dx*dy
        # normalise units and return the mutual
        # capacitance as a positive value
        c = -integral*EPS0
        # done
        return c

    def incrementResolution(self):

        """ increase the resolution by a factor two. It can be noticed
        empirically that with larger resolution, the Jacobi iteration
        "propagates" the scalar potential "wave" at a slower "velocity".
        At lower resolutions the convergence to an approximate solution
        is quicker. In order to benefit from that observation, a gradual
        incrementation of the resolution is performed. The convergence
        to an accurate solution is preserved without compromising to much
        the computing time. """

        if VERBOSE: print("<-", end = "")

        # get current values
        nx, ny = self.nn        

        # update local
        nx, ny = nx<<1, ny<<1

        # update values
        self.nn = nx, ny

        # debug
        if VERBOSE: print("map.D-", end = "")

        # loop through maps
        for k in self.M.keys():

            # save current data for the upgrade
            D = copy(self.M[k].D[1:-1, 1:-1])
            """ edges are ignored. """

            # reserve memory for the new data set
            self.M[k].D = full([nx+2, ny+2], ZV, _DT)
            """ the new data set is initialised to zeros.
            This is clearing up the edges. The rest of the
            array is defined during the 'quadruplication' """

            if VERBOSE: print("quad-", end = "")
            # quadruplicate the data set
            self.M[k].D[1:-1:2, 1:-1:2] = D
            self.M[k].D[2:-1:2, 1:-1:2] = D
            self.M[k].D[1:-1:2, 2:-1:2] = D
            self.M[k].D[2:-1:2, 2:-1:2] = D

        # debug
        if VERBOSE: print("map.A-", end = "")

        # loop through maps
        for k in self.M.keys():
            # re-build ANCHOR mask
            self.M[k].refreshMask(self.nn, self.ll)

        # debug
        if VERBOSE: print("map.C-", end = "")

        # loop through maps
        for k in self.M.keys():
            # reserve memory for the CLEAR mask
            self.M[k].C = full([nx+2, ny+2], MV, _DT)
            # loop though all the other maps
            for l in self.M.keys():
                if k == l: continue
                # merge masks
                self.M[k].C &= invert(self.M[l].A)

        # debug
        if VERBOSE: print("mesh-", end = "")

        # re-compute mesh
        self._computemesh()

        # debug
        if VERBOSE: print(">")

        # done
        return

class _map():

    def __init__(self, nn):
        if VERBOSE: print(f"instantiate '_map'")
        # get grid size
        nx, ny = nn
        if VERBOSE: print(f"nx = {nx}, ny = {ny}")
        # new data map (potential between conductors)
        self.D = full([nx+2, ny+2], ZV, _DT)
        # clear map (zero potential from other conductors)
        self.C = full([nx+2, ny+2], MV, _DT)
        # anchor map (unit potential from this conductor)
        self.A = full([nx+2, ny+2], ZV, _DT)
        # debug
        if VERBOSE_MAP: print(f"map:")
        if VERBOSE_MAP: print(f"{self.D}")
        # declare masks list
        self.M = None
        # declare gradients list
        self.dD = None
        # done
        return

    """ maybe the two next
    methods can be merged...
    or made independent... """

    """ This adds a new solid shape
    that is registered by merging it
    to the existing ANCHOR map. There
    is only one mask implemented for
    the moment.
    """
    def setMask(self, mask, nn, ll):
        # register the class instance
        self.M = mask
        # and setup anchor map
        self.refreshMask(nn, ll)
        # done
        return

    """ when there is more than one
    function to build the mask, a set
    of masks that must be merged together.
    This is done through refreshMask()
    """
    def refreshMask(self, nn, ll):
        
        # go through all map building function HERE!
        # go through all map building function HERE!
        # go through all map building function HERE!

        # use the registered mask methods
        # self.A |= self.M.mask(nn, ll)
        self.A = self.M.mask(nn, ll)
        # debug        
        if VERBOSE_MAP: print(f"anchor map:")
        if VERBOSE_MAP: print(f"{self.A}:")
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
        to "ZER". A logical AND instead of a product is preferred
        for speed. """

        # apply the set mask
        self.D |= self.A
        """ The "S" mask force data point to the value "MAX".
        The technique used is to apply an bitwise logical "OR"
        to the data array using the currently defined "set" mask.
        Mask regions with the ZER value leave the data unchanged,
        while mask regions with the MAX value will force the data
        to "MAX". A logical OR instead of a product is preferred
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
        To do: investigate on that lambda factor to
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

class DiskSolid():

    # parameters
    def __init__(self, r):

        if VERBOSE: print("instanciate 'DiskSolid'")

        # record parameters
        self.l, self.r = -r, +r
        self.t, self.b = +r, -r
        if VERBOSE: print(f"l, r, t, b = ", end = '')
        if VERBOSE: print(f"{self.l}, {self.r}, {self.t}, {self.b}")

        # done
        return

    # mask
    def mask(self, nn, ll):
        if VERBOSE: print(f"create disk solid mask")
        # get array size (without edges)
        nx, ny = nn
        # get array length
        lx, ly = ll
        # compute index origin
        ox = (nx + 1.0) / 2.0
        oy = (ny + 1.0) / 2.0
        # compute interval
        ax = lx / nx
        ay = ly / ny
        if VERBOSE: print(f"nx ny lx ly ox oy = ", end ='')
        if VERBOSE: print(nx, ny, lx, ly, ox, oy)
        # boundary box
        l, r =  ceil(self.l / ax + ox), floor(self.r / ax + ox)
        t, b = floor(self.t / ay + oy),  ceil(self.b / ay + oy)
        if VERBOSE: print(f"l, r, t, b = ", end = '')
        if VERBOSE: print(f"{l}, {r}, {t}, {b}")
        # create image
        i = Image.new(mode = "L", size = (nx+2, ny+2), color = 0)
        # link handle "d" to image "i"
        d = ImageDraw.Draw(i)
        # punch aperture
        d.ellipse([l, b, r, t], fill = 1, width = 0)
        # convert to numpy array
        M = asarray(i, _DT) * MV
        # done
        return M

    # graph for 2D display
    def decor(self):
        c = Circle(
            (0.0, 0.0),                     # origin
            self.r,                         # size
            edgecolor = (0.7, 0.7, 0.7),    # colour
            linestyle = "--",               # line style
            fill = False,                   # fill style
            )
        # done
        return c

class DiskAperture():

    # parameters
    def __init__(self, r):

        if VERBOSE: print("instantiate 'DiskAperture'")

        # record parameters
        self.l, self.r = -r, +r
        self.t, self.b = +r, -r
        if VERBOSE: print(f"l, r, t, b = ", end = '')
        if VERBOSE: print(f"{self.l}, {self.r}, {self.t}, {self.b}")

        # done
        return

    # mask
    def mask(self, nn, ll):
        if VERBOSE: print(f"create disk aperture mask")
        # get array size (without edges)
        nx, ny = nn
        # get array length
        lx, ly = ll
        # compute index origin
        ox = (nx + 1.0) / 2.0
        oy = (ny + 1.0) / 2.0
        # compute interval
        ax = lx / nx
        ay = ly / ny
        if VERBOSE: print(f"nx ny lx ly ox oy = ", end ='')
        if VERBOSE: print(nx, ny, lx, ly, ox, oy)
        # boundary box
        l, r =  ceil(self.l / ax + ox), floor(self.r / ax + ox)
        t, b = floor(self.t / ay + oy),  ceil(self.b / ay + oy)
        if VERBOSE: print(f"l, r, t, b = ", end = '')
        if VERBOSE: print(f"{l}, {r}, {t}, {b}")
        # create image
        i = Image.new(mode = "L", size = (nx+2, ny+2), color = 1)
        # link handle "d" to image "i"
        d = ImageDraw.Draw(i)
        # punch aperture
        d.ellipse([l, b, r, t], fill = 0, width = 0)
        # convert to numpy array
        M = asarray(i, _DT) * MV
        # done
        return M

    # graph for 2D display
    def decor(self):
        c = Circle(
            (0.0, 0.0),                     # origin
            self.r,                         # size
            edgecolor = (0.7, 0.7, 0.7),    # colour
            linestyle = "--",               # line style
            fill = False,                   # fill style
            )
        # done
        return c

class PlateSolid():

    # parameters
    def __init__(self, x, y, w, h):

        if VERBOSE: print("instantiate 'PlateSolid'")

        # record parameters
        self.l, self.r = x - w / 2.0, x + w / 2.0
        self.t, self.b = y + h / 2.0, y - h / 2.0
        if VERBOSE: print(f"l, r, t, b = ", end = '')
        if VERBOSE: print(f"{self.l}, {self.r}, {self.t}, {self.b}")

        # done
        return

    # mask
    def mask(self, nn, ll):
        if VERBOSE: print(f"create plate solid mask")
        # get array size (without edges)
        nx, ny = nn
        # get array length
        lx, ly = ll
        # compute index origin
        ox = (nx + 1.0) / 2.0
        oy = (ny + 1.0) / 2.0
        # compute interval
        ax = lx / nx
        ay = ly / ny
        if VERBOSE: print(f"nx ny lx ly ox oy = ", end ='')
        if VERBOSE: print(nx, ny, lx, ly, ox, oy)
        # boundary box
        l, r =  ceil(self.l / ax + ox), floor(self.r / ax + ox)
        t, b = floor(self.t / ay + oy),  ceil(self.b / ay + oy)
        if VERBOSE: print(f"l, r, t, b = ", end = '')
        if VERBOSE: print(f"{l}, {r}, {t}, {b}")
        # create image
        i = Image.new(mode = "L", size = (nx+2, ny+2), color = 0)
        # link handle "d" to image "i"
        d = ImageDraw.Draw(i)
        # punch aperture
        d.rectangle([l, b, r, t], fill = 1, width = 0)
        # convert to numpy array
        M = asarray(i, _DT) * MV
        # done
        return M

    def decor(self):
        # get geometry
        l, r = self.l, self.r
        t, b = self.t, self.b
        # patch
        R = Rectangle((l, b), r-l, t-b,
            edgecolor = (0.7, 0.7, 0.7),    # colour
            linestyle = "--",               # line style
            fill        = False,            # no filling
            )
        # done
        return R

class PlateAperture():

    # parameters
    def __init__(self, x, y, w, h):

        if VERBOSE: print("instantiate 'PlateAperture'")

        # record parameters
        self.l, self.r = x - w / 2.0, x + w / 2.0
        self.t, self.b = y + h / 2.0, y - h / 2.0
        if VERBOSE: print(f"l, r, t, b = ", end = '')
        if VERBOSE: print(f"{self.l}, {self.r}, {self.t}, {self.b}")

        # done
        return

    # mask
    def mask(self, nn, ll):
        if VERBOSE: print(f"create plate aperture mask")
        # get array size (without edges)
        nx, ny = nn
        # get array length
        lx, ly = ll
        # compute index origin
        ox = (nx + 1.0) / 2.0
        oy = (ny + 1.0) / 2.0
        # compute interval
        ax = lx / nx
        ay = ly / ny
        if VERBOSE: print(f"nx ny lx ly ox oy = ", end ='')
        if VERBOSE: print(nx, ny, lx, ly, ox, oy)
        # boundary box
        l, r =  ceil(self.l / ax + ox), floor(self.r / ax + ox)
        t, b = floor(self.t / ay + oy),  ceil(self.b / ay + oy)
        if VERBOSE: print(f"l, r, t, b = ", end = '')
        if VERBOSE: print(f"{l}, {r}, {t}, {b}")
        # create image
        i = Image.new(mode = "L", size = (nx+2, ny+2), color = 1)
        # link handle "d" to image "i"
        d = ImageDraw.Draw(i)
        # punch aperture
        d.rectangle([l, b, r, t], fill = 0, width = 0)
        # convert to numpy array
        M = asarray(i, _DT) * MV
        # done
        return M

    def decor(self):
        # get geometry
        l, r = self.l, self.r
        t, b = self.t, self.b
        # patch
        R = Rectangle((l, b), r-l, t-b,
            edgecolor = (0.7, 0.7, 0.7),    # colour
            linestyle = "--",               # line style
            fill        = False,            # no filling
            )
        # done
        return R

class Document():

    def __init__(self, pathname = None):
        if pathname is not None:
            self._DOC = self.opendocument(pathname)
        return

    def opendocument(self, pathname):
        self._DOC = PdfPages(pathname)
        return self._DOC

    def exportfigure(self, name):
        if VERBOSE: print(f"export '{name}'")
        args = selectfigure(name)
        self._DOC.savefig(args[0])
        return

    def closedocument(self):
        self._DOC.close()
        return

def _getTickPositions(start, stop, ticks):

    def _getTickIntervals(start, stop, ticks):

        ln10 = 2.3025850929940459

        # trial table
        T = [0.010, 0.020, 0.025, 0.050,
             0.100, 0.200, 0.250, 0.500,
             1.000, 2.000, 2.500, 5.000]

        # corresponding tick sub division intervals
        S = [5.0,   4.0,   5.0,   5.0,
             5.0,   4.0,   5.0,   5.0,
             5.0,   4.0,   5.0,   5.0]

        span = stop - start                         # get span
        d = exp(ln10 * floor(log10(span)))          # find decade
        span /= d                                   # re-scale

        # find number of ticks below and closest to n
        i, m = 0, floor(span / T[0])                # start up
        while m > ticks:                            # next try?
            i, m = i + 1, floor(span / T[i + 1])    # try again 

        # re-scale
        mi =  d * T[i]   # main tick intervals
        si = mi / S[i]   # sub tick intervals

        # done
        return mi, si

    # get intervals
    mi, si = _getTickIntervals(start, stop, ticks)

    # main ticks (round is the built-in python version)
    ns = ceil(start / mi - 0.001) * mi  # start
    ne = floor(stop / mi + 0.001) * mi  # end
    p  = round((ne - ns) / mi) + 1      # fail safe
    M  = linspace(ns, ne, p)            # main positions

    # sub ticks (round is the built-in python version)
    ns = ceil(start / si + 0.001) * si  # start
    ne = floor(stop / si - 0.001) * si  # end
    p  = round((ne - ns) / si) + 1      # fail safe
    S  = linspace(ns, ne, p)            # sub positions

    # done
    return M, S

def selectfigure(name):
    if not fignum_exists(name):
        # create figure
        fg = figure(name)
        # set A4 paper dimensions
        fg.set_size_inches(8.2677, 11.6929)
        # create square axis
        w, h = array([1, 1 / 1.4143])*0.7
        x, y = (1-w)/2, (1-h)/2
        ax = fg.add_axes([x, y, w, h])
    else:
        # select figure
        # (here the figure can be of any type)
        fg = figure(name)
        # get axes
        ax = fg.get_axes()[0]
    # done
    return fg, ax

def selectmapfigure(name):
    # map figure is the same as selectfigure
    # but adds an extra axis for a side colour bar
    # it also returns an extra axis
    if not fignum_exists(name): 
        # create figure
        fg = figure(name)
        # set A4 paper dimensions
        fg.set_size_inches(8.2677, 11.6929)
        # create square axis
        w, h, k = [*list(array([1, 1 / 1.4143])*0.7), 0.03]
        x, y = (1-w-k)/3, (1-h)/2
        ax = fg.add_axes([x, y, w, h]) 
        # plus side axis for the colour bar 
        bx = fg.add_axes([x+w+x, y, k, h])
    else:
        # select an existing figure
        # (here the figure has to be a map figure)
        fg = figure(name)
        # get axes
        ax, bx = fg.get_axes()[:2]
    # done
    return fg, ax, bx

def showResolution(solver, ax):
    # get the mesh coordinates
    X, Y = solver.mesh
    # show resolution visually
    nx, ny = solver.nn
    lx, ly = solver.ll
    dx, dy = lx / nx, ly / ny
    # for i, j in [(0, 0), (0, -1), (-1, 0), (-1, -1)]:
    for i, j in [(1, 1), (1, -2), (-2, 1), (-2, -2)]:
        x, y = X[i, j], Y[i, j]
        l, r, t, b = x-dx/2, x+dx/2, y+dy/2, y-dy/2
        ax.add_patch(Rectangle((l, b), dx, dy,
            edgecolor   = "gainsboro",
            facecolor   = "gainsboro",
            zorder      = 3,
            ))
    # done
    return

def vplot(solver, figname, mapname, *args, **kwargs):
    # create/select a map figure
    fg, ax, bx = selectmapfigure(figname)
    # get the mesh coordinates
    X, Y = solver.mesh
    # get normalised data without the edges
    D = solver.M[mapname].D[1:-1, 1:-1] / MV
    # create filled contour map
    QCS = ax.contourf(X, Y, D, 32, *args, **kwargs)
    # adjust ticks
    lx, ly = solver.ll
    MX, SX = _getTickPositions(-solver.lx/2.0, +solver.lx/2.0, 9)
    MY, SY = _getTickPositions(-solver.ly/2.0, +solver.ly/2.0, 9)
    ax.set_xticks(MX)
    ax.set_yticks(MY)
    # create the colour bar
    cb = fg.colorbar(QCS, cax = bx)
    cb.set_ticks(list(arange(0.0, 1.1, 0.1)))
    # set colour map
    QCS.set_cmap(cm.inferno)
    # add decors
    for k in solver.M.keys():
        ax.add_artist(solver.M[k].M.decor())
    showResolution(solver, ax)
    # done
    return fg, ax, bx

def mplot(solver, figname, mapname):
    # create/select a map figure
    fg, ax, bx = selectmapfigure(figname)
    # get the mesh coordinates
    X, Y = solver.mesh
    # get normalised data without the edges
    D = solver.M[mapname].D[1:-1, 1:-1] / MV
    # create filled contour map
    QCS = ax.pcolormesh(X, Y, D, shading = 'auto', rasterized = True)
    # adjust ticks
    lx, ly = solver.ll
    MX, SX = _getTickPositions(-solver.lx/2.0, +solver.lx/2.0, 9)
    MY, SY = _getTickPositions(-solver.ly/2.0, +solver.ly/2.0, 9)
    ax.set_xticks(MX)
    ax.set_yticks(MY)
    # create the colour bar
    cb = fg.colorbar(QCS, cax = bx)
    cb.set_ticks(list(arange(0.0, 1.1, 0.1)))
    # set colour map
    QCS.set_cmap(cm.inferno)
    # add decors
    for k in solver.M.keys():
        ax.add_artist(solver.M[k].M.decor())
    showResolution(solver, ax)
    # done
    return fg, ax, bx

def headerText(text, fg):
    w, h = array([1, 1 / 1.4143])*0.7
    x, y = (1-w)/2, (1-h)/2
    tx = fg.text(x+w/2, 3*y/2+h, text)
    tx.set_fontfamily('monospace')
    tx.set_horizontalalignment('center')
    tx.set_verticalalignment('center')
    tx.set_fontsize('large')
    return tx

def footerText(text, fg):
    w, h = array([1, 1 / 1.4143])*0.7
    x, y = (1-w)/2, (1-h)/2
    tx = fg.text(x+w/2, y/2, text)
    tx.set_fontfamily('monospace')
    tx.set_horizontalalignment('center')
    tx.set_verticalalignment('center')
    tx.set_fontsize('large')
    return tx

###############################################################################

SELECT = [
    "CONVERGENCE",  # 0
    "COMPARISON",   # 1
    ][0]

if __name__ == "__main__":

    if SELECT == "CONVERGENCE":
        
        from matplotlib.artist import Artist

        D = Document()
        D.opendocument("../local/plot.pdf")

        w, h, g, e, r = 0.3, 0.05, 0.05, 0.0, 0.45

        NN = []

        # for n in [1, 5, 10, 50, 100, 200]:
        for n in [3, 10, 30, 100, 300]:
            NN.append(
                ([r],                       # shield (radius)
                [+e, (+h+g)/2.0, w, h],     # top    (x,y,w,h)
                [-e, (-h-g)/2.0, w, h],     # bottom (x,y,w,h)
                [n]*8)                      # series (steps, ...)
            )

        X, Y = [], []

        for SH, C1, C2, N in NN:

            S = SolverTwoDimensions()

            S.addMap("SH", DiskAperture(*SH)) # SHIELD
            S.addMap("C1",   PlateSolid(*C1)) # PLATE1
            S.addMap("C2",   PlateSolid(*C2)) # PLATE2

            K, C = S.jacobiStepSeries(N, "C1", "C2")

            fg, ax = selectfigure("PLOT")
            ax.semilogx(K, array(C)*1E12, ".-")
            ax.set_xlabel("iterations")
            ax.set_ylabel("Capacitance computed [pF/m]")

            X.append(K[-1])
            Y.append(C[-1])
            print(K[-1], C[-1])

        fg, ax = selectfigure("RESULTS")
        ax.semilogx(X, array(Y)*1E12, ".-")

        ax.set_xlabel("iterations")
        ax.set_ylabel("Capacitance computed [pF/m]")
        ax.set_ylim(0.0, 100.0)

        D.exportfigure("PLOT") 
        D.exportfigure("RESULTS") 

        D.closedocument()

    if SELECT == "COMPARISON":

        DATA = array([
            [
            3, 3.0436270563981022e-12,
            10, 5.1216678036887545e-12,
            30, 5.388440551520647e-12,
            100, 5.384585983847195e-12,
            300, 5.3845859348425295e-12,
            ],
            [
            7, 1.3545469497385262e-11,
            21, 1.4124039507603287e-11,
            61, 1.3673979676305713e-11,
            200, 1.3576555499625977e-11,
            600, 1.3576477010394754e-11,
            ],
            [
            10, 2.2063907663026516e-11,
            32, 2.3459003810195627e-11,
            92, 2.3688678652095264e-11,
            300, 2.386603938645882e-11,
            900, 2.3916223989714314e-11,
            ],
            [
            15, 3.4069942354760546e-11,
            43, 3.634845120054825e-11,
            123, 3.638886271302465e-11,
            401, 3.6355243645236484e-11,
            1200, 3.637257009179239e-11,
            ],
            [
            19, 5.0055143884454275e-11,
            54, 5.148642724515242e-11,
            154, 5.182449505402855e-11,
            500, 5.176298878944347e-11,
            1500, 5.1766003171494677e-11,
            ],
            [
            23, 5.893946748795823e-11,
            65, 5.780413511747902e-11,
            184, 5.78213629295024e-11,
            600, 5.784758349382356e-11,
            1800, 5.785180109663687e-11,
            ],
            [
            27, 6.04384253286411e-11,
            76, 5.851355714867663e-11,
            215, 5.817602289962751e-11,
            700, 5.7952018021336755e-11,
            2100, 5.7888629596266815e-11,
            ],
            [
            31, 6.332311209754091e-11,
            87, 6.048215759261286e-11,
            247, 5.992852197701082e-11,
            801, 5.971289191392915e-11,
            2400, 5.967733448169437e-11,
            ],
            [
            35, 6.832368567458584e-11,
            98, 6.308824560138485e-11,
            275, 6.205336337416299e-11,
            907, 6.169605168399679e-11,
            2700, 6.158444576691678e-11,
            ],
            [
            39, 7.113133410465773e-11,
            100, 6.415771689528503e-11,
            308, 6.27469380545123e-11,
            1000, 6.227336547622629e-11,
            3000, 6.209861777818127e-11,
            ],
            ])

        from numpy import size, shape

        D = Document()
        D.opendocument("../local/plot.pdf")

        fg, ax = selectfigure("RESULTS1")
        N = 8
        for j in range(shape(DATA)[0]):
            X, Y = DATA[j, 0::2], DATA[j, 1::2]
            ax.semilogx(X, Y*1E12, ".-", label = f"{N}x{N}")
            N <<= 1

        ax.set_xlabel("iterations")
        ax.set_ylabel("Computed capacitance[pF/m]")
        ax.legend()
        D.exportfigure("RESULTS1") 

        fg, ax = selectfigure("RESULTS2")
        X, Y = DATA[:, -2], DATA[:, -1]
        ax.semilogx(X, Y*1E12, ".-")

        ax.set_xlabel("iterations")
        ax.set_ylabel("Computed capacitance[pF/m]")
        D.exportfigure("RESULTS2") 

        D.closedocument()
