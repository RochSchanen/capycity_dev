
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
unsigned integer type is not sufficient but can be used for debugging.
Select the required length of bits by un-commenting one of following three
lines: """

# from numpy import uint8  as _DT
# from numpy import uint16  as _DT
from numpy import uint32 as _DT

# EXTREMA (zero value and maximum value)
ZV, MV = array([0, -1], _DT)

# DEBUG DISPLAY
VERBOSE_MAP = False
VERBOSE = False

class SolverTwoDimensions():

    def __init__(self, px = 3, py = None):
        if VERBOSE: print(f"instantiate 'SolverTwoDimensions'")

        # LENGTH

        """ the length is set for a square shape only so far.
        therefore, just use square grids for the moment. it will
        be extended to two independent lengths later. The issue
        is that the grid ratio and length ratio must be kept
        identical """
        
        self.ll = 1.0
        if VERBOSE: print(f"length = {self.ll}")

        # GRID
        
        """ the grid size is set for all the maps and masks.
        When the grid resolution is changed, maps and masks
        are upgraded to fit the new grid size. """
        
        # compute initial size
        if py is None: py = px
        if VERBOSE: print(f"px = {px}, py = {py}")
        nx, ny = 2**px, 2**py
        if VERBOSE: print(f"nx = {nx}, ny = {ny}")

        # locals
        self.pp = px, py     # size in power of two
        self.nn = nx, ny     # size in unit cell

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

        """ create a new map to hold a new component to the capacitance.
        A CLEAR map and a ANCHOR map are automatically added for the
        computation. The ANCHOR map should contain the shape of the solid.
        Where the maps is solid, the value is VM, otherwise, the value is ZV.
        A CLEAR map is build from all the other maps already available. """
        
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
        # intervals (intervals should be identical)
        dx, dy = self.ll / nx, self.ll / ny
        # limits
        lx, ly = (self.ll - dx) / 2.0, (self.ll - dy) / 2.0
        # domains
        Dx, Dy = arange(-lx, +lx + dx, dx), arange(-ly, +ly + dy, dy)
        # mesh
        self.mesh = meshgrid(Dx, Dy)
        # done
        return

    def jacobiSteps(self, n):

        """ using the same number of steps on every map is justified if
        the convergence rate depends mostly on the size/resolution of the
        grid and not so much on the geometry of the system. This should 
        also allow to devise an alternating series of incrementing the
        resolution followed by a steps series that could be used for all
        cases: see the next method """

        for i in range(n):
            # n Jacobi step iterations over every maps
            for k in self.M.keys():
                self.M[k].jacobiStep()
        # increment step counter
        self.k += n
        # done
        return

    def jacobiStepSeries(self, S, C1, C2, n = 100):

        """ S is a list of step series numbers. After each step series,
        the resolution of the grid is incremented by a factor two. This
        is why the resolution of the grid is always a power of two. C1
        and C2 are the names of the two conductors for which we compute
        the mutual capacitance. The other conductors affect the computation
        only by their presence (geometry) and the value of their potential
        is kept at zero (See theoretical considerations from the report).
        At the end of each series, the total number of steps and the value
        of the capacitance is computed and recorded in list K and R. This
        allows to monitor the convergence the series. """

        if VERBOSE:
            nx, ny = self.nn
            print(f"resolution: {nx}x{ny}")

        i, a = 0, S[0]
        K, R = [], []
        # record first point
        K.append(self.k)  
        R.append(self.computeCapacitance(C1, C2))
        # compute steps series
        m = ceil(log(sum(S)) / log(10))
        I = diff(logspace(1, m, n, dtype = 'int'))
        # loop through series
        for n in I:
            self.jacobiSteps(n)
            # record new point
            K.append(self.k)  
            R.append(self.computeCapacitance(C1, C2))
            # check if new resolution triggered
            if self.k < a: continue
            if i+1 == len(S): break
            if VERBOSE: print(f"{K[-1]:6}, {R[-1]*1E12:.3f}")
            i, a = i+1, a + S[i+1]
            self.incrementResolution()
            if VERBOSE:
                nx, ny = self.nn
                print(f"resolution: {nx}x{ny}")
            K.append(self.k)
            R.append(self.computeCapacitance(C1, C2))
        # last point
        K.append(self.k)
        R.append(self.computeCapacitance(C1, C2))
        # done 
        return K, R

    def computegradient(self, name):

        """ The electric field is directly computed
        as the gradient of the electric scalar potential """
        
        nx, ny = self.nn
        # compute interval (intervals should be identical)
        dx, dy = self.ll / nx, self.ll / ny
        # compute gradient without edges
        dDy, dDx = gradient(self.M[name].D[1:-1, 1:-1] / MV)
        # fix units and record gradient without edges
        self.M[name].dD = dDx / dx, dDy / dy
        # done
        return

    def computeCapacitance(self, name1, name2):
        nx, ny = self.nn
        # compute interval (intervals should be identical)
        dx, dy = self.ll / nx, self.ll / ny
        # compute map gradient (electric field)
        self.computegradient(name1)
        # get component fields
        EX1, EY1 = self.M[name1].dD
        self.computegradient(name2)
        # get component fields
        EX2, EY2 = self.M[name2].dD
        # compute 'field' dot products
        integrand = MLT(EX1, EX2) + MLT(EY1, EY2)
        # compute numeric integral
        integral = SUM(integrand)*dx*dy
        # normalise units and return the mutual
        # capacitance as a positive value
        C = -integral*EPS0
        # done
        return C

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
        px, py = self.pp
        # update values
        px, py = px + 1, py + 1
        nx, ny = 2**px, 2**py
        # update local parameters
        self.pp = px, py
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
            array is defined during the "quadruplication """
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

        # done`
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
        # compute index origin
        ox = (nx + 1.0) / 2.0
        oy = (ny + 1.0) / 2.0
        # compute interval
        ax = ll / nx
        ay = ll / ny
        if VERBOSE: print(f"nx ny ll ox oy = ", end ='')
        if VERBOSE: print(nx, ny, ll, ox, oy)
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
        # compute index origin
        ox = (nx + 1.0) / 2.0
        oy = (ny + 1.0) / 2.0
        # compute interval
        ax = ll / nx
        ay = ll / ny
        if VERBOSE: print(f"nx ny ll ox oy = ", end ='')
        if VERBOSE: print(nx, ny, ll, ox, oy)
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
        # compute index origin
        ox = (nx + 1.0) / 2.0
        oy = (ny + 1.0) / 2.0
        # compute interval
        ax = ll / nx
        ay = ll / ny
        if VERBOSE: print(f"nx ny ll ox oy = ", end ='')
        if VERBOSE: print(nx, ny, ll, ox, oy)
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
        # compute index origin
        ox = (nx + 1.0) / 2.0
        oy = (ny + 1.0) / 2.0
        # compute interval
        ax = ll / nx
        ay = ll / ny
        if VERBOSE: print(f"nx ny ll ox oy = ", end ='')
        if VERBOSE: print(nx, ny, ll, ox, oy)
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
    dx, dy = solver.ll / nx, solver.ll / ny
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
    M, S = _getTickPositions(-solver.ll/2.0, +solver.ll/2.0, 9)
    ax.set_xticks(M)
    ax.set_yticks(M)
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
    M, S = _getTickPositions(-solver.ll/2.0, +solver.ll/2.0, 9)
    ax.set_xticks(M)
    ax.set_yticks(M)
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

if __name__ == "__main__":

    from matplotlib.artist import Artist

    D = Document()
    D.opendocument("../local/plot.pdf")

    w, h, g = 0.3, 0.05, 0.05

    NN = []
    CC = []

    sh = 0.0
    # for sh in arange(0.0, w/2, w/4/16):
    for r in range(1, 8):    
        NN.append(
        ([0.45],                 # shield (radius)
        [+sh, (+h+g)/2.0, w, h], # top    (x,y,w,h)
        [-sh, (-h-g)/2.0, w, h], # bottom (x,y,w,h)
        [100]*r)                 # series (steps, ...)
        )

    for SH, C1, C2, N in NN:

        S = SolverTwoDimensions(3)

        S.addMap("SH", DiskAperture(*SH)) # shield
        S.addMap("C1",   PlateSolid(*C1)) # plate 1
        S.addMap("C2",   PlateSolid(*C2)) # plate 2

        K, R = S.jacobiStepSeries(N, "C1", "C2")
        CC.append(R[-1])

        fg, ax, bx = vplot(S, "FIGURE", "C1")
    
        nx, ny = S.nn

        txt  = f"width = {w}mm\n"
        txt += f"height = {h}mm\n"
        txt += f"gap = {g}mm\n"
        txh = headerText(txt, fg)

        txt  = f"total steps = {K[-1]:6}\n"
        txt += f"capacity  = {R[-1]*1E12:.3f}pF\n"
        txf = footerText(txt, fg)

        ax.set_title(f"{nx} X {ny}\n")

        D.exportfigure("FIGURE") 
        ax.clear()

        Artist.remove(txh)
        Artist.remove(txf)

    D.closedocument()



#########################################################################
###                 SETS OF SELECTED RESULTS                          ###
#########################################################################



# SET 1 #######################################################################

"""

8 X 8 -> 4096 X 4096

NN = [
    (3, [2]*10),
    (3, [5]*10),
    (3, [10]*10),
    (3, [50]*10),
    (3, [100]*10),
    (3, [500]*10),
    (3, [1000]*10),
]

 steps, capacitance

    20, 0.352   0.35
    50, 29.103  29.1
    91, 37.478  37.5
   512, 37.646  37.6
   991, 37.448  37.4
  5327, 37.016  37.0
  9991, 36.981  37.0

"capycity_dev/outputs/CoaxialCapacitorConvergenceSeries.pdf"

"""

# SET 2 #######################################################################

"""

NN = [
    #     8  16   32   64   128   256   512  1024  2048
    (3, [20, 75, 200, 400, 2000, 1200, 2000, 2000, 10000]),
]

for x1, x2, y1, y2, tt in [
    (   1,    50, -5.0, 30.0,    "8 x 8"),
    (  10,   150, 25.0, 38.0,   "16 x 16"),
    (  80,   400, 29.0, 38.0,   "32 x 32"),
    ( 250,   800, 29.0, 33.0,   "64 x 64"),
    ( 500,  4000, 32.0, 37.0,  "128 x 128"),
    (2500,  5000, 35.5, 37.0,  "256 x 256"),
    (3500,  7000, 36.4, 37.0,  "512 x 512"),
    (5000, 10000, 36.5, 37.0, "1024 x 1024"),
    (5000, 20000, 36.6, 37.0, "2048 x 2048"),
    ]:
    ax.set_xlim(x1, x2)
    ax.set_ylim(y1, y2)
    ax.set_title(tt)
    D.exportfigure("FIGURE") 

 steps, capacitance

   8x8,       21,  26.128
  16x16,     103,  33.220
  32x32,     303,  29.493
  64x64,     713,  32.648
 128x128,   2906,  35.712
 256x256,   4220,  36.706
 512x512,   6126,  36.558
1024x1024,  8102,  36.699
2048x2048, 18729, 36.912

[Finished in 3665.9s]

One minute past an hour: mostly the last series: 2048x2048: 10000 points

"capycity_dev/outputs/CoaxialCapacitorConvergenceSeries-2.pdf"

"""

# SET 3 #######################################################################

"""

NN = [
    #     8  16  32   64  128  256  512 1024 2048
    (3, [20, 20, 20, 100, 100, 200, 200, 200, 200]),
]

for x1, x2, y1, y2, tt in [
    (   1,    50, -5.0, 30.0,    "8 x 8"),
    (  10,   150, 25.0, 38.0,   "16 x 16"),
    (  80,   400, 29.0, 38.0,   "32 x 32"),
    ( 250,   800, 29.0, 33.0,   "64 x 64"),
    ( 500,  4000, 32.0, 37.0,  "128 x 128"),
    (2500,  5000, 35.5, 37.0,  "256 x 256"),
    (3500,  7000, 36.4, 37.0,  "512 x 512"),
    (5000, 10000, 36.5, 37.0, "1024 x 1024"),
    (5000, 20000, 36.6, 37.0, "2048 x 2048"),
    ]:
    ax.set_xlim(x1, x2)
    ax.set_ylim(y1, y2)
    ax.set_title(tt)
    D.exportfigure("FIGURE") 

    21, 26.128
    40, 33.570
    61, 30.552
   165, 32.713
   275, 35.732
   488, 36.788
   696, 36.623
   860, 36.774
  1063, 37.024

[Finished in 74.0s]

only fourteen seconds past one minute!

"capycity_dev/outputs/CoaxialCapacitorConvergenceSeries-3.pdf"

"""

# SET 4 #######################################################################

"""

NN = [
    (3, [20]),       #   8
    (4, [50]),       #   16
    (5, [250]),      #   32
    (6, [1000]),     #   64
    (7, [4000]),     #   128
    (8, [15000]),    #   256
    (9, [60000]),    #   512
    (10, [100000]),  #  1024
]

    20, 26.131  8
    50, 33.247  16
   250, 29.481  32
   991, 32.644  64
  4028, 35.708  128
 15547, 36.701  256
 62793, 36.548  512
 99991, 34.802 1024 - clearly not converging yet (see graph)

[Finished in 7567.8s]

Very slow convergence for high resolutions: 5 mins past two hours.

"capycity_dev/outputs/CoaxialCapacitorConvergenceSeries-4.pdf"

"""
