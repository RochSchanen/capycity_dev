
# LIBRARIES

# from package: "https://numpy.org/"
from numpy import full
from numpy import ceil
from numpy import copy
from numpy import floor
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
from matplotlib.pyplot import Circle
from matplotlib.pyplot import figure
from matplotlib.pyplot import Rectangle
from matplotlib.pyplot import fignum_exists
from matplotlib.backends.backend_pdf import PdfPages

# DATA TYPE
""" The data type precision is fixed here. Empirically, the algorithm converges
nicely when using 32 bits unsigned integers. A lower 8 bits unsigned integer
type can be used when debugging."""
# from numpy import uint8  as _DT
# from numpy import uint16  as _DT
from numpy import uint32 as _DT

# EXTREMA (zero value and maximum value)
ZV, MV = array([0, -1], _DT)

# DEBUG DISPLAY
VERBOSE = False

# SOLVER CLASS
class SolverTwoDimensions():

    def __init__(self, px = 3, py = None):
        if VERBOSE: print(f"instanciate 'SolverTwoDimensions'")

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
        # compute intial size
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

    """ create a new map to hold a new component to the capacitance. A CLEAR
    map and a ANCHOR map are automatically added for the computation. The
    ANCHOR map should contain the shape of the solid. where the maps is solid
    the value is VM, otherwise, the value is ZV. A clear map is build from
    the other maps. """
    def addMap(self, name, mask = None):
        if VERBOSE: print(f"new map '{name.upper()}'")
        # instanciate new map class
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
        # n jacobi step iterations over every maps
        for i in range(n):
            for k in self.M.keys():
                self.M[k].jacobiStep()
        """ this approach is justified if convegence rate mostly depends on
        the size of the grid and not the gemometry of the system. This should
        allow to devise a fixed 'step-increment' series to be used for all
        cases """
        # incremenb step counter
        self.k += n
        # done
        return

    def computegradient(self, name):
        nx, ny = self.nn
        # compute interval (intervals shoud be identical)
        dx, dy = self.ll / nx, self.ll / ny
        # compute gradient without edges
        dDy, dDx = gradient(self.M[name].D[1:-1, 1:-1] / MV)
        # fix units and record gradient without edges
        self.M[name].dD = dDx / dx, dDy / dy
        # done
        return

    def computeCapacitance(self, name1, name2):
        nx, ny = self.nn
        # compute interval (intervals shoud be identical)
        dx, dy = self.ll / nx, self.ll / ny
        # compute map gradient (electric field)
        self.computegradient(name1)
        EX1, EY1 = self.M[name1].dD
        self.computegradient(name2)
        EX2, EY2 = self.M[name2].dD
        # compute dor products
        integrand = MLT(EX1, EX2) + MLT(EY1, EY2)
        # compute numeric integral
        integral = SUM(integrand)*dx*dy
        # normalise capacitance 
        CN = -integral*EPS0
        # done
        return CN

    def incrementResolution(self):
        """ increase the resolution by a factor two. One can notice that
        with large resolution, the jacobi iteration "propagates" the potential
        "wave" at a slower "velocity". At lower resolutions the convergence 
        to a solution is quicker. In order to benefit from that fact, an  """

        if VERBOSE: print("<-", end = "")

        # get current values
        px, py = self.pp

        # update values
        px, py = px + 1, py + 1
        nx, ny = 2**px, 2**py

        # update local parameters
        self.pp = px, py
        self.nn = nx, ny

        if VERBOSE: print("map.D-", end = "")

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

        if VERBOSE: print("map.A-", end = "")

        for k in self.M.keys():
            # re-build ANCHOR mask
            self.M[k].refreshMask(self.nn, self.ll)

        if VERBOSE: print("map.C-", end = "")

        for k in self.M.keys():
            # reserve memory for the CLEAR mask
            self.M[k].C = full([nx+2, ny+2], MV, _DT)
            # merge inverse of all other maps
            for l in self.M.keys():
                if k == l: continue
                self.M[k].C &= invert(self.M[l].A)

        if VERBOSE: print("mesh-", end = "")

        # re-compute mesh
        self._computemesh()

        if VERBOSE: print(">")

        # done
        return

class _map():

    def __init__(self, nn):
        if VERBOSE: print(f"instanciate '_map'")
        # get grid size
        nx, ny = nn
        if VERBOSE: print(f"nx = {nx}, ny = {ny}")
        # new maps
        self.D = full([nx+2, ny+2], ZV, _DT) # DATA MAP (BETWEEN COMPONENTS)
        self.C = full([nx+2, ny+2], MV, _DT) # CLEAR MAP (OTHER COMPONENTS)
        self.A = full([nx+2, ny+2], ZV, _DT) # ANCHOR MAP (THIS COMPONENT)
        if VERBOSE: print(f"map:")
        if VERBOSE: print(f"{self.D}")
        # masks
        self.M = None
        # gradients
        self.dD = None
        # done
        return

    """ maybe the two next methods can be merged... or made independent... """

    """ This adds a new solid shape that is registered a merged to the existing
    ANCHOR map.There is only one mask for the moment. """
    def setMask(self, mask, nn, ll):
        # register the class instance
        self.M = mask
        # and setup anchor map
        self.refreshMask(nn, ll)
        # done
        return

    """ when there is more than one function to build the mask, a set of masks
    must be merged together. This is done by refreshMask() """
    def refreshMask(self, nn, ll):
        
        # go through all map building function HERE!
        # go through all map building function HERE!
        # go through all map building function HERE!

        # use the registered mask methods
        # self.A |= self.M.mask(nn, ll)
        self.A = self.M.mask(nn, ll)
        
        if VERBOSE: print(f"anchor map:")
        if VERBOSE: print(f"{self.A}:")
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
        self.D |= self.A
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
    def decorTwoDimensions(self):
        c = Circle(
            (0.0, 0.0),                     # origin
            self.r,                         # size
            edgecolor = (0.7, 0.7, 0.7),    # color
            linestyle = "--",               # line style
            fill = False,                   # fill style
            )
        # done
        return c

class DiskAperture():

    # parameters
    def __init__(self, r):

        if VERBOSE: print("instanciate 'DiskAperture'")

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
    def decorTwoDimensions(self):
        c = Circle(
            (0.0, 0.0),                     # origin
            self.r,                         # size
            edgecolor = (0.7, 0.7, 0.7),    # color
            linestyle = "--",               # line style
            fill = False,                   # fill style
            )
        # done
        return c

class PlateSolid():

    # parameters
    def __init__(self, x, y, w, h):

        if VERBOSE: print("instanciate 'PlateSolid'")

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

    def decorTwoDimensions(self):
        # get geometry
        l, r = self.l, self.r
        t, b = self.t, self.b
        # patch
        R = Rectangle((l, b), r-l, t-b,
            edgecolor = (0.7, 0.7, 0.7),    # color
            linestyle = "--",               # line style
            fill        = False,            # no filling
            )
        # done
        return R

class PlateAperture():

    # parameters
    def __init__(self, x, y, w, h):

        if VERBOSE: print("instanciate 'PlateAperture'")

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

class Document():

    def __init__(self, pathname = None):
        if pathname is not None:
            self._DOC = opendocument(pathname)
        return

    def opendocument(self, pathname):
        self._DOC = PdfPages(pathname)
        return self._DOC

    def exportfigure(self, name):
        if VERBOSE: print(f"export '{name}'")
        fg, ax = selectfigure(name)
        self._DOC.savefig(fg)
        return

    def closedocument(self):
        self._DOC.close()
        return

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
        fg = figure(name)
        # select axes
        ax = fg.get_axes()[0]
    # done
    return fg, ax

S = SolverTwoDimensions(3)
S.addMap("core", DiskSolid(0.1))
S.addMap("shield", DiskAperture(0.45))

R = []
K = []

nx, ny = S.nn
print(f"resolution: {nx}x{ny}")

# first point
R.append(S.computeCapacitance("shield", "core"))
K.append(S.k)

# first series
for n in [1]*30:
    S.jacobiSteps(n)
    R.append(S.computeCapacitance("shield", "core"))
    K.append(S.k)

# resolution series
for s in [[5]*30]*8:
    S.incrementResolution()
    nx, ny = S.nn
    print(f"resolution: {nx}x{ny}")
    R.append(S.computeCapacitance("shield", "core"))
    K.append(S.k)
    for n in s:
        S.jacobiSteps(n)
        R.append(S.computeCapacitance("shield", "core"))
        K.append(S.k)

fg, ax = selectfigure("FIGURE")
# ax.plot(K, array(R)*1E12, ".-")
ax.semilogx(K, array(R)*1E12, ".-")
ax.grid("True")

D = Document()
D.opendocument("../local/plot.pdf")
D.exportfigure("FIGURE") 
D.closedocument()
