
# from package: "https://numpy.org/"
from numpy import array
from numpy import full
from numpy import ceil
from numpy import floor
from numpy import asarray

# from package "https://python-pillow.org/"
from PIL import Image
from PIL import ImageDraw

# DATA TYPE
""" The data type precision is fixed here. Empirically, the algorithm converges
nicely when using 32 bits unsigned integers. A lower 8 bits unsigned integer
type can be used when debugging."""

from numpy import uint8  as _DT
# from numpy import uint16  as _DT
# from numpy import uint32 as _DT

# EXTREMA (zero value and maximum value)

ZV, MV = array([0, -1], _DT)

# SOLVER CLASS

class SolverTwoDimensions():

    def __init__(self, px = 3, py = None):

        if VERBOSE: print(f"instanciate 'SolverTwoDimensions'")

        # LENGTH

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

        # MAPS AND MASKS

        """ Masks define the geometry of the conducting materials
        that compose the the capacitance to calculate. There is
        one potential field map to compute per conductor/mask. New
        maps are automatically created when new masks are added.
        Masks are build using pre-defined or user-defined methods.
        A single mask can be obtained by merging several mask."""

        self.M = {}

        return

    def addMap(self, name, mask = None):
        # instanciate new map class
        self.M[name] = _map(self.nn)
        # no mask, done
        if mask is None: return
        # add mask
        self.addMask(name, mask)
        # done
        return

    def addMask(self, name, mask):
        # use the map method
        self.M[name].setMask(mask, self.nn, self.ll)
        # done
        return

class _map():

    def __init__(self, nn):
        
        if VERBOSE: print(f"instanciate '_map'")

        # get grid size
        nx, ny = nn
        if VERBOSE: print(f"nx = {nx}, ny = {ny}")
        
        # new maps
        self.D = full([nx+2, ny+2], ZV, _DT) # DATA MAP
        self.C = full([nx+2, ny+2], MV, _DT) # CLEAR MAP
        self.A = full([nx+2, ny+2], ZV, _DT) # ANCHOR MAP
        if VERBOSE: print(f"map:")
        if VERBOSE: print(f"{self.D}")

        # masks
        self.M = None

        # done
        return

    def setMask(self, mask, nn, ll):
        # register the class instance
        self.M = mask
        # and setup anchor map
        self.refreshMask(nn, ll)
        # done
        return

    def refreshMask(self, nn, ll):
        # use the mask method
        # self.A |= self.M.mask(nn, ll)
        self.A = self.M.mask(nn, ll)
        if VERBOSE: print(f"anchor map:")
        if VERBOSE: print(f"{self.A}:")
        # done
        return

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

VERBOSE = True

S = SolverTwoDimensions()
S.addMap("shell", DiskAperture(0.5))

# S.addMap("core", Disk(0.0, 0.0, 0.1))
# S.jacobiSteps(100)
# S.capacitance("core", "shell")
# cummulate steps
