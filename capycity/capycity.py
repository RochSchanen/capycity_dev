
# DEBUG DISPLAY
VERBOSE         = True
VERBOSE_MAPDATA = False

# TIME RESSOURCE
from time import time
time_start = time()

# (TBI) means "to be implemented"

# VERSION
version = 0.02
"""

    implement multi grid:
        keep the square shape. use concentric increase of the
        resolution with a list of concentric sizes: it is a
        re-implementation of the JacobiStepSeries where the size
        of the higher resolution map is redefined as a subset of
        the original map. The outer coarser map works where the
        potential gradients are small (towards the edges with a
        shield of large size). Maybe, in the future, the size and
        number of sub-map can be arbitrary. This is why a recursive
        scheme is now chosen for future improvements

    improve next:

        - use diff instead of gradient to prevent usage of float
        and keep the use of integers for keeping calculations
        fast and prevent from using too much memory.

        - decrease the maximum value by a factor 4. Sum the arrays
            first and then divide (double shift) the sum by four
            only once instead of four times as it is implemented
            now.

        - allow rectangular shape

        - merge anchor masks

        - check on convergence factor technique and Gauss-Seidel

"""
if VERBOSE: print(f"version {version}")

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
from numpy.random import rand

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
from numpy import uint32 as _DT
# from numpy import uint16 as _DT
# from numpy import uint8  as _DT
""" The data type precision is determined here. Empirically, the
algorithm converges well using the 32 bits unsigned integer type.
The lower 8 bits unsigned integer type is not sufficient but can
be used for debugging purposes. Select the required length of bits
by un-commenting one of previous declarations lines. """

# EXTREMA (zero value and maximum value)
ZV, MV = array([0, -1], _DT)
""" Extrema are computed here. The maximum value might have to be
decreased by a factor 4 in order to prevent overflowing events
during the evaluation of some optimised algorithms. """

def mesh(nn, ll):                                                 # +
    # get grid sizes
    nx, ny = nn
    # get grid lengths
    lx, ly = ll
    # intervals
    dx, dy = lx / nx, ly / ny
    # boundaries
    bx, by = (lx - dx) / 2.0, (ly - dy) / 2.0
    # domains
    Dx, Dy = linspace(-bx, +bx, nx), linspace(-by, +by, ny)
    # done
    return meshgrid(Dx, Dy)

class SolverTwoDimensions():

    def __init__(self, n = 8, l = 1.0):                           # +
        """ GENERAL DESCRIPTION HERE """
        if VERBOSE: print(f"instantiate 'SolverTwoDimensions'")

        # ALL MAPS ARE SQUARE
        nx, ny = n, n
        lx, ly = l, l

        # GRID SIZE (number of grid points)
        # if ny is None: ny = nx
        self.ll = lx, ly
        """ the top grid size is set for all the maps and masks.
        When the grid resolution is changed, maps and masks
        are upgraded to fit the new grid size. Sub-maps can be
        later added on, with higher resolutions and smaller sizes,
        in order to improve accuracy without depleting our computer
        resources too fast. """
        if VERBOSE: print(f"nx = {nx}, ny = {ny}")
        
        # GRID LENGTH (physical length in meters)
        # if ly is None: ly = ny*lx/nx
        self.nn = nx, ny
        """ the length is set for a square shape at the moment.
        Two independent lengths will be added later. The issue
        here is that the grid size and length ratio must be kept
        identical at all time. """
        if VERBOSE: print(f"lx = {lx}, ly = {ly}")

        # MAPS AND MASKS
        self.M = {}
        """ Masks define the geometry of the conducting parts
        that compose the the capacitance to compute. There is
        one potential field map to be evaluated per part. New
        masks are automatically created when a new map is added.
        The ANCHOR mask defines the shape of the part. The masks
        can be build using pre-defined or user-defined methods.
        In the future, masks should be obtained also by merging
        existing masks. """

        # ITERATIONS
        self.k = 1
        """  an iteration index is regularly updated at each step
        of the computation or each increase of the resolution, or
        each creation of a sub-map. """

        # done
        return

    def addPart(self, name, mask = None):                         # +
        """ this adds a new conductor part to the system. A specific
        map is instantiated to solve the Poisson equation for this
        particular part. An ANCHOR mask that defines the shape of the
        part can be provided now or later. A CLEAR mask is created
        automatically from the parts already existing. Some sub-maps
        can be added later, with smaller size and higher resolution,
        when necessary. Where the part is solid, the value of the
        anchor mask should be VM, otherwise, the value should be ZV.
        """
        if VERBOSE: print(f"new map '{name.upper()}'")
        # instantiate new map for this part
        self.M[f"{name}"] = _map(self.nn, self.ll)
        # build a CLEAR mask for this part
        for k in self.M.keys():
            if k == name: continue
            # merge CLEAR maps from other parts
            self.M[f"{name}"].C &= invert(self.M[k].A)
        # add an anchor mask if already provided
        if mask == None: return
        self.addMask(name, mask)
        # done
        return

    def addMask(self, name, mask):                                # +
        # add mask to the part
        self.M[f"{name}"].addNewAnchor(mask)
        # update the CLEAR masks of other parts:
        I = invert(self.M[f"{name}"].A)
        for k in self.M.keys():
            if k == name: continue
            self.M[k].C &= I
        # done
        return

    def addSubnet(self, L):                                       # !

        if VERBOSE: print("subnet: \n --- from:")

        # build the bottom subnet layer
        S = []
        for k in self.M.keys():
            # initialise loop with main map
            t, s = self.M[k], self.M[k].S
            # recurrently loop until no more subnet
            while s: t, s = s, s.S
            # register last subnet
            S.append(t)
            #debug
            if VERBOSE: print(k, t.nn, t.ll)

        if VERBOSE: print(" --- to:")

        # scan through top and bottom
        # layer and create new layer
        for s, k in zip(S, self.M.keys()):
            # get bottom geometry
            nx, ny = s.nn
            lx, ly = s.ll
            # center position
            ox = (nx + 1.0) / 2.0
            oy = (ny + 1.0) / 2.0
            # scale
            ax = lx / nx
            ay = ly / ny
            # new binding box
            l = int(ceil(-L/ax/2+ox))
            b = int(ceil(-L/ay/2+oy))
            r, t = nx+1-l, ny+1-b
            """ forcing the right and top
            values to be centre-symmetric
            in relation to the left and
            bottom values. """
            if VERBOSE: print("VERBOSE_000:", l, r, t, b)
            # compute new subnet geometry
            nx, ny = r-l+1, t+1-b
            lx, ly = ax*nx, ax*ny
            if VERBOSE: print("VERBOSE_001:", nx, ny)
            if VERBOSE: print("VERBOSE_002:", lx, ly)
            # increase resolution
            nx, ny = 2*nx, 2*ny
            """ the resolution in increased
            by increasing the grid sizes and
            keeping the same lengths. """
            # create the new subnet map
            s.S = _map([nx, ny], [lx, ly])
            # reference mask instance from main map
            s.S.addNewAnchor(self.M[k].M)
            # save current data for the quadruplication
            D = copy(s.D[b:t+1, l:r+1])
            # quadruplicate the data set
            s.S.D[1:-1:2, 1:-1:2] = D
            s.S.D[2:-1:2, 1:-1:2] = D
            s.S.D[1:-1:2, 2:-1:2] = D
            s.S.D[2:-1:2, 2:-1:2] = D

        # build subnets clear mask
        for s in S:
            # loop though all the parts but s
            for t in S:
                if s == t: continue
                # merge masks
                s.S.C &= invert(t.S.A)
        """ the clear masks can only be
        build after all the anchor masks
        have been created. """

        # done
        return

    def jacobiSteps(self, n, name1, name2, SavePattern = None):   # +
        """ Perform a series of "n" Jacobi steps on the
        maps associated with part "name1" and part "name2".
        It returns the capacitance value at the end of the
        series or it can record a series of capacitance values
        following the optional "SavePattern" argument. This
        argument is a list of indices for which the capacitance
        has to be calculated. The same number of steps is
        applied to all maps. This is justified if the convergence
        rate depends only/mostly on the size of the grid and
        not so much on the geometry (masks definitions) of the
        parts. This assumption also should allow to devise
        a fixed series of "resolution increments"/"Jacobi steps"
        that should be valid for all geometries and should only
        depends on the grid size. """
        K, C = [], []
        for i in range(n):
            # n Jacobi step iterations over C1 and C2
            self.M[name1].jacobiStep()
            self.M[name2].jacobiStep()
            # save if step match pattern element
            if SavePattern:
                if (self.k+i) in SavePattern:
                    K.append(self.k+i)
                    C.append(self.computeCapacitance(name1, name2))
        # increment step counter
        self.k += n
        # return a single value
        if SavePattern is None:
            return self.computeCapacitance(name1, name2)
        # return saved data
        return K, C

    def jacobiStepSeries(self, S, C1, C2, n = 500):               # !
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

    def computegradient(self, name):                              # !
        """ The electric field is directly computed
        as the gradient of the electric scalar potential """
        nx, ny = self.nn
        lx, ly = self.ll
        # compute intervals
        dx, dy = lx / nx, ly / ny
        # compute gradient without edges
        dDy, dDx = gradient(self.M[f"{name}"].D[1:-1, 1:-1] / MV)
        # fix units and record gradient without edges
        self.M[f"{name}"].dD = dDx / dx, dDy / dy
        # done
        return

    def computeCapacitance(self, name1, name2):                   # !
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

    def incrementResolution(self):                                # !
        """ increases the resolution of the main map
        by a factor two. The length is left unchanged.
        Therefore the resolution is increased by the
        same factor. It can be noticed empirically that
        with larger resolution, the Jacobi iteration
        "propagates" the potential values at a slower
        "pace". At lower resolutions the convergence
        to an approximate solution is a lot quicker.
        In order to increase the convergence speed,
        a gradual incrementation of the resolution
        is performed between series of Jacobi steps.
        The convergence to an accurate solution is
        preserved without compromising to much the
        computing time. This scheme is not used with
        sub-maps, only with the main map. """

        # debug
        if VERBOSE: print("<-", end = "")
        # get current values
        nx, ny = self.nn        
        # update to new values
        nx, ny = nx<<1, ny<<1
        # register the new values
        self.nn = nx, ny        
        # debug
        if VERBOSE: print("map.D-", end = "")
        # QUADRUPLICATE ALL DATA MAPS AND UPDATE LOCALS
        for k in self.M.keys():
            # update grid size
            self.M[k].nn = nx, ny
            # save current data for the upgrade
            D = copy(self.M[k].D[1:-1, 1:-1])
            """ edges are ignored. """
            # reserve memory for the new data set
            self.M[k].D = full([nx+2, ny+2], ZV, _DT)
            """ the new data set is initialised to
            zeros. This is clearing up the edges.
            The rest of the array is defined during
            the 'quadruplication' """
            # debug
            if VERBOSE: print("quad-", end = "")
            # quadruplicate the data set
            self.M[k].D[1:-1:2, 1:-1:2] = D
            self.M[k].D[2:-1:2, 1:-1:2] = D
            self.M[k].D[1:-1:2, 2:-1:2] = D
            self.M[k].D[2:-1:2, 2:-1:2] = D
        # debug
        if VERBOSE: print("map.A-", end = "")
        # RE-BUILD ALL ANCHOR MASKS
        for k in self.M.keys():
            # re-build ANCHOR mask
            self.M[k].buildAnchorMask()
        # debug
        if VERBOSE: print("map.C-", end = "")
        # RE-BUILD ALL CLEAR MASKS
        for k in self.M.keys():
            # reserve new memory for the new CLEAR mask
            self.M[k].C = full([nx+2, ny+2], MV, _DT)
            # loop though all the other parts
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

    def __init__(self, nn, ll):                                   # !
        """ GENERAL DESCRIPTION """
        if VERBOSE: print(f"instantiate '_map'")
        # save parameters locally
        self.nn, self.ll = nn, ll
        # get detailed geometry for declarations
        nx, ny = nn
        # debug
        if VERBOSE: print(f"nx = {nx}, ny = {ny}")
        # create new data map (potential between conductors)
        self.D = full([ny+2, nx+2], ZV, _DT)
        # fill the initial value with random values
        # self.D[1:-1, 1:-1] = (MV*rand(nx, ny)).astype(_DT)
        # create new CLEAR map (zero potential from other parts)
        self.C = full([ny+2, nx+2], MV, _DT)
        # create new anchor map (unit potential from this part)
        self.A = full([ny+2, nx+2], ZV, _DT)
        # debug
        if VERBOSE_MAPDATA: print(f"map:")
        if VERBOSE_MAPDATA: print(f"{self.D}")
        # sub-map
        self.S = None
        # mask generator methods
        self.M = None
        # done
        return

    def addNewAnchor(self, mask):                                 # !
        # register a new mask
        self.M = mask
        # self.M.append(mask)
        """ so far, only a single mask can be used
        later, more of them could be added and
        merged together using the buildAnchorMask()
        method """
        self.buildAnchorMask()
        # done
        return

    def buildAnchorMask(self):                                    # !
        # MERGE ANCHOR MAPS
        self.A = self.M.mask(self.nn, self.ll)
        # self.A |= self.M.mask(nn, ll)
        """ when there is more than one mask constructor
        to build the mask, all the masks generated must be
        merged together here."""
        if VERBOSE_MAPDATA: print(f"anchor map:")
        if VERBOSE_MAPDATA: print(f"{self.A}:")
        # done
        return

    def getboundary(self):
        """ the outer layers of the
        current map are averaged and
        returned to the caller """
        # left
        l  = self.D[1:-1:2, +1]//4
        l += self.D[2:-1:2, +1]//4
        l += self.D[1:-1:2, +2]//4
        l += self.D[2:-1:2, +2]//4
        # right
        r  = self.D[1:-1:2, -2]//4
        r += self.D[2:-1:2, -2]//4
        r += self.D[1:-1:2, -3]//4
        r += self.D[2:-1:2, -3]//4
        # top
        t  = self.D[-2, 1:-1:2]//4
        t += self.D[-2, 2:-1:2]//4
        t += self.D[-3, 1:-1:2]//4
        t += self.D[-3, 2:-1:2]//4
        # bottom
        b  = self.D[+1, 1:-1:2]//4
        b += self.D[+1, 2:-1:2]//4
        b += self.D[+2, 1:-1:2]//4
        b += self.D[+2, 2:-1:2]//4

        t = full(len(self.D[1:-1:2, +1]), MV, _DT)

        # done
        return l, r, t, b

    def setboundary(self, l, r, t, b):
        """ the given boundaries data are automatically centred
        and set to the outer boundary of this map. If this method
        is never called, the outer boundary is left at the the
        zero value ZV by default. """
        # left
        self.D[1:-1:2,  0] = l
        self.D[2:-1:2,  0] = l
        # right
        self.D[1:-1:2, -1] = r
        self.D[2:-1:2, -1] = r
        # top
        self.D[-1, 1:-1:2] = t
        self.D[-1, 2:-1:2] = t
        # bottom
        self.D[ 0, 1:-1:2] = b
        self.D[ 0, 2:-1:2] = b
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

        """
        Also, the algorithm of Gauss-Seidel could be
        implemented instead (explicit loop). However,
        this requires an implementation in cpython or
        some pre-compilation with optimisation.
        """        
        
        """
        Another thing. It may be worth investigate the
        use of the graphic card ability of processing
        these operations very fast: like operation on
        images. All the computation are linear and most
        likely can be handled very fast by a graphics
        card.
        """

        """ check for sub map and setup boundaries between 
        map interfaces """

        # link to sub subnet
        if self.S:

            # get geometry
            mx, my = self.nn
            nx, ny = self.S.nn
            px = (mx - (nx // 2)) // 2 + 1
            py = (my - (ny // 2)) // 2 + 1
            """ the + 1 in previous lines is
            due to the presence of the edge
            layer. """

            # set inner boundary values from subnet
            l, r, t, b = self.S.getboundary()
            self.D[py:-py,  px     ] = l
            self.D[py:-py,  mx-px+1] = r
            self.D[my-py+1, px:-px ] = t
            self.D[py,      px:-px ] = b

            # set outer boundary values of subnet
            l = self.D[py:-py,  px-1   ]
            r = self.D[py:-py,  mx-px+2]
            t = self.D[my-py+2, px:-px ]
            b = self.D[py-1,    px:-px ]
            self.S.setboundary(l, r, t, b)

        # JACOBI STEP:

        # compute shifted arrays
        l = copy(self.D[1:-1, 0:-2]) # left
        r = copy(self.D[1:-1, 2:  ]) # right
        t = copy(self.D[0:-2, 1:-1]) # top
        b = copy(self.D[2:  , 1:-1]) # bottom

        # compute neighbour average
        self.D[1:-1, 1:-1] = l//4 + r//4 + t//4 + b//4
        """ should be able to reduce the number of divisions
        by summing first and dividing (double shifting) only
        once after. However, be aware of the possible overflow
        of our unsigned integer type. This could work if we
        divided MV by a factor four. That would prevent a
        possible overflow: the idea needs to be tested.
        """

        # reset boundaries
        self.applyMasks()

        if self.S:
            # compute all subnet jacobi step:
            self.S.jacobiStep()
            """ All subnet maps will be
            computed recurrently. """

        # done
        return

    def decor(self):
        # get array length
        lx, ly = self.ll
        # get boundaries
        l, r, t, b = -lx/2, +lx/2, +ly/2, -ly/2
        # make decor
        c = Rectangle(
                (l, b), r-l, t-b,
                edgecolor = (0.8, 0.4, 0.4),
                linestyle = "-.",# line style
                fill      = False,
                linewidth = 1,
            )
        return c

class DiskSolid():

    # parameters
    def __init__(self, r):

        if VERBOSE: print("instantiate 'DiskSolid'")

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
        # debug
        if VERBOSE: print(R)
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
    # map figure is the same as selectfigure()
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

def showResolution(ax, m):
    # get geometry
    X, Y = mesh(m.nn, m.ll)
    nx, ny = m.nn
    lx, ly = m.ll
    dx, dy = lx / nx, ly / ny
    # loop through four corners
    for i, j in [(1, 1), (1, -2), (-2, 1), (-2, -2)]:
        # compute coordinates of a rectangle of size dx, dy
        x, y = X[i, j], Y[i, j]
        l, r, t, b = x-dx/2, x+dx/2, y+dy/2, y-dy/2
        # display the patch
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
    X, Y = mesh(solver.nn, solver.ll)
    # get normalised data without the edges
    D = solver.M[mapname].D[1:-1, 1:-1] / MV
    # create filled contour map
    QCS = ax.contourf(X, Y, D, *args, **kwargs)
    # adjust ticks
    lx, ly = solver.ll
    MX, SX = _getTickPositions(-lx/2.0, +lx/2.0, 9)
    MY, SY = _getTickPositions(-ly/2.0, +ly/2.0, 9)
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
    X, Y = mesh(solver.nn, solver.ll)
    # get normalised data without the edges
    D = solver.M[mapname].D[1:-1, 1:-1] / MV
    # create filled contour map
    QCS = ax.pcolormesh(
        X, Y, D, 
        vmin = 0.0, vmax = 1.0,
        shading = 'auto',
        rasterized = True)
    # adjust ticks
    lx, ly = solver.ll
    MX, SX = _getTickPositions(-lx/2.0, +lx/2.0, 9)
    MY, SY = _getTickPositions(-ly/2.0, +ly/2.0, 9)
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
    if VERBOSE: showResolution(ax, solver.M[k])
    # done
    return fg, ax, bx

def splot(ax, m):
    VX, VY = 0.5, 0.5
    # VX, VY = 0.0, 0.0
    while m.S:
        # get the mesh coordinates
        X, Y = mesh(m.S.nn, m.S.ll)
        # get normalised data without the edges
        D = m.S.D[1:-1, 1:-1] / MV
        # create filled contour map
        QCS = ax.pcolormesh(
            X+VX, Y+VY, D,
            vmin = 0.0, vmax = 1.0,
            shading = 'auto', 
            rasterized = True)
        # set colour map
        QCS.set_cmap(cm.inferno)
        # add decor
        if VERBOSE: ax.add_artist(m.S.decor())
        if VERBOSE: showResolution(ax, m.S)
        # get next subnet
        m = m.S
    # done
    return

def headerText(text, fg):
    w, h = array([1, 1 / 1.4143])*0.7
    x, y = (1-w)/2, (1-h)/2
    tx = fg.text(x+w/2, 3*y/2+h, text)
    tx.set_fontfamily('monospace')
    tx.set_horizontalalignment('centre')
    tx.set_verticalalignment('centre')
    tx.set_fontsize('large')
    return tx

def footerText(text, fg):
    w, h = array([1, 1 / 1.4143])*0.7
    x, y = (1-w)/2, (1-h)/2
    tx = fg.text(x+w/2, y/2, text)
    tx.set_fontfamily('monospace')
    tx.set_horizontalalignment('centre')
    tx.set_verticalalignment('centre')
    tx.set_fontsize('large')
    return tx

#####################################################################

SELECT = [
    "SERIES-CONVERGENCE",  # 0
    "SERIES-COMPARISON",   # 1
    "MULTI-GRID",          # 2
    ][2]

if __name__ == "__main__":

    time_compute = time()

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

    if SELECT == "MULTI-GRID":

        SZ = 32

        NN = [
            (500, 2.0),
            (100, 1.0),
            # (100, 2.0),
            # (100, 2.0),
            # (100, 2.0),
            ]

        # INIT SOLVER
        N, L = NN[0]
        S = SolverTwoDimensions(n = SZ, l = L)
        S.addPart("C", PlateSolid(0.0, 0.0, 0.2, 0.01))
        S.addPart("S", DiskAperture(0.95))
        # main series
        for i in range(N):
            S.M["C"].jacobiStep()
            S.M["S"].jacobiStep()
        # subnet series
        for N, L in NN[1:]:
            # ADD SUBNET
            S.addSubnet(L)
            for i in range(N):
                S.M["C"].jacobiStep()
                S.M["S"].jacobiStep()

        # DISPLAY
        fg, ax, bx = mplot(S, "P1", "C")
        if S.M["C"].S: splot(ax, S.M["C"])

        fg, ax, bx = mplot(S, "P2", "S")
        if S.M["S"].S: splot(ax, S.M["S"])

        # BUILD PDF
        D = Document()
        D.opendocument("../local/plot.pdf")
        for p in ["P1", "P2"]:
            D.exportfigure(p)
        D.closedocument()

    # done
    time_done = time()
    print(f"load in {time_compute - time_start}s")
    print(f"compute in {time_done - time_compute}s")
