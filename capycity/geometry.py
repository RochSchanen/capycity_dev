# file: geometry.py
# author: Roch Schanen
# date: 2022-11-20
# content: definition classes of the geometries used by the solver

# from package: "https://matplotlib.org/"
from matplotlib.pyplot import Circle
from matplotlib.pyplot import vlines

# from package: "https://numpy.org/"
from numpy import array
from numpy import empty
from numpy import shape
from numpy import square

class _geometry():
    """ generic class:
    1) Thge "__init__()" method invoques the user method "init()"
    2) The "mask()" method invoques the user method "discr()" to detect
    wether a point is whitin the mask boundaries.
    """

    def __init__(self, *args, **kwargs):
        # call user init
        self.init(*args, **kwargs)
        # done
        return

    def mask(self, D):
        # get array size
        n, n = shape(D)
        # get array type
        t = D.dtype
        # get extrema from data type
        ZER, MAX = array([0, -1], t)
        # init transparent masks
        M = empty((n, n), dtype = t)
        # compute grid origin
        o = (n - 1) / 2.0
        # compute scaling factor
        a = 1.0 / (n - 2)
        # build masks:
        for i in range(n):
            # compute coordinate
            y = (i-o)*a
            for j in range(n):
                # compute coordinate
                x = (j-o)*a
                # compute dicriminant
                d = self.discr(x, y)
                # apply discriminant
                M[i, j] = [ZER, MAX][int(d)]
        # done
        return M

    def init(self):
        """ the init method is used to pass the parameters that define
        the mask geometry (height, width, thickness, etc...)
        """
        pass

    def discr(self, x, y):
        """ the discr method is used to discriminate the points that belong
        or not to the mask boundaries (disks, rectangles, etc...)
        """
        return False

#####################################################################    DISK
###                            DISK                              ####
#####################################################################

class Disk(_geometry):

    # parameters
    def init(self, r):
        self.r = r
        self.R = square(r)
        return

    # discriminant
    def discr(self, x, y):
        T = square(x) + square(y)
        return (T < self.R)

    ##########################################################
    """ decors are created only on the current figure. It must
    be created on each figure that wants to include that decor.
    """

    def contour2D_decor(self):
        
        # build decor
        c = Circle(
            (0.0, 0.0),                     # origin
            self.r,                         # size
            edgecolor = (0.7, 0.7, 0.7),    # color
            linestyle = "--",               # line style
            fill = False,                   # fill style
            )   
        
        # done
        return c

    def section_decor(self):

        l = vlines(
            [-self.r, +self.r],             # horizontal positions
            0.0, 1.0,                       # vertical span
            edgecolor = (0.7, 0.7, 0.7),    # color
            linestyle = "--",               # line style
            linewidth = 0.25,               # line thickness
            )

        # done
        return l
