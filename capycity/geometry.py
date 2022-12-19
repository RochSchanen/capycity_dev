# file: geometry.py
# author: Roch Schanen
# date: 2022-11-20
# content: definition classes for the geometries used by the solver

# from package: "https://matplotlib.org/"
from matplotlib.pyplot import Circle, Rectangle
# from matplotlib.pyplot import vlines, hlines

# from package: "https://numpy.org/"
from numpy import array, empty, shape, asarray
from numpy import ceil, floor

# from package "https://python-pillow.org/"
from PIL import Image, ImageDraw

# from package "capycity" (in heavy development)
from capycity.config import DATATYPE, MAX

#####################################################################      DISK
###                            DISK                              ####
#####################################################################

class Disk():

    # parameters
    def __init__(self, r):
        # record parameters
        self.l, self.r = -r, +r
        self.t, self.b = +r, -r
        # done
        return

    # mask
    def mask(self, S):
        # get array size (without edges)
        n = S.n
        # compute index origin
        o = (n + 1.0) / 2.0
        # compute interval
        a = S.l / n
        # boundary box
        l, r =  ceil(self.l / a + o), floor(self.r / a + o)
        t, b = floor(self.t / a + o),  ceil(self.b / a + o)
        # create image
        i = Image.new(mode = "L", size = (n+2, n+2), color = 0)
        # link handle "d" to image "i"
        d = ImageDraw.Draw(i)
        # punch aperture
        d.ellipse([l, b, r, t], fill = 1, width = 0)
        # convert to numpy array
        M = asarray(i, DATATYPE) * MAX
        # done
        return M

    # graph
    def contour2D_decor(self):
        c = Circle(
            (0.0, 0.0),                     # origin
            self.r,                         # size
            edgecolor = (0.7, 0.7, 0.7),    # color
            linestyle = "--",               # line style
            fill = False,                   # fill style
            )
        # done
        return c

#####################################################################  APERTURE
###                            APERTURE                          ####
#####################################################################
    
class Aperture():

    # parameters
    def __init__(self, r):
        # record parameters
        self.l, self.r = -r, +r
        self.t, self.b = +r, -r
        # done
        return

    # mask
    def mask(self, S):
        # get array size (without edges)
        n = S.n
        # compute index origin
        o = (n + 1.0) / 2.0
        # compute interval
        a = S.l / n
        # boundary box
        l, r =  ceil(self.l / a + o), floor(self.r / a + o)
        t, b = floor(self.t / a + o),  ceil(self.b / a + o)
        # create image
        i = Image.new(mode = "L", size = (n+2, n+2), color = 1)
        # link handle "d" to image "i"
        d = ImageDraw.Draw(i)
        # punch aperture
        d.ellipse([l, b, r, t], fill = 0, width = 0)
        # convert to numpy array
        M = asarray(i, DATATYPE) * MAX
        # done
        return M

    # graph
    def contour2D_decor(self):
        c = Circle(
            (0.0, 0.0),                     # origin
            self.r,                         # size
            edgecolor = (0.7, 0.7, 0.7),    # color
            linestyle = "--",               # line style
            fill = False,                   # fill style
            )
        # done
        return c


#####################################################################     PLATE
###                            PLATE                             ####
#####################################################################

class Plate():

    # parameters
    def __init__(self, x, y, w, h):
        # record parameters
        self.l, self.r = x - w / 2.0, x + w / 2.0
        self.t, self.b = y + h / 2.0, y - h / 2.0
        # done
        return

    # mask
    def mask(self, S):
        # get array size (without edges)
        n = S.n
        # compute index origin
        o = (n + 1.0) / 2.0
        # compute interval
        a = S.l / n
        # boundary box
        l, r =  ceil(self.l / a + o), floor(self.r / a + o)
        t, b = floor(self.t / a + o),  ceil(self.b / a + o)
        # create image
        i = Image.new(mode = "L", size = (n+2, n+2), color = 0)
        # link handle "d" to image "i"
        d = ImageDraw.Draw(i)
        # punch aperture
        d.rectangle([l, b, r, t], fill = 1, width = 0)
        # convert to numpy array
        M = asarray(i, DATATYPE) * MAX
        # done
        return M

    def contour2D_decor(self):
        
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
