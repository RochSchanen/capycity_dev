# file: capacitance_measurement.py

#######################################################
### DATA: MEASURED CAPACITANCE VERSUS OVERLAP ANGLE ###
#######################################################
# (the angle origin is arbitrary in these measurments)

from numpy import array

D = array([

    [-16.00,  79.56],
    [-15.75,  76.20],
    [-15.50,  74.70],
    [-15.25,  71.78],
    [-15.00,  69.37],
    [-14.75,  66.71],
    [-14.50,  63.67],
    [-14.25,  61.15],

    [-14.00,  56.96],
    [-13.75,  54.35],
    [-13.50,  52.16],
    [-13.25,  48.77],
    [-13.00,  45.94],
    [-12.75,  43.11],
    [-12.50,  39.922],
    [-12.25,  36.692],

    [-12.00,  34.103],
    [-11.75,  32.888],
    [-11.50,  31.598],
    [-11.25,  30.843],
    [-11.00,  30.528],
    [-10.75,  30.306],
    [-10.50,  29.998],
    [-10.25,  29.870],
    [-10.00,  29.680],

    [ -9.75,  29.629],
    [ -9.50,  29.540],
    [ -9.25,  29.509],
    [ -9.00,  29.450],
    [ -8.75,  29.440],
    [ -8.50,  29.424],
    [ -8.25,  29.460],
    [ -8.00,  29.545],

    [ -7.75,  29.615],
    [ -7.50,  29.735],
    [ -7.25,  30.026],
    [ -7.00,  30.158],
    [ -6.75,  30.474],
    [ -6.50,  30.835],
    [ -6.25,  31.290],
    [ -6.00,  32.103],

    [ -5.75,  33.126],
    [ -5.50,  35.055],
    [ -5.25,  38.615],
    [ -5.00,  41.21],
    [ -4.75,  43.32],
    [ -4.50,  47.01],
    [ -4.25,  50.74],
    [ -4.00,  53.89],

    [ -3.75,  56.74],
    [ -3.50,  59.66],
    [ -3.25,  63.26],
    [ -3.00,  66.56],
    [ -2.75,  68.95],
    [ -2.50,  72.08],
    [ -2.25,  75.65],
    [ -2.00,  78.01],

    [ -1.75,  80.51],
    [ -1.50,  83.96],
    [ -1.25,  87.82],
    [ -1.00,  89.66],
    [ -0.75,  92.38],
    [ -0.50,  95.69],
    [ -0.25,  98.61],
    [  0.00, 101.94],

    [  0.25, 104.92],
    [  0.50, 108.43],
    [  0.75, 110.97],
    [  1.00, 113.75],
    [  1.25, 116.92],
    [  1.50, 120.03],
    [  1.75, 123.81],
    [  2.00, 125.29],

    [  2.25, 128.98],
    [  2.50, 130.97],
    [  2.75, 134.32],
    [  3.00, 137.58],
    [  3.25, 140.30],
    [  3.50, 142.75],
    [  3.75, 145.16],
    [  4.00, 148.41],

    [  4.25, 150.85],
    [  4.50, 153.69],
    [  4.75, 156.82],
    [  5.00, 159.14],
    [  5.25, 161.95],
    [  5.50, 163.84],
    [  5.75, 167.34],
    [  6.00, 168.66],

    [  6.25, 170.43],
    [  6.50, 170.62],
    [  6.75, 169.21],
    [  7.00, 167.03],
    [  7.25, 164.24],
    [  7.50, 161.37],
    [  7.75, 158.47],
    [  8.00, 155.18],

    [  8.25, 153.15],
    [  8.50, 149.97],
    [  8.75, 148.09],
    [  9.00, 145.16],
    [  9.25, 142.75],
    [  9.50, 139.91],
    [  9.75, 136.98],
    [ 10.00, 134.57],

    [ 10.25, 132.80],
    [ 10.50, 129.05],
    [ 10.75, 126.20],
    [ 11.00, 123.46],
    [ 11.25, 120.73],
    [ 11.50, 116.65],
    [ 11.75, 113.27],
    [ 12.00, 110.90],

    [ 12.25, 108.74],
    [ 12.50, 105.45],
    [ 12.75, 101.71],
    [ 13.00,  99.22],
    [ 13.25,  96.02],
    [ 13.50,  91.83],
    [ 13.75,  87.92],
    [ 14.00,  85.42],

    [ 14.25,  83.40],
    [ 14.50,  80.30],
    [ 14.75,  77.35],
    [ 15.00,  75.97],
    [ 15.25,  72.07],
    [ 15.50,  68.01],
    [ 15.75,  65.75],
    [ 16.00,  62.97],
])

####################
### USEFUL TOOLS ###
####################

def selectfigure(name):
    from matplotlib.pyplot import figure, fignum_exists
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

class Document():

    def __init__(self, pathname = None):
        if pathname is not None:
            self._DOC = self.opendocument(pathname)
        return

    def opendocument(self, pathname):
        from matplotlib.backends.backend_pdf import PdfPages
        self._DOC = PdfPages(pathname)
        return self._DOC

    def exportfigure(self, name):
        args = selectfigure(name)
        self._DOC.savefig(args[0])
        return

    def closedocument(self):
        self._DOC.close()
        return

#########################
### IDEAL CAPACITANCE ###
#########################

# the following function compute the overlap angle
# from the angle of rotation (the origin is defined
# when the overlap is maximum)
def overlap(a):
    # a: angle in degrees (origin at maximum)
    
    from numpy import mod, fabs
    # mod(x,1.0) returns the fractional part of x
    # fabs(x) returns the absolute value of x

    # define boundaries, M: maximum overlap, T: period
    M, T = 12.0, 30.0
    
    # compute overlap from angle
    O = M-T*abs(mod(a/T+0.5, 1.0)-0.5)
    
    # nullify negative values
    O[O<0] = 0
    
    # done
    return O 

def capacitance(angle, gap):
    # angle: angle in degrees
    # gap: gap in metres

    from scipy.constants import epsilon_0 as E0
    from numpy import pi, square, fabs

    # overlap "eta" in radiants
    eta = overlap(angle)*pi/180.0

    # define star shape capacitance parameters
    R1, R2 = 40E-3/2, 15E-3/2
    N, L, l = 12, R1-R2, (R1+R2)/2.0

    # compute star shape capacitance
    C0 = N*E0*L*l/gap*eta

    # define ring capacitance parameters
    R3 = 13E-3/2

    # compute ring capacitance
    C1 = E0*(pi*square(R2)-pi*square(R3))/gap

    # done
    return (C0+C1)*1E12 # [pF]

########################
### PDF PLOT RESULTS ###
########################

A = D[:, 0] # angle in degrees
C = D[:, 1] # measured capacitance

# 40µm plastic sheet parameters:
# e: measured dielectric constant relative to vacuum (epsilon)
# g: measured thickness (capacitance gap)
e, g = 1.81, 40E-6

fg, ax = selectfigure("capa-vs-angle")

# plot C versus A (shifted by 6.5 degree to set the origin)
ax.plot(A-6.5, C/e, "k.-")

# same plot shifed by "nearly" one period (30.5 degerees)
ax.plot(A-6.5 + 30.5, C/e, "k.-")

# analytical plot
ax.plot(A, capacitance(A, g), "r-")
ax.plot(A + 30.5, capacitance(A, g), "r-")

# add labels
ax.set_xlabel("angle [degrees]")
ax.set_ylabel("Capacitance [pF]")

# add vertical and horizontal guide lines
ax.vlines(  [list(array([-12, -6, 0, +6, +12]))],
            -10, 110, "b", "dashed", linewidth = 0.50)
ax.vlines(  [list(array([-15, +15]))],
            -10, 110, "r", "dashed", linewidth = 0.50)
# ax.hlines( [list(array([0]))],
#             -30, 50, "b", "dashed", linewidth = 0.50)

# modify ticks
from numpy import linspace
# ax.set_xticks(linspace(-24, 48, 1+(48+24)//12))
# ax.set_xticks(linspace(-24, 48, 1+(48+24)//6), minor = True)
ax.set_yticks(linspace(0, 100, 1+100//20))
ax.set_yticks(linspace(0, 100, 1+100//10), minor = True)

# show grid on both axes
ax.grid("True", axis = "y", which = "both")

# export figure to pdf document
doc = Document("../result.pdf") # file path + name
doc.exportfigure("capa-vs-angle") # figure name
doc.closedocument()

