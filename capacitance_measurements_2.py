# file: capacitance_measurements_2.py
# created: 2023 01 07
# author: Roch Schanen

# libraries, function and classes

from numpy import pi, array, linspace
from matplotlib.pyplot import figure, fignum_exists
from matplotlib.backends.backend_pdf import PdfPages

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
        self._DOC = PdfPages(pathname)
        return self._DOC

    def exportfigure(self, name):
        args = selectfigure(name)
        self._DOC.savefig(args[0])
        return

    def closedocument(self):
        self._DOC.close()
        return

#####################################################################
###                          DATA                                 ###
#####################################################################

# overlap data
# COL#1: relative angle in degrees (no absolute zero)
# COL#2: measured capacitance using the ET4502 LCR

D1 = array([ 

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

# gap data (one of the two electrodes)
# COL#1: 1/12 turn units (12 means one turn, pitch is 350µm/turn)
# COL#2: measured capacitance using the ET4502 LCR

D2 = array([

    # [ 0.0,   inf ],
    [ 1.0,  66.86 ],
    [ 2.0,  40.68 ],
    [ 3.0,  30.931],
    [ 4.0,  25.475],
    [ 5.0,  22.116],

    [ 6.0,  19.338],
    [ 7.0,  17.261],
    [ 8.0,  16.427],
    [ 9.0,  14.375],
    [10.0,  13.372],
    [11.0,  12.476],

    [12.0,  11.723],
    [13.0,  11.121],
    [14.0,  10.729],
    [15.0,  10.270],
    [16.0,   9.842],
    [17.0,   9.524],

    [18.0,   9.203],
])

# gap data (the other of the two electrodes)
# COL#1: 1/12 turn units (12 means one turn, pitch is 350µm/turn)
# COL#2: measured capacitance using the ET4502 LCR

D3 = array([

    # [ 0.0,   inf ],
    [ 1.0,  96.82 ],
    [ 2.0,  47.80 ],
    [ 3.0,  34.19 ],
    [ 4.0,  26.446],
    [ 5.0,  21.818],

    [ 6.0,  18.868],
    [ 7.0,  16.723],
    [ 8.0,  15.262],
    [ 9.0,  14.531],
    [10.0,  13.481],
    [11.0,  12.721],

    [12.0,  12.001],
    [13.0,  11.356],
    [14.0,  10.810],
    [15.0,  10.365],
    [16.0,   9.909],
    [17.0,   9.510],

    [18.0,   9.092],
])

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

#####################################################################
###                          FIGURES                              ###
#####################################################################

# figure overlap dependence
fg, ax = selectfigure("Angle dependence")

ax.plot(D1[:, 0], D1[:, 1], "r.-")
ax.plot(D1[:, 0]+30, D1[:, 1], "b.-")

# add vertical lines
ax.vlines([list(array([-12, -6, 0, +6, +12])+6.40)], 0, 200, "k", "dashed")

# set decorum
ax.set_ylim(0, 200)
ax.set_xlabel("angle [degrees]")
ax.set_ylabel("Capacitance [pF]")
ax.grid("True")

# theoretical capacitance value as a function of the gap:
def C(d):
    e0 = 8.854E-12
    N, r, R, eta_deg = 12, 7.5E-3, 20E-3, 6
    L, l, eta = R-r, (R+r)/2, pi/180*eta_deg
    C = N*e0*L*l*eta/d
    return C

# figure gap dependence
fg, ax = selectfigure("Gap dependence")

X, Y = D2[:, 0]*350/12, D2[:, 1]
ax.plot(X, Y, "b.-.", linewidth = 0.50)

X, Y = D3[:, 0]*350/12, D3[:, 1]
ax.plot(X, Y, "r.-.", linewidth = 0.50)

# add theoretical values
X = linspace(20.0, 500.0, 50)
ax.plot(X, C(X*1E-6)*1E12, "k-.", linewidth = 0.50)
ax.plot(X, capacitance(array([6.0]), X*1E-6), "k-", linewidth = 0.50)

# set decorum
ax.set_xlim(0, 600)
ax.set_ylim(0, 100)
ax.set_xlabel("gap [µm]")
ax.set_ylabel("Capacitance [pF]")
ax.grid("True")

# export document
doc = Document("../result.pdf")
doc.exportfigure("Angle dependence")
doc.exportfigure("Gap dependence")
doc.closedocument()
