# file: PlateCapacitor.py
# author: Roch Schanen
# updated: 2022 12 12

# from package: "https://matplotlib.org/"
from matplotlib.pyplot import title

# from package: "https://numpy.org/"
from numpy import arange, double, empty, shape, array
from numpy import gradient
from numpy import load, save
from numpy import sum as Sum
from numpy import multiply as Multiply
# from numpy import square, sqrt, log, meshgrid, average

# 
from scipy.constants import epsilon_0 as e0

# from package "capycity" (in development)
from capycity.solvers  import laplace2DSolver
from capycity.geometry import Disk, Aperture, Plate
from capycity.graphics import opendocument, closedocument
from capycity.graphics import selectfigure, exportfigure
# from capycity.graphics import plotcontour, plotmesh
# from capycity.config   import ZER, MAX

# standard library
from os.path import exists
from os import remove
from sys import exit

# ################################################################### CAPACITANCE

def computecapacitance(S1, S2):

    # potential gradients
    S1.computegradient()
    S2.computegradient()

    # electric fields
    EX1, EY1 = S1.dD
    EX2, EY2 = S2.dD

    # "field" dot product
    integrand = Multiply(EX1, EX2) + Multiply(EY1, EY2)

    # get integration intervals
    n, n = shape(EX1)
    d = S1.l / n 

    # numeric integral
    integrale = Sum(integrand)*d*d
    
    # capacitance per unit length (third dimension)
    CN = -integrale*e0

    # done
    return CN

#################################################################### PARAMETERS

# local parameters
r = 0.9 / 2 # shield radius
w = 0.4     # width
h = 0.1     # height
g = 0.05    # gap 

################################################################# configuration

SUBSET = [
    [1]*20+[10]*15, [50]*5, [100]*5, [150]*5,
    [100]*5,  # 512
    [100]*5,  # 1024
    [100]*5,  # 2048
    [100]*5,  # 4096
]

SERIESLIST = [
    ( 0.0 * w/2.0, 5, SUBSET),
    ( 0.1 * w/2.0, 5, SUBSET),
    ( 0.2 * w/2.0, 5, SUBSET),
    ( 0.3 * w/2.0, 5, SUBSET),
    ( 0.4 * w/2.0, 5, SUBSET),
    ( 0.5 * w/2.0, 5, SUBSET),
    ( 0.6 * w/2.0, 5, SUBSET),
    ( 0.7 * w/2.0, 5, SUBSET),
    ( 0.8 * w/2.0, 5, SUBSET),
    ( 0.9 * w/2.0, 5, SUBSET),
    ( 1.0 * w/2.0, 5, SUBSET),
    ( 1.1 * w/2.0, 5, SUBSET),
    ]

pm = 0
fg, ax = selectfigure("FIGURE")

for sh, P0, SERIES in SERIESLIST:
    
    # top electrode potential
    S1 = laplace2DSolver(P0)
    S1.mergeMask(Aperture(r), maskType = "0")
    S1.mergeMask(Plate(-sh, (+g+h)/2, w, h), maskType = "1")
    S1.mergeMask(Plate(+sh, (-g-h)/2, w, h), maskType = "0")

    # bottom electrode potential
    S2 = laplace2DSolver(P0)
    S2.mergeMask(Aperture(r), maskType = "0")
    S2.mergeMask(Plate(0, (+g+h)/2, w, h), maskType = "0")
    S2.mergeMask(Plate(0, (-g-h)/2, w, h), maskType = "1")

    p, P, C = 0, [], []
    
    # step 1
    p += 1
    S1.jacobiSteps(1)
    S2.jacobiSteps(1)
    P.append(p)
    C.append(computecapacitance(S1, S2))

    S = SERIES[0]

    # step blocks
    for n in S:
                    
        p += n
        S1.jacobiSteps(n)
        S2.jacobiSteps(n)
        P.append(p)
        C.append(computecapacitance(S1, S2))

    # resolutions
    for S in SERIES[1:]: 

        # increment resolution
        p += 1
        S1.incrementResolution()
        S2.incrementResolution()    
        P.append(p)
        C.append(computecapacitance(S1, S2))

        # step 1
        p += 1
        S1.jacobiSteps(1)
        S2.jacobiSteps(1)
        P.append(p)
        C.append(computecapacitance(S1, S2))

        # step blocks
        for n in S:
                        
            p += n
            S1.jacobiSteps(n)
            S2.jacobiSteps(n)
            P.append(p)
            C.append(computecapacitance(S1, S2))

    pm = max(pm, p)
    X, Y = array(P), array(C)*1E12
    # ax.plot(X, Y, ".-")
    ax.semilogx(X, Y, ".-")

    print(f"{Y[-1]}")

# Ideal - Analytical value
C = e0* w / g *1E12
ax.hlines([C], 0, pm, colors = "r", linestyle = "dashed")

ax.grid(True, "both")
# ax.legend([])

opendocument("./local/plots.pdf")
exportfigure("FIGURE") 
closedocument()
