# file: PlateCapacitorGap.py
# author: Roch Schanen
# updated: 2022 12 20

# from package: "https://matplotlib.org/"
from matplotlib.pyplot import title

# from package: "https://numpy.org/"
from numpy import shape, array
from numpy import gradient
from numpy import load, save
from numpy import sum as Sum
from numpy import multiply as Multiply

# from package: "https://scipy.org/"
from scipy.constants import epsilon_0 as e0

# from package "capycity" (in development)
from capycity.solvers  import laplace2DSolver
from capycity.geometry import Disk, Aperture, Plate
from capycity.graphics import opendocument, closedocument
from capycity.graphics import selectfigure, exportfigure

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
    d = S1.l / S1.n
    # numeric integral
    integrale = Sum(integrand)*d*d
    # capacitance per unit length (infinite third dimension)
    CN = -integrale*e0
    # done
    return CN

#################################################################### PARAMETERS

# local parameters
r = 0.9 / 2 # shield radius
w = 0.5     # width
h = 0.1     # height
g = 0.1     # gap 

################################################################# configuration

SERIESLIST = [
    ( 3, [[  1]* 10]*7+[[ 10]*30]),
    ( 3, [[  3]* 10]*7+[[ 10]*30]),
    ( 3, [[ 10]* 10]*7+[[ 10]*30]),
    ( 4, [[100]*100]*1),
    ( 5, [[100]*100]*1),
    ( 6, [[100]*100]*1),
    ( 7, [[100]*100]*1),
    ( 8, [[100]*100]*1),
    ( 9, [[100]*100]*1),
    ( 10, [[100]*100]*1),
    ]

pm = 0
fg, ax = selectfigure("FIGURE")

for P0, SERIES in SERIESLIST:
    
    # top electrode potential
    S1 = laplace2DSolver(P0)
    S1.mergeMask(Aperture(r), maskType = "0")
    S1.mergeMask(Plate(0, (+g+h)/2, w, h), maskType = "1")
    S1.mergeMask(Plate(0, (-g-h)/2, w, h), maskType = "0")

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

# Ideal - Analytical value
C = e0* w / g *1E12
ax.hlines([C], 0, pm, colors = "r", linestyle = "dashed")

ax.grid(True, "both")
# ax.legend([])

opendocument("./local/plots.pdf")
exportfigure("FIGURE") 
closedocument()
