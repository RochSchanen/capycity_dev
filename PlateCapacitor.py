# file: PlateCapacitor.py
# author: Roch Schanen
# updated: 2022 12 12

# from package: "https://matplotlib.org/"
from matplotlib.pyplot import plot, quiver, scatter
from matplotlib.pyplot import title

from matplotlib.patches import Rectangle

# from package: "https://numpy.org/"
from numpy import arange, double, empty, shape
from numpy import square, sqrt, log
from numpy import gradient, meshgrid, average
from numpy import load, save
from numpy import sum as Sum
from numpy import multiply as Multiply

# from package "capycity" (in development)
from capycity.config   import ZER, MAX
from capycity.solvers  import laplace2DSolver
from capycity.geometry import Disk, Aperture, Plate
from capycity.graphics import opendocument, closedocument
from capycity.graphics import selectfigure, exportfigure
from capycity.graphics import plotcontour, plotmesh

# standard library
from os.path import exists
from os import remove
from sys import exit

#################################################################### PARAMETERS

# file paths
FP0 = "plate_PHI0.npy"
FP1 = "plate_PHI1.npy"
FP2 = "plate_PHI2.npy"

# use previous computations
USEFILES = {
    FP0: True,
    FP1: True,
    FP2: True,
}

# clear unused files
for fp in USEFILES.keys():
    if not USEFILES[fp]:
        if exists(fp):
            remove(fp)

# intial resolution
P0 = 3 # 8 X 8

# computation series 8, 16, 32, 64, 128, 256, 512, 1024
SERIES = [500, 500, 500, 500, 1000, 2000, 5000, 10000]

# computation series 8, 16, 32, 64, 128
# SERIES = [100, 100, 100, 100, 100]

# computation series 8, 16
# SERIES = [100, 100]

# local parameters
r = 0.9 / 2 # shield radius
w = 0.5     # width
h = 0.1     # height
g = 0.1     # gap 

MAKEPLOTS = {
    "MESH0":        False,
    "CONTOUR0":     True,
    "GRADIENT0":    True,
    "CONTOUR1":     True,
    "GRADIENT1":    True,
    "CONTOUR2":     True,
    "GRADIENT2":    True,
}

####################################################################### COMPUTE

def computefields(S, F):
    """ S is the solver, F is the data file path """

    # compute or load potential field
    if exists(F):
        print(f"load {F}")
        S.load(F)
        print(f"resolution is {S.n}X{S.n}.")

    else:
        m, n, p = f"STEP_0", SERIES[0], 0
        print(f"{m}: {n} jacobi steps for resolution ({S.n}X{S.n}).")
        S.jacobiSteps(SERIES[0])
        for n in SERIES[1:]:
            m, p = f"STEP_{p+1}", p+1
            S.incrementResolution()
            print(f"{m}: {n} jacobi steps for resolution ({S.n}X{S.n}).")
            S.jacobiSteps(n)
        S.save(F)

    # compute electric field
    S.computegradient()

#                                                                          PHI0
S0 = laplace2DSolver(P0)
S0.mergeMask(Aperture(r), maskType = "1")
S0.mergeMask(Plate(0, (+g+h)/2, w, h), maskType = "0")
S0.mergeMask(Plate(0, (-g-h)/2, w, h), maskType = "0")
computefields(S0, FP0)

#                                                                          PHI1
S1 = laplace2DSolver(P0)
S1.mergeMask(Aperture(r), maskType = "0")
S1.mergeMask(Plate(0, (+g+h)/2, w, h), maskType = "1")
S1.mergeMask(Plate(0, (-g-h)/2, w, h), maskType = "0")
computefields(S1, FP1)

#                                                                          PHI2
S2 = laplace2DSolver(P0)
S2.mergeMask(Aperture(r), maskType = "0")
S2.mergeMask(Plate(0, (+g+h)/2, w, h), maskType = "0")
S2.mergeMask(Plate(0, (-g-h)/2, w, h), maskType = "1")
computefields(S2, FP2)

####################################################################### FIGURES

#                                                                         MESH0
if MAKEPLOTS["MESH0"]:
    print(f"make plot MESH0")
    plotmesh(S0, name = "MESH0")
    title(r"$\phi_0$ [V]")
#                                                                      CONTOUR0
if MAKEPLOTS["CONTOUR0"]:
    print(f"make plot CONTOUR0")
    plotcontour(S0, 32, name = "CONTOUR0")
    title(r"$\phi_0$ [V]")

#                                                                      CONTOUR1
if MAKEPLOTS["CONTOUR1"]:
    print(f"make plot CONTOUR1")
    plotcontour(S1, 32, name = "CONTOUR1")
    title(r"$\phi_1$ [V]")

#                                                                      CONTOUR2
if MAKEPLOTS["CONTOUR2"]:
    print(f"make plot CONTOUR2")
    plotcontour(S2, 32, name = "CONTOUR2")
    title(r"$\phi_2$ [V]")
#                                                                     GRADIENTS
def coarsenresolution(X):
    """ X is a floating point square array of dimension N X N
    the array size is halved along both dimension (interlaced)
    and the four output submatrices are averaged into one
    output matrix of dimension N/2 X N/2. Of course N must
    be an even value for the operation to be defined """ 
    # copy submatrices
    NW, NE = X[0::2, 0::2], X[0::2, 1::2]
    SW, SE = X[1::2, 0::2], X[1::2, 1::2]
    # return average
    return (NW + NE + SW + SE) / 4.0

def plotgradient(S, index):
    print(f"make plot gradient: GRADIENT{index}")
    fg, ax = selectfigure(f"GRADIENT{index}")
    # load electric field
    ex, ey = S.dD
    # reduced grid size by averaging for better display
    n, n = shape(ex)
    # loop
    while n > 32: # choose here the grid size
        ex = coarsenresolution(ex)
        ey = coarsenresolution(ey)
        n, n = shape(ex)
    # new mesh
    d = S.l / n               # interval
    l = (S.l - d) / 2.0       # limit
    D = arange(-l, +l + d, d) # domain
    X, Y = meshgrid(D, D)
    # plot
    quiver(X, Y, ex, ey,
        color="C0",
        scale_units='inches',
        scale=0.5*1000/25.4,
        width=.005,
        pivot = 'mid',
        )
    # more cosmetics
    scatter(X, Y, 5.0, c = 'k')
    ax.set_xticks(list(arange(-0.5, 0.5 + 0.1, 0.1)))
    ax.set_yticks(list(arange(-0.5, 0.5 + 0.1, 0.1)))
    ax.set_xlim(-0.6, +0.6)
    ax.set_ylim(-0.6, +0.6)
    ax.grid(True)
    # display volume element at point (0, 0) and other corners
    for i, j in [(0, 0), (0, -1), (-1, 0), (-1, -1)]:
        x, y = X[i, j], Y[i, j]
        l, r, t, b = x-d/2, x+d/2, y+d/2, y-d/2
        ax.add_patch(Rectangle((l, b), d, d,
            edgecolor   = None,
            facecolor   = "gainsboro",
            zorder      = 0,
            ))
    # add boundaries
    for g in S.G:
        ax.add_artist(g.contour2D_decor())
    # add title
    title(f"Averaged Electric field $E_{index}$ [V/m]")

#                                                                     GRADIENT0
if MAKEPLOTS["GRADIENT0"]: plotgradient(S0, 0)
#                                                                     GRADIENT1
if MAKEPLOTS["GRADIENT1"]: plotgradient(S1, 1)
#                                                                     GRADIENT2
if MAKEPLOTS["GRADIENT2"]: plotgradient(S2, 2)

# ################################################################### CAPACITANCE

# e0  = 8.854187817E-12

# from numpy import pi

# integrand = Multiply(EX1, EX2) + Multiply(EY1, EY2)

# # get interval
# n, n = shape(EX2)
# d = S2.l / n 

# # compute integrale numerically
# integrale = Sum(integrand)*d*d
# CN = integrale*e0

# print(f"integrale C/l = {CN:.3e} F/m")

# theory = -2*pi/log(r2/r1)
# CT = theory*e0

# print(f"theory C/l = {CT:.3e} F/m")

# # print error [%]
# print(f"error = {(CN-CT)/CT*100:.2f} %")

# ################################################################# EXPORT TO PDF

opendocument("./local/output.pdf")
for k in MAKEPLOTS.keys():
    if MAKEPLOTS[k]:
        exportfigure(k) 
closedocument()

# Run Laplace 2D Solver
# STEP_0: 500 jacobi steps for resolution (8X8).
# STEP_1: 500 jacobi steps for resolution (16X16).
# STEP_2: 500 jacobi steps for resolution (32X32).
# STEP_3: 500 jacobi steps for resolution (64X64).
# STEP_4: 1000 jacobi steps for resolution (128X128).
# STEP_5: 2000 jacobi steps for resolution (256X256).
# STEP_6: 5000 jacobi steps for resolution (512X512).
# STEP_7: 10000 jacobi steps for resolution (1024X1024).
# Run Laplace 2D Solver
# STEP_0: 500 jacobi steps for resolution (8X8).
# STEP_1: 500 jacobi steps for resolution (16X16).
# STEP_2: 500 jacobi steps for resolution (32X32).
# STEP_3: 500 jacobi steps for resolution (64X64).
# STEP_4: 1000 jacobi steps for resolution (128X128).
# STEP_5: 2000 jacobi steps for resolution (256X256).
# STEP_6: 5000 jacobi steps for resolution (512X512).
# STEP_7: 10000 jacobi steps for resolution (1024X1024).
# Run Laplace 2D Solver
# STEP_0: 500 jacobi steps for resolution (8X8).
# STEP_1: 500 jacobi steps for resolution (16X16).
# STEP_2: 500 jacobi steps for resolution (32X32).
# STEP_3: 500 jacobi steps for resolution (64X64).
# STEP_4: 1000 jacobi steps for resolution (128X128).
# STEP_5: 2000 jacobi steps for resolution (256X256).
# STEP_6: 5000 jacobi steps for resolution (512X512).
# STEP_7: 10000 jacobi steps for resolution (1024X1024).
# make plot CONTOUR0
# make plot CONTOUR1
# make plot CONTOUR2
# export 'CONTOUR0'
# export 'CONTOUR1'
# export 'CONTOUR2'
# [Finished in 1112.2s]
