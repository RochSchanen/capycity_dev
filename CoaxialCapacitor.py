# file: coaxialCapacitor.py
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
from capycity.geometry import Disk, Aperture
from capycity.graphics import opendocument, closedocument
from capycity.graphics import selectfigure, exportfigure
from capycity.graphics import plotcontour, plotmesh

# standard library
from os.path import exists
from os import remove
from sys import exit

#################################################################### PARAMETERS

# VERBOSE = False

# force re-computation when CLEAR = True
# use previous results when CLEAR = False
CLEAR = True

# intial resolution
P0 = 3 # 8 X 8

# computation series 8, 16, 32, 64, 128, 256, 512, 1024
SERIES = [100, 100, 100, 100, 100, 100, 100, 100]

# local parameters
r1 = 0.10   # center radius
r2 = 0.40   # shield radius

USEPLOT = {
    "CROSSSECTION": True,
    "ERROR":        True,

    "MESH 1":       False,
    "CONTOUR 1":    True,
    "GRADIENT 1":   True,

    "MESH 2":       False,
    "CONTOUR 2":    True,
    "GRADIENT 2":   True,
}

# outputfiles
datapath_theory = "coaxial_theory.npy"
datapath_PHI1   = "coaxial_PHI1.npy"
datapath_PHI2   = "coaxial_PHI2.npy"

# clear files to force computation
if CLEAR:
    if exists(datapath_theory): remove(datapath_theory)
    if exists(datapath_PHI1): remove(datapath_PHI1)
    if exists(datapath_PHI2): remove(datapath_PHI2)

######################################################################## PHI1

# instanciate the solver with boundary masks
S1, p = laplace2DSolver(P0), 0
S1.mergeMask(Disk(r1), maskType = "1")
S1.mergeMask(Aperture(r2), maskType = "0")

if exists(datapath_PHI1):
    print(f"load previous computation")
    # re-load previous computation
    S1.load(datapath_PHI1)
    print(f"resolution is {S1.n}X{S1.n}.")
else:
    # iterations at initial resolution
    print(f"step name_{p}, {SERIES[0]} jacobi steps ({S1.n}X{S1.n}).")
    S1.jacobiSteps(SERIES[0]) # 8 X 8
    # iterations through increasing resolution
    for n in SERIES[1:]:
        # new graph name
        m, p = f"name_{p+1}", p+1
        # increase resolution, run jacobi steps and plot
        S1.incrementResolution()
        # display info
        print(f"step {m}, {n} jacobi steps ({S1.n}X{S1.n}).")
        S1.jacobiSteps(n)
    S1.save(datapath_PHI1)

######################################################################## THEORY

if exists(datapath_theory):
    print(f"load previous theory data")
    # re-load previous computation
    M = load(datapath_theory)
else:
    # compute theoretical distribution
    n, n = shape(S1.D)
    # reserve memory
    M = empty((n, n), dtype = double)
    # get grid origin
    o = (n - 1) / 2.0
    # get scaling factor, length = 1 unit
    a = 1.0 / (n - 2)
    # compute theoretical values
    for i in range(n):
        # compute coordinate
        y = (i-o)*a
        for j in range(n):
            # compute coordinate
            x = (j-o)*a
            # compute radius
            r = sqrt(square(x) + square(y))
            # center
            if r < r1:
                M[i, j] = 1.0
                continue
            # shell
            if r > r2:
                M[i, j] = 0.0
                continue
            # vacuum
            M[i, j] = log(r / r2) / log(r1 / r2)
    # save results
    save(datapath_theory, M)

####################################################################### FIGURES

###############################                                         CONTOUR
""" contour is mostly visual and
inherently smoothes out the data """

if USEPLOT[r"CONTOUR 1"]:
    # print(f"plot potential contour")
    plotcontour(S1, 32, name = "CONTOUR 1")
    title(r"Potential $\phi_1$ [V]")

###############################                                            MESH
""" mesh smooth less the data but is too
slow to export for high resolution """

if USEPLOT["MESH 1"]:
    # print(f"plot potential mesh")
    plotmesh(S1, name = "MESH 1")
    title(r"Potential $\phi_1$ [V]")

###############################                                   CROSS SECTION
""" plotting a cross section is the most
quantitative form of information display """

if USEPLOT["CROSSSECTION"]:
    # print(f"plot potential cross section")
    n, n = shape(S1.D)
    # center line index
    cline = int((n -  2) // 2)
    # limits and domain
    inter = S1.l / (n - 2)
    limit = (S1.l - inter) / 2.0
    X = arange(-limit, +limit + inter, inter)
    # create figure
    fg, ax = selectfigure("CROSSSECTION")
    # select solver cross section and plot
    Y = S1.D[cline, 1:-1] / MAX
    plot(X, Y, ".", color = "C0")
    # select theory cross section and plot
    T = M[cline, +1:-1]
    plot(X, T, "r-", linewidth = 0.75)
    ax.set_xticks(list(arange(-0.5, 0.5 + 0.1, 0.1)))
    ax.set_yticks(list(arange( 0.0, 1.0 + 0.1, 0.1)))
    ax.grid(True, 'both')
    ax.legend(["Numerical", "Analytical"])
    title(r"Potential $\phi_1$ [V]")

###############################                                       ERROR [%]
""" typical error from the theoretical
potential distribution. typically on
the order of 1% is acceptable """

if USEPLOT["ERROR"]:
    # print(f"plot Error")
    fg, ax = selectfigure("ERROR")
    plot(X, 100*(Y-T), "r.--")
    ax.set_xticks(list(arange(-0.5, 0.5 + 0.1, 0.1)))
    ax.set_yticks(list(arange(-5, 5 + 1, 1)))
    ax.grid(True, 'both')
    ax.set_title(r"Error $\phi_1$ [%]")

###############################                                        GRADIENT
""" this figure is useful to check
consistency of the results """

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

if USEPLOT["GRADIENT 1"]:
    # print(f"plot gradient")
    fg, ax = selectfigure("GRADIENT 1")
    # normalise and remove edges
    D = S1.D[1:-1, 1:-1] / MAX
    # get size
    n, n = shape(D)
    # compute gradient components [a.u.]
    dy, dx = gradient(D)
    # compute intervals value [m]
    d = S1.l / n
    # scale to S.I. units [V/m]
    EX1, EY1 = dx / d, dy / d
    # reduced data size by averaging for display clarity
    ex, ey = EX1, EY1
    n, n = shape(ex)
    while n > 16: # choose here the number of arrows
        ex = coarsenresolution(ex)
        ey = coarsenresolution(ey)
        n, n = shape(ex)
    # limits and domain
    d = S1.l / n
    l = (S1.l - d) / 2.0
    D = arange(-l, +l + d, d)
    X, Y = meshgrid(D, D)
    # plot
    quiver(X, Y, ex, ey,
        color="C0",
        # angles='xy',
        # scale_units='xy',
        scale_units='inches',
        scale=0.2*1000/25.4,
        # scale=0.1*1000/25.4,
        width=.005,
        pivot = 'mid',
        )
    # plot
    scatter(X, Y, 5.0, c = 'k')
    # cosmetics
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
            # fill        = False,
            # linewidth   = 0.5,
            zorder      = 0,
            ))
    # add boundaries
    for g in S1.G:
        ax.add_artist(g.contour2D_decor())
    # add title
    title(r"Averaged Electric field $E_1$ [V/m]")

# ######################################################################## PHI2

# instanciate the solver with boundary masks
S2, p = laplace2DSolver(P0), 0
S2.mergeMask(Disk(r1), maskType = "0")
S2.mergeMask(Aperture(r2), maskType = "1")

if exists(datapath_PHI2):
    print(f"load previous computation")
    # re-load previous computation
    S2.load(datapath_PHI2)
    print(f"resolution is {S2.n}X{S2.n}.")
else:
    # iterations at initial resolution
    print(f"step name_{p}, {SERIES[0]} jacobi steps ({S2.n}X{S2.n}).")
    S2.jacobiSteps(SERIES[0]) # 8 X 8
    # iterations through increasing resolution
    for n in SERIES[1:]:
        # new graph name
        m, p = f"name_{p+1}", p+1
        # increase resolution, run jacobi steps and plot
        S2.incrementResolution()
        # display info
        print(f"step {m}, {n} jacobi steps ({S2.n}X{S2.n}).")
        S2.jacobiSteps(n)

    S2.save(datapath_PHI2)

####################################################################### FIGURES

###############################                                         CONTOUR
""" contour is mostly visual and
inherently smoothes out the data """

if USEPLOT[r"CONTOUR 2"]:
    # print(f"plot potential contour")
    plotcontour(S2, 16, name = "CONTOUR 2")
    title(r"Potential $\phi_2$ [V]")

###############################                                            MESH
""" mesh smooth less the data but is too
slow to export for high resolution """

if USEPLOT["MESH 2"]:
    # print(f"plot potential mesh")
    plotmesh(S2, name = "MESH 2")
    title(r"Potential $\phi_2$ [V]")

###############################                                        GRADIENT
""" this figure is useful to check
consistency of the results """

if USEPLOT["GRADIENT 2"]:
    # print(f"plot gradient")
    fg, ax = selectfigure("GRADIENT 2")
    # normalise and remove edges
    D = S2.D[1:-1, 1:-1] / MAX
    # get size
    n, n = shape(D)
    # compute gradient components [a.u.]
    dy, dx = gradient(D)
    # compute intervals value [m]
    d = S2.l / n
    # scale to S.I. units [V/m]
    EX2, EY2 = dx / d, dy / d
    # reduced data size by averaging for display clarity
    ex, ey = EX2, EY2
    n, n = shape(ex)
    while n > 16: # choose here the number of arrows
        ex = coarsenresolution(ex)
        ey = coarsenresolution(ey)
        n, n = shape(ex)
    # limits and domain
    d = S2.l / n
    l = (S2.l - d) / 2.0
    D = arange(-l, +l + d, d)
    X, Y = meshgrid(D, D)
    # plot
    quiver(X, Y, ex, ey,
        color="C0",
        # angles='xy',
        # scale_units='xy',
        scale_units='inches',
        scale=0.2*1000/25.4,
        # scale=0.1*1000/25.4,
        width=.005,
        pivot = 'mid',
        )
    # plot
    scatter(X, Y, 5.0, c = 'k')
    # cosmetics
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
            # fill        = False,
            # linewidth   = 0.5,
            zorder      = 0,
            ))
    # add boundaries
    for g in S2.G:
        ax.add_artist(g.contour2D_decor())
    # add title
    title(r"Averaged Electric field $E_2$ [V/m]")

################################################################### CAPACITANCE

e0  = 8.854187817E-12

from numpy import pi

integrand = Multiply(EX1, EX2) + Multiply(EY1, EY2)

# get interval
n, n = shape(EX2)
d = S2.l / n 

# compute integrale numerically
integrale = Sum(integrand)*d*d
CN = integrale*e0

print(f"integrale C/l = {CN:.3e} F/m")

theory = -2*pi/log(r2/r1)
CT = theory*e0

print(f"theory C/l = {CT:.3e} F/m")

# print error [%]
print(f"error = {(CN-CT)/CT*100:.2f} %")

################################################################# EXPORT TO PDF

opendocument("./local/output.pdf")
for k in USEPLOT.keys():
    if USEPLOT[k]:
        exportfigure(k) 
closedocument()
