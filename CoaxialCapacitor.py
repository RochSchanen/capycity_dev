# file: coaxialCapacitor.py
# author: Roch Schanen
# updated: 2022 12 12

# from package: "https://matplotlib.org/"
from matplotlib.pyplot import plot
from matplotlib.pyplot import title
from matplotlib.pyplot import quiver

from matplotlib.patches import Rectangle

# from package: "https://numpy.org/"
from numpy import arange, double, empty, shape
from numpy import square, sqrt, log
from numpy import gradient, meshgrid, average
from numpy import load, save

# from package "capycity" (in development)
from capycity.config   import ZER, MAX
from capycity.solvers  import laplace2DSolver
from capycity.geometry import Disk
from capycity.graphics import opendocument, closedocument
from capycity.graphics import selectfigure, exportfigure
from capycity.graphics import plotcontour, plotmesh

# standard library
from os.path import exists
from os import remove
from sys import exit

#################################################################### PARAMETERS

VERBOSE = False

# force computation (True)
CLEAR = False

# local parameters
r1 = 0.10   # center radius
r2 = 0.40   # shield radius
p0 = 3      # intial resolution 8 X 8

USEPLOT = {
    "MESH":         False,
    "CONTOUR":      True,
    "CROSSSECTION": True,
    "ERROR":        True,
    "GRADIENT":     True,
}

# outputfiles
datapath_solver = "coaxial_solver.npy"
datapath_theory = "coaxial_theory.npy"

######################################################################## SOLVER

# clear files to force computation
if CLEAR:
    if exists(datapath_solver): remove(datapath_solver)
    if exists(datapath_theory): remove(datapath_theory)

# instanciate the solver with boundary masks
solver, p = laplace2DSolver(p0), 0
solver.mergeMask(Disk(r1), "S")     # CENTER
solver.mergeMask(Disk(r2), "C")     # SHELL

if exists(datapath_solver):
    print(f"load previous computation")
    # re-load previous computation
    solver.load(datapath_solver)
    print(f"resolution is {solver.n}X{solver.n}.")
else:
    # compute initial resolution
    solver.jacobiSteps(500) # 8 X 8
    # iterate through increasing resolution
    for n in [100, 100, 100, 1000, 5000]: # 16, 32, 64, 128, 256
        # new graph name
        m, p = f"name_{p}", p+1
        # display info
        print(f"new graph name {m}, {n} jacobi steps ({solver.n}X{solver.n}).")
        # increase resolution, run jacobi steps and plot
        solver.incrementResolution()
        solver.jacobiSteps(n)
    solver.save(datapath_solver)

######################################################################## THEORY

if exists(datapath_theory):
    print(f"load previous theory data")
    # re-load previous computation
    M = load(datapath_theory)
else:
    # compute theoretical distribution
    n, n = shape(solver.D)
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

if USEPLOT["CONTOUR"]:
    # print(f"plot potential contour")
    plotcontour(solver, 8, name = "CONTOUR")
    title(" Potential [V]")

###############################                                            MESH
""" mesh smooth less the data but is too
slow to export for high resolution """

if USEPLOT["MESH"]:
    # print(f"plot potential mesh")
    plotmesh(solver, name = "MESH")
    title(" Potential [V]")

###############################                                   CROSS SECTION
""" plotting a cross section is the most
quantitative form of information display """

if USEPLOT["CROSSSECTION"]:
    # print(f"plot potential cross section")
    n, n = shape(solver.D)
    # center line index
    cline = int((n -  2) // 2)
    # limits and domain
    inter = solver.l / (n - 2)
    limit = (solver.l - inter) / 2.0
    X = arange(-limit, +limit + inter, inter)
    # create figure
    fg, ax = selectfigure("CROSSSECTION")
    # select solver cross section and plot
    Y = solver.D[cline, 1:-1] / MAX
    plot(X, Y, ".", color = "C0")
    # select theory cross section and plot
    T = M[cline, +1:-1]
    plot(X, T, "r-", linewidth = 0.75)
    ax.set_xticks(list(arange(-0.5, 0.5 + 0.1, 0.1)))
    ax.set_yticks(list(arange( 0.0, 1.0 + 0.1, 0.1)))
    ax.grid(True, 'both')
    ax.legend(["Numerical", "Analytical"])
    title(" Potential [V]")

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
    ax.set_title("Error [%]")

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

if USEPLOT["GRADIENT"]:
    # print(f"plot gradient")
    fg, ax = selectfigure("GRADIENT")
    # normalise and remove edges
    D = solver.D[1:-1, 1:-1] / MAX
    # compute gradient components [a.u.]
    dy, dx = gradient(D)
    # compute intervals value [m]
    d = solver.l / n
    # scale to S.I. units [V/m]
    EX, EY = dx / d, dy / d
    # reduced data size by averaging for display clarity
    n, n = shape(EY)
    while n > 16: # choose here the number of arrows
        EX = coarsenresolution(EX)
        EY = coarsenresolution(EY)
        n, n = shape(EY)
    # limits and domain
    d = solver.l / n
    l = (solver.l - d) / 2.0
    D = arange(-l, +l + d, d)
    X, Y = meshgrid(D, D)
    # plot
    quiver(X, Y, EX, EY,
        color="C0",
        # angles='xy',
        # scale_units='xy',
        scale_units='inches',
        scale=0.2*1000/25.4,
        width=.005,
        pivot = 'mid',
        )
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

################################################################# EXPORT TO PDF

opendocument("./local/output.pdf")
for k in USEPLOT.keys():
    if USEPLOT[k]:
        exportfigure(k) 
closedocument()
