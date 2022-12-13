# file: coaxialCapacitor.py
# author: Roch Schanen
# updated: 2022 12 12

# from package: "https://matplotlib.org/"
from matplotlib.pyplot import plot
from matplotlib.pyplot import title
from matplotlib.pyplot import quiver

# from package: "https://numpy.org/"
from numpy import arange, double, empty, shape
from numpy import square, sqrt, log
from numpy import gradient, meshgrid
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
    solver.jacobiSteps(100) # 8 X 8
    # iterate through increasing resolution
    for n in [100, 100, 100, 1000, 5000]: # 16, 32, 64, 128, 256
    # for n in [100, 100]: # 16, 32
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
    n = solver.n
    # reserve memory
    M = empty((n, n), dtype = double)
    # get grid origin
    o = (n - 1) / 2.0
    # get scaling factor, length = 1 unit
    a = 1.0 / n
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
    plotcontour(solver, 16, name = "CONTOUR")
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
    # center line index
    cline = int(solver.n // 2 + 1)
    # limits
    inter = solver.l / (solver.n - 1)
    limit = solver.l / 2
    # domain
    X = arange(-limit, +limit + inter, inter)
    # create figure
    fg, ax = selectfigure("CROSSSECTION")
    # select solver cross section and plot
    Y = solver.D[cline, 1:-1] / MAX
    plot(X, Y, ".", color = "C0")
    # select theory cross section and plot
    T = M[cline, :]
    plot(X, T, "r--")
    ax.set_xticks(list(arange(-0.5, 0.5 + 0.1, 0.1)))
    ax.set_yticks(list(arange( 0.0, 1.0 + 0.1, 0.1)))
    title(" Potential [V]")

###############################                                       ERROR [%]
""" typical error from the theoretical
potential distribution. typically on
the order of 1% is acceptable """

if USEPLOT["ERROR"]:
    # print(f"plot Error")
    fg, ax = selectfigure("ERROR")
    plot(X, 100*(Y-T), "r--")
    ax.set_xticks(list(arange(-0.5, 0.5 + 0.1, 0.1)))
    ax.set_yticks(list(arange(-5, 5 + 1, 1)))
    ax.grid(True, 'both')
    ax.set_title("Error [%]")

###############################                                        GRADIENT
""" this figure is useful to check
the numerics consistency """

if USEPLOT["GRADIENT"]:
    fg, ax = selectfigure("GRADIENT")
    dY, dX = gradient(solver.D/MAX)
    # quiver intervals
    p = 2**(max(0, solver.p-4))
    # limits
    inter = solver.l / p
    limit = solver.l / 2
    # domain
    D = arange(-limit, +limit, inter)
    X, Y = meshgrid(D, D)
    quiver(
        X, Y, 
        dX[1:-1:p, 1:-1:p],
        dY[1:-1:p, 1:-1:p],
        color="C0",
        # angles='xy',
        # scale_units='xy',
        scale_units='inches',
        scale=2/25.4,
        width=.005,
        pivot = 'mid',
        )
    ax.set_xticks(list(arange(-0.5, 0.5 + 0.1, 0.1)))
    ax.set_yticks(list(arange(-0.5, 0.5 + 0.1, 0.1)))
    ax.set_xlim(-0.6, +0.6)
    ax.set_ylim(-0.6, +0.6)

################################################################# EXPORT TO PDF

opendocument("plots.pdf")
for k in USEPLOT.keys():
    if USEPLOT[k]:
        exportfigure(k) 
closedocument()
