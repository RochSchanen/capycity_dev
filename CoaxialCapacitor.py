# file: coaxialCapacitor.py

# from package: "https://matplotlib.org/"
from matplotlib.pyplot import plot
from matplotlib.pyplot import title
from matplotlib.pyplot import quiver

# from package: "https://numpy.org/"
from numpy import arange
from numpy import empty
from numpy import square
from numpy import sqrt
from numpy import log
from numpy import double
from numpy import gradient

# from package "capycity"
from capycity.config   import ZER, MAX
from capycity.solvers  import laplace2DSolver
from capycity.geometry import Disk

from capycity.graphics import opendocument, closedocument
from capycity.graphics import selectfigure, exportfigure
from capycity.graphics import plotcontour, plotmesh

opendocument("results.pdf")

# local parameter
r1 = 0.10
r2 = 0.40

# instanciate solver
solver, p = laplace2DSolver(3), 0   # 8 X 8
solver.mergeMask(Disk(r1), "S")     # CENTER
solver.mergeMask(Disk(r2), "C")     # SHELL

# inital computation
solver.jacobiSteps(100)

# iterate through increasing resolution
for n in [100, 100, 100, 1000, 5000]:
# for n in [100, 100, 100, 100, 100]:
    # new graph name
    m, p = f"name_{p}", p+1
    # display info
    print(f"new graph name {m}, {n} jacobi steps ({solver.n}X{solver.n}).")
    # increase resolution, run jacobi steps and plot
    solver.incrementResolution()
    solver.jacobiSteps(n)

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

# center line index of the cross section
cline = int(solver.n // 2 + 1)
# compute parameters
inter = solver.l / (solver.n - 1)
limit = solver.l / 2
# compute domain
X = arange(-limit, +limit + inter, inter)

# FIGURE                                                      Potential contour

plotcontour(solver, 16, name = "Potential Contour")
title(" Potential [V]")

# FIGURE                                                         Potential mesh

plotmesh(solver, name = "Potential Mesh")
title(" Potential [V]")

# FIGURE                                               Potentials Cross Section
 
fg, ax = selectfigure("Potential Cross Section")
# select solver cross section and plot
Y = solver.D[cline, 1:-1] / MAX
plot(X, Y, "g.")
# select theory cross section and plot
T = M[cline, :]
plot(X, T, "b--")
ax.set_xticks(list(arange(-0.5, 0.5 + 0.1, 0.1)))
ax.set_yticks(list(arange( 0.0, 1.0 + 0.1, 0.1)))
title(" Potential [V]")

# FIGURE                                                             Difference

# compute difference and plot
fg, ax = selectfigure("Difference")
plot(X, 100*(Y-T), "r--")
ax.set_xticks(list(arange(-0.5, 0.5 + 0.1, 0.1)))
# ax.set_yticks(list(arange(-5, 5 + 1, 1)))
title("Error [%]")

# FIGURE

fg, ax = selectfigure("gradient")
dY, dX = gradient(solver.D)
quiver([X, X], dX[1:-1,1:-1], dY[1:-1,1:-1])

# EXPORT FIGURES

# export and done
exportfigure("Potential Mesh")
exportfigure("Potential Contour")
exportfigure("gradient")
exportfigure("Potential Cross Section")
exportfigure("Difference")
closedocument()
