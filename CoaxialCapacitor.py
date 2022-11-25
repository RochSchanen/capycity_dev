# file: coaxialCapacitor.py

# from package: "https://matplotlib.org/"
from matplotlib.backends.backend_pdf import PdfPages

# from package: "https://numpy.org/"
# from numpy import save
# from numpy import load

# from package "capycity"
from capycity.solvers  import laplace2DSolver
from capycity.geometry import Disk
from capycity.graphics import plotPotentialContourFill
from capycity.graphics import plotCrossSection
from capycity.graphics import plotTheory

#####################################################################
###                            SOLVER                            ####
#####################################################################

document = PdfPages("results.pdf")

# instanciate solver
solver, p = laplace2DSolver(3), 0   # 8 X 8
solver.mergeMask(Disk(0.10), "S")   # CENTER
solver.mergeMask(Disk(0.40), "C")   # SHELL

# inital computation
solver.jacobiSteps(100)

# iterate through increasing resolution
for n in [100, 100, 100, 1000, 1000, 1000]:
    # new graph name
    m, p = f"name_{p}", p+1
    # debug
    print(f"new graph name {m}, {n} jacobi steps.")
    # increase resolution, run jacobi steps and plot
    solver.incrementResolution()
    solver.jacobiSteps(n)
    
    plotCrossSection(
        solver,
        "k.", linewidth = 1.0,
        name = m)
    
    plotTheory(
        0.1, 0.4, 0.0, 1.0,
        "r--", linewidth = 1.0,
        name = m, pdfdoc = document)

# the first section takes about 45s

################

# new graph

plotCrossSection(
    solver,
    "k-", linewidth = 1.0,
    name = "compare")

plotTheory(
    0.1, 0.4, 0.0, 1.0,
    "r--", linewidth = 1.0,
    name = "compare")

# new solver
solver, p = laplace2DSolver(solver.p), 0
solver.mergeMask(Disk(0.10), "S")   # CENTER
solver.mergeMask(Disk(0.40), "C")   # SHELL

# inital computation
solver.jacobiSteps(10000)

plotCrossSection(
    solver,
    "g-", linewidth = 1.0,
    name = "compare", pdfdoc = document)

solver.jacobiSteps(20000)

plotCrossSection(
    solver,
    "b-", linewidth = 1.0,
    name = "compare", pdfdoc = document)

solver.jacobiSteps(20000)

plotCrossSection(
    solver,
    "b-", linewidth = 1.0,
    name = "compare", pdfdoc = document)

solver.jacobiSteps(20000)

plotCrossSection(
    solver,
    "b-", linewidth = 1.0,
    name = "compare", pdfdoc = document)

# done
document.close()

# [Finished in 946.5s]
