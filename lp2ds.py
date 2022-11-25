# file: c

# from package: "https://matplotlib.org/"
from matplotlib.backends.backend_pdf import PdfPages

# from package: "https://numpy.org/"
# from numpy import save
# from numpy import load

# from package "capycity"
from capycity.geometry import Disk
from capycity.solvers import laplace2DSolver
from capycity.graphics import plotPotentialContourFill
from capycity.graphics import plotCrossSection
from capycity.graphics import plotTheory

#####################################################################
###                            SOLVER                            ####
#####################################################################


p = PdfPages("results.pdf")

# instanciate solver
l = laplace2DSolver(3) # 8 X 8

l.mergeMask(Disk(0.10), "S") # CENTER 0.2 diameter
l.mergeMask(Disk(0.40), "C") # SHELL 0.8 diameter 

# apply boundaries and display
l.applyMasks()

plotCrossSection(l, "r-", name = "8")
l.jacobiSteps(100) # 8 X 8
plotTheory(0.1, 0.4, 1.0, 0.0, "b--", name = "8")
plotCrossSection(l, "k.--", name = "8", pdfdoc = p)

l.incrementResolution()
plotCrossSection(l, "r-", name = "16") 
l.jacobiSteps(100) # 16 X 16
plotTheory(0.1, 0.4, 1.0, 0.0, "b--", name = "16")
plotCrossSection(l, "k.--", name = "16", pdfdoc = p) 

l.incrementResolution()
plotCrossSection(l, "r-", name = "32")
l.jacobiSteps(100) # 32 X 32
plotTheory(0.1, 0.4, 1.0, 0.0, "b--", name = "32")
plotCrossSection(l, "k.--", name = "32", pdfdoc = p)

l.incrementResolution()
plotCrossSection(l, "r-", name = "64")
l.jacobiSteps(1000) # 64 X 64
plotTheory(0.1, 0.4, 1.0, 0.0, "b--", name = "64")
plotCrossSection(l, "k.--", name = "64", pdfdoc = p)

l.incrementResolution()
plotCrossSection(l, "r-", name = "128")
l.jacobiSteps(1000) # 128 X 128
plotTheory(0.1, 0.4, 1.0, 0.0, "b--", name = "128")
plotCrossSection(l, "k.--", name = "128", pdfdoc = p)

###############################################################################

l.incrementResolution()
# plotCrossSection(l, "r-", name = "256")
l.jacobiSteps(2000) # 256 X 256
# plotCrossSection(l, "k-", linewidth = 0.35, name = "256", pdfdoc = p)
plotCrossSection(l, "k-", linewidth = 1.0, name = "256")
plotTheory(0.1, 0.4, 1.0, 0.0, "b--", linewidth = 1.0, name = "256", pdfdoc = p)

# l.incrementResolution()
# # plotCrossSection(l, "r-", name = "512")
# l.jacobiSteps(500) # 512 X 512
# # plotTheory(0.1, 0.4, 1.0, 0.0, "b--", linewidth = 1.0, name = "256")
# plotCrossSection(l, "r-", linewidth = 0.66, name = "256")

# l.incrementResolution()
# # plotCrossSection(l, "r-", name = "512")
# l.jacobiSteps(500) # 1024 X 1024
# # plotTheory(0.1, 0.4, 1.0, 0.0, "b--", linewidth = 1.0, name = "256")
# plotCrossSection(l, "g-", linewidth = 0.33, name = "256", pdfdoc = p)

# instanciate new solver
# m = laplace2DSolver(8) # 256 X 256

# m.mergeMask(Disk(0.10), "S") # CENTER 0.2 diameter
# m.mergeMask(Disk(0.40), "C") # SHELL 0.8 diameter 

# m.jacobiSteps(15000) # 256 X 256
# plotCrossSection(m, "g-", pdfdoc = p)

p.close()

# to do: add the analytical solution to compare convergences.
