from capycity.capycity import selectfigure
from capycity.capycity import Document

# from package: "https://numpy.org/"
from numpy import arange

# from package "https://scipy.org/"
from scipy.constants import epsilon_0 as EPS0

Y = [41.656, 41.591, 41.443, 41.211, 40.892,
     40.525, 40.058, 39.516, 38.935, 38.301,
     37.662, 36.958, 36.192, 35.419, 34.614,
     33.827, 32.993, 32.103, 31.222, 30.321,
     29.452, 28.541, 27.582, 26.642, 25.687,
     24.772, 23.818, 22.834, 21.871, 20.908,
     20.000, 19.068,]

fg, ax = selectfigure("FIG")

w, h, g = 0.4, 0.2, 0.1
s = arange(0.0, w/2, w/4/16)

ax.plot(s/w*100, Y, ".-")

C = (w-2*s)/g*EPS0
ax.plot(s/w*100, C*1E12, "r-")


D = Document("local/plot.pdf")
D.exportfigure("FIG")
D.closedocument()
