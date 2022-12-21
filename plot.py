from numpy import array

D1 = array([
	76.92875155399936,
	72.79131261907658,
	67.65944001449182,
	62.09613653469326,
	56.05349182231327,
	50.04582315179711,
	43.61266302375045,
	37.48049098214236,
	30.93866978670963,
	24.69720674405142,
	18.50252884903883,
])

D2 = array([
	79.5182424456959,
	75.1451644976990,
	69.7963495393127,
	64.0731095943159,
	57.8109967635211,
	51.6321205973807,
	44.9632507226511,
	38.6506382171436,
	31.8713274916639,
	25.4464547011499,
	19.0142997189037,
	13.7869992292058,
])

D3 = array([
80.0949495762547,
75.6089870718837,
70.2337676860025,
64.4566225631840,
58.2201243959246,
51.9825056593432,
45.1251091385027,
38.7367325838053,
31.9600148426267,
25.4619917644131,
19.0110355881078,
13.7595076767761,
])

from capycity.graphics import selectfigure, exportfigure
from capycity.graphics import opendocument, closedocument

fg, ax = selectfigure("FIGURE")

ax.plot(D1,"b.-")
ax.plot(D2,"r.-")
ax.plot(D3,"g.-")

# local parameters
r = 0.9 / 2 # shield radius
w = 0.4     # width
h = 0.1     # height
g = 0.05    # gap 

from scipy.constants import epsilon_0 as e0
from numpy import arange

C = e0 * w / g * arange(1.0, -0.1, -0.1)
ax.plot(C*1E12,"k--")

opendocument("./local/plot2.pdf")
exportfigure("FIGURE") 
closedocument()
