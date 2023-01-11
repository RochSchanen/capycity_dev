from capycity.capycity import Document
from capycity.capycity import SolverTwoDimensions
from capycity.capycity import DiskAperture
from capycity.capycity import PlateSolid
from capycity.capycity import vplot
from capycity.capycity import headerText
from capycity.capycity import footerText
from capycity.capycity import selectfigure

from matplotlib.artist import Artist

from numpy import arange
from numpy import array
from numpy import append

SKIP_NUMERICS, SKIP_RESULTS = True, False

if not SKIP_NUMERICS:

    D = Document()
    D.opendocument("./local/maps.pdf")

    w, h, g = 0.15, 0.03, 0.03

    NN = []
    CC = []
    SR = []

    sh = 0.0
    # for sh in arange(0.0, w/2, w/4/16):
    # for r in range(1, 9):    
    # for r in [0.20, 0.25, 0.30, 0.35, 0.40, 0.45]:    
    for r in arange(0.10, 0.95, 0.025):    
        NN.append(
            ([r],                    # shield (radius)
            [+sh, (+h+g)/2.0, w, h], # top    (x,y,w,h)
            [-sh, (-h-g)/2.0, w, h], # bottom (x,y,w,h)
            [100]*10)                # series (steps, ...)
        )

    for SH, C1, C2, SERIES in NN:

        S = SolverTwoDimensions(3)
        S.ll = 2.0

        S.addMap("SH", DiskAperture(*SH)) # shield
        S.addMap("C1",   PlateSolid(*C1)) # plate 1
        S.addMap("C2",   PlateSolid(*C2)) # plate 2

        K, R = S.jacobiStepSeries(SERIES, "C1", "C2")
        CC.append(R[-1]*1E12) # [pF]
        SR.append(SH[0])

        fg, ax, bx = vplot(S, "FIGURE", "C1")

        nx, ny = S.nn

        txt  = f"width = {w}mm\n"
        txt += f"height = {h}mm\n"
        txt += f"gap = {g}mm\n"
        txh = headerText(txt, fg)

        txt  = f"total steps = {K[-1]:6}\n"
        txt += f"capacity  = {R[-1]*1E12:.3f}pF\n"
        txf = footerText(txt, fg)

        print(f"capacity  = {R[-1]*1E12:.3f}pF\n")

        ax.set_title(f"{nx} X {ny}\n")

        D.exportfigure("FIGURE")
        ax.clear()

        Artist.remove(txh)
        Artist.remove(txf)

    fg, ax = selectfigure("PLOT")
    ax.plot(SR, CC, "r.-")
    ax.set_ylim(0, max(CC)*1.1)
    D.exportfigure("PLOT")

    D.closedocument()

if not SKIP_RESULTS:

    D = Document()
    D.opendocument("./local/results.pdf")

    ###########################################################################
    # 4096X4096 [500]*10

    D1 = [
        0.100, 44.699,
        0.125, 48.025,
        0.150, 49.811,
        0.175, 50.968,
        0.200, 51.784,
        0.225, 52.412,
        0.250, 52.910,
        0.275, 53.321,
        0.300, 53.671,
        0.325, 53.969,
        0.350, 54.232,
        0.375, 54.463,
        0.400, 54.669,
        0.425, 54.853,
        ]

    # computation time was 19723.4s, about 5h30, about 24min/point

    ###########################################################################
    # 1024X1024 [100]*8

    D2 = [
        0.100, 43.658,
        0.125, 46.998,
        0.150, 48.776,
        0.175, 49.875,
        0.200, 50.698,
        0.225, 51.301,
        0.250, 51.813,
        0.275, 52.236,
        0.300, 52.575,
        0.325, 52.863,
        0.350, 53.129,
        0.375, 53.374,
        0.400, 53.584,
        0.425, 53.784,
        ]

    # computation time was 214.9s, about 3.5min

    ###########################################################################
    # 4096X4096 [100]*10

    D3 = [
        0.100, 44.809,
        0.125, 48.155,
        0.150, 49.949,
        0.175, 51.105,
        0.200, 51.924,
        0.225, 52.554,
        0.250, 53.061,
        0.275, 53.483,
        0.300, 53.835,
        0.325, 54.120,
        0.350, 54.396,
        0.375, 54.639,
        0.400, 54.847,
        0.425, 55.055,
        ]

    # computation time was 3190.8s, about 53mins

    ###########################################################################
    # 8192X8192 [100]*11

    D4 = [
        0.100, 44.801,
        0.125, 48.148,
        0.150, 49.945,
        0.175, 51.100,
        0.200, 51.922,
        0.225, 52.548,
        0.250, 53.057,
        0.275, 53.481,
        0.300, 53.833,
        0.325, 54.117,
        0.350, 54.393,
        0.375, 54.638,
        0.400, 54.848,
        0.425, 55.056,
        0.450, 55.195,
        0.475, 55.384,
        0.500, 55.523,
        0.525, 55.642,
        0.550, 55.804,
        0.575, 55.905,
        ]

    # fcomputation failed at lack of memory: 
    # numpy.core._exceptions.MemoryError: Unable to allocate 512. MiB for
    # an array with shape (8192, 8190) and data type float64 in 'gradient':
    # out[tuple(slice1)] = (f[tuple(slice4)] - f[tuple(slice2)]) / (2. * ax_dx)

    # ###########################################################################
    # # 4096X4096 [100]*10  double length -> 2.0m

    D5 = [
        0.100, 43.468,
        0.125, 46.791,
        0.150, 48.574,
        0.175, 49.714,
        0.200, 50.549,
        0.225, 51.169,
        0.250, 51.674,
        0.275, 52.094,
        0.300, 52.441,
        0.325, 52.733,
        0.350, 53.004,
        0.375, 53.247,
        0.400, 53.455,
        0.425, 53.660,
        0.450, 53.805,
        0.475, 53.994,
        0.500, 54.134,
        0.525, 54.254,
        0.550, 54.415,
        0.575, 54.524,
        0.600, 54.661,
        0.625, 54.752,
        0.650, 54.835,
        0.675, 54.939,
        0.700, 55.040,
        0.725, 55.095,
        0.750, 55.170,
        0.775, 55.244,
        0.800, 55.268,
        0.825, 55.332,
        0.850, 55.293,
        0.875, 55.371,
        0.900, 55.448,
        0.925, 55.434,
        ]

# computation time was 7986.2s, about 2h10

#     ###########################################################################

    fg, ax = selectfigure("PLOT")

    X, Y = array(D1[0::2]), array(D1[1::2])
    ax.plot(X, Y, "b-", label = "4096X4096 [500]*10, 5h30min")
    m = max(Y)

    X, Y = array(D2[0::2]), array(D2[1::2])
    ax.plot(X, Y, "r-", label = "1024X1024 [100]*8 , 2.5min")
    m = max(m, max(Y))

    X, Y = array(D3[0::2]), array(D3[1::2])
    ax.plot(X, Y, "g-", label = "4096X4096 [100]*10 , 53min")
    m = max(m, max(Y))

    X, Y = array(D4[0::2]), array(D4[1::2])
    ax.plot(X, Y, "c-", label = "8192X8192 [100]*11 , ...")
    m = max(m, max(Y))

    X, Y = array(D5[0::2]), array(D5[1::2])
    ax.plot(X, Y, "y.-", label = "4096X4096 [100]*10 , double size")
    m = max(m, max(Y))

    if not SKIP_NUMERICS:
        ax.plot(SR, CC, "k.-", label = "last")
        m = max(m, max(CC))

    ax.set_ylim(0, m*1.1)
    ax.set_xlim(0, 1.0)
    
    ax.grid(True)
    
    ax.legend(loc = 'center')
    
    ax.set_xlabel("shield radius [m]")
    ax.set_ylabel("capacitance [pF]")

    D.exportfigure("PLOT")
    D.closedocument()
