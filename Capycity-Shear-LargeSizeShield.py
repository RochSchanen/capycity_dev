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
from numpy import linspace
from numpy import array
from numpy import append

SKIP_NUMERICS, SKIP_RESULTS = True, False

if not SKIP_NUMERICS:

    D = Document()
    D.opendocument("./local/maps.pdf")

    w, h, g = 0.15, 0.03, 0.03

    shieldRadius = 0.95

    NN = []
    CC = []
    XX = []

    sh = 0.0
    # for sh in arange(0.0, w/2, w/4/16):
    # for r in range(1, 9):    
    # for r in arange(0.10, 0.95, 0.025):    
    for sh in list(linspace(0.0, +w/2*1.25, 25)):    
        NN.append(
            ([0.95],                 # shield (radius)
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
        x = C1[0]
        XX.append(x)

        fg, ax, bx = vplot(S, "FIGURE", "C1")

        nx, ny = S.nn

        txt  = f"width      = {w:.3f}mm\n"
        txt += f"height     = {h:.3f}mm\n"
        txt += f"gap        = {g:.3f}mm\n"
        txt += f"shield     = {shieldRadius:.3f}mm\n"
        txt += f"\n"
        txt += f"resolution = {nx} X {ny}\n"
        txh = headerText(txt, fg)

        txt += f"SERIES: {SERIES}\n"        
        txt  = f"total steps = {K[-1]:6}\n"
        txt += f"capacity  = {R[-1]*1E12:.3f}pF\n"
        txf = footerText(txt, fg)

        print(f"{x:.6f}, {R[-1]*1E12:.6f},")

        ax.set_title(f"\n shear = {x:.3f}mm\n")

        D.exportfigure("FIGURE")
        ax.clear()

        Artist.remove(txh)
        Artist.remove(txf)

    fg, ax = selectfigure("PLOT")
    ax.plot(XX, CC, "r.-")
    ax.set_ylim(0, max(CC)*1.1)
    D.exportfigure("PLOT")

    D.closedocument()

if not SKIP_RESULTS:

    D = Document()
    D.opendocument("./local/results.pdf")

    ###########################################################################
    # 512X512

    D1 = [
        -0.150000, 9.346740,
        -0.143878, 9.543127,
        -0.137755, 10.204468,
        -0.131633, 11.025836,
        -0.125510, 11.369121,
        -0.119388, 12.374366,
        -0.113265, 12.785803,
        -0.107143, 14.068920,
        -0.101020, 14.670027,
        -0.094898, 16.280173,
        -0.088776, 17.222303,
        -0.082653, 19.428593,
        -0.076531, 22.287511,
        -0.070408, 23.721748,
        -0.064286, 27.026371,
        -0.058163, 28.419182,
        -0.052041, 31.854536,
        -0.045918, 33.257858,
        -0.039796, 36.248778,
        -0.033673, 39.374058,
        -0.027551, 40.388254,
        -0.021429, 43.251108,
        -0.015306, 43.876153,
        -0.009184, 45.814562,
        -0.003061, 45.570014,
         0.003061, 45.570014,
         0.009184, 45.814562,
         0.015306, 43.876153,
         0.021429, 43.251108,
         0.027551, 40.388254,
         0.033673, 39.374058,
         0.039796, 36.248778,
         0.045918, 33.257858,
         0.052041, 31.854536,
         0.058163, 28.419182,
         0.064286, 27.026371,
         0.070408, 23.721748,
         0.076531, 22.287511,
         0.082653, 19.428593,
         0.088776, 17.222303,
         0.094898, 16.280173,
         0.101020, 14.670027,
         0.107143, 14.068920,
         0.113265, 12.785803,
         0.119388, 12.374366,
         0.125510, 11.369121,
         0.131633, 11.025836,
         0.137755, 10.204468,
         0.143878, 9.543127,
         0.150000, 9.346740,
    ]

    ###########################################################################
    # 1024X1024

    D2 = [
        0.000000, 49.635726,
        0.003906, 49.439299,
        0.007812, 48.977297,
        0.011719, 48.205808,
        0.015625, 47.504734,
        0.019531, 46.333720,
        0.023438, 45.031983,
        0.027344, 43.596278,
        0.031250, 42.146224,
        0.035156, 40.561581,
        0.039062, 38.954155,
        0.042969, 37.259403,
        0.046875, 35.628752,
        0.050781, 33.875141,
        0.054688, 32.124028,
        0.058594, 30.315632,
        0.062500, 28.603156,
        0.066406, 26.797351,
        0.070312, 25.081515,
        0.074219, 23.381471,
        0.078125, 21.804877,
        0.082031, 20.342445,
        0.085938, 19.063096,
        0.089844, 17.887263,
        0.093750, 16.879851,
    ]

    ###########################################################################
    # 2048X2048

    D3 = [
        0.000000, 55.527413,
        0.003906, 55.321288,
        0.007812, 54.745673,
        0.011719, 53.849607,
        0.015625, 52.925323,
        0.019531, 51.588836,
        0.023438, 50.059562,
        0.027344, 48.430625,
        0.031250, 46.735181,
        0.035156, 44.956064,
        0.039062, 43.100944,
        0.042969, 41.204526,
        0.046875, 39.314958,
        0.050781, 37.367003,
        0.054688, 35.362007,
        0.058594, 33.341975,
        0.062500, 31.362174,
        0.066406, 29.334221,
        0.070312, 27.347733,
        0.074219, 25.411025,
        0.078125, 23.568429,
        0.082031, 21.903531,
        0.085938, 20.420819,
        0.089844, 19.093616,
        0.093750, 17.938006,
    ]

    ###########################################################################
    # 4096X4096

    D4 = [
        0.000000, 55.393884,
        0.003906, 55.196020,
        0.007812, 54.626532,
        0.011719, 53.747880,
        0.015625, 52.824213,
        0.019531, 51.509775,
        0.023438, 49.993409,
        0.027344, 48.388664,
        0.031250, 46.702551,
        0.035156, 44.948786,
        0.039062, 43.106684,
        0.042969, 41.235873,
        0.046875, 39.352120,
        0.050781, 37.429144,
        0.054688, 35.436850,
        0.058594, 33.441669,
        0.062500, 31.468594,
        0.066406, 29.463363,
        0.070312, 27.485356,
        0.074219, 25.563519,
        0.078125, 23.710943,
        0.082031, 22.049805,
        0.085938, 20.559798,
        0.089844, 19.229046,
        0.093750, 18.060443,
    ]

    fg, ax = selectfigure("PLOT")

    X, Y = array(D1[0::2]), array(D1[1::2])
    ax.plot(X, Y, "b.-", label = "512X512")
    m = max(Y)

    X, Y = array(D2[0::2]), array(D2[1::2])
    ax.plot(X, Y, "c.-", label = "1024X1024")
    m = max(m, max(Y))

    X, Y = array(D3[0::2]), array(D3[1::2])
    ax.plot(X, Y, "r.-", label = "2048X2048")
    m = max(m, max(Y))

    X, Y = array(D4[0::2]), array(D4[1::2])
    ax.plot(X, Y, "k-", label = "4096X4096")
    m = max(m, max(Y))

    if not SKIP_NUMERICS:
        ax.plot(XX, CC, "k.-", label = "last")
        m = max(m, max(CC))

    ax.set_ylim(0, m*1.1)
    # ax.set_xlim(0, 1.0)
    
    ax.grid(True)
    
    ax.legend(loc = 'center')
    
    ax.set_xlabel("shield radius [m]")
    ax.set_ylabel("capacitance [pF]")

    D.exportfigure("PLOT")
    D.closedocument()
