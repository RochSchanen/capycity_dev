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

    w, h, g, sh = 0.60, 0.03, 0.03, 0.0

    sr = 0.95 # shield radius

    NN = []
    CC = []
    XX = []

    # for sh in arange(0.0, w/2, w/4/16):
    # for r in range(1, 9):    
    # for sr in arange(0.10, 0.95, 0.025):    
    # for sh in list(linspace(0.0, +w/2*1.25, 25)):    
    for g in list(linspace(0.03, 0.30, 20)):    
        NN.append(
            ([sr],                   # shield (radius)
            [+sh, (+h+g)/2.0, w, h], # top    (x,y,w,h)
            [-sh, (-h-g)/2.0, w, h], # bottom (x,y,w,h)
            [100]*9)                 # series (steps, ...)
        )

    for SH, C1, C2, SERIES in NN:

        S = SolverTwoDimensions(3)
        S.ll = 2.0

        S.addMap("SH", DiskAperture(*SH)) # shield
        S.addMap("C1",   PlateSolid(*C1)) # plate 1
        S.addMap("C2",   PlateSolid(*C2)) # plate 2

        K, R = S.jacobiStepSeries(SERIES, "C1", "C2")
        CC.append(R[-1]*1E12) # [pF]
        x = C1[1]
        XX.append(x)

        fg, ax, bx = vplot(S, "FIGURE", "C1")

        nx, ny = S.nn

        txt  = f"width      = {w:.3f}mm\n"
        txt += f"height     = {h:.3f}mm\n"
        txt += f"gap        = {g:.3f}mm\n"
        txt += f"shield     = {sr:.3f}mm\n"
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

    #################################################################
    # 2048X2048, W = 0.150

    D1 = [
        0.030000, 55.527413,
        0.037105, 39.108748,
        0.044211, 31.605161,
        0.051316, 26.897743,
        0.058421, 23.377755,
        0.065526, 20.323287,
        0.072632, 18.429379,
        0.079737, 16.835642,
        0.086842, 15.241471,
        0.093947, 14.165333,
        0.101053, 13.301016,
        0.108158, 12.462907,
        0.115263, 11.551141,
        0.122368, 10.947341,
        0.129474, 10.382846,
        0.136579, 9.882326,
        0.143684, 9.273518,
        0.150789, 8.905575,
        0.157895, 8.473548,
        0.165000, 8.033780,
    ]

    #################################################################
    # 2048X2048, W = 0.300

    D2 = [
        0.030000, 99.407535,
        0.037105, 68.393053,
        0.044211, 54.299423,
        0.051316, 45.282089,
        0.058421, 38.931713,
        0.065526, 33.513150,
        0.072632, 30.110139,
        0.079737, 27.296282,
        0.086842, 24.573709,
        0.093947, 22.700808,
        0.101053, 21.169655,
        0.108158, 19.764682,
        0.115263, 18.255978,
        0.122368, 17.224038,
        0.129474, 16.267942,
        0.136579, 15.455744,
        0.143684, 14.449421,
        0.150789, 13.823681,
        0.157895, 13.136891,
        0.165000, 12.445509,
    ]   
    
    #################################################################
    # 2048X2048, W = 0.600

    D3 = [
        0.030000, 185.362063,
        0.037105, 125.408468,
        0.044211, 98.252341,
        0.051316, 81.064990,
        0.058421, 69.099389,
        0.065526, 59.059190,
        0.072632, 52.652582,
        0.079737, 47.434690,
        0.086842, 42.494395,
        0.093947, 39.073734,
        0.101053, 36.231308,
        0.108158, 33.719799,
        0.115263, 31.087680,
        0.122368, 29.216517,
        0.129474, 27.521349,
        0.136579, 26.098961,
        0.143684, 24.337495,
        0.150789, 23.205562,
        0.157895, 22.051027,
        0.165000, 20.877991,
        ]

    PLOTS = [
        (D1, "2048X2048, width = 0.150", 1.0), 
        (D2, "2048X2048, width = 0.300", 1.0),
        (D3, "2048X2048, width = 0.600", 1.0),
        ]


    fg, ax = selectfigure("PLOT")

    m = max(array(D1[1::2]))
    for DT, LB, RT in PLOTS:
        X, Y = array(DT[0::2]), array(DT[1::2])/RT
        ax.plot(X, Y, ".-", label = LB)
        m = max(Y)

    if not SKIP_NUMERICS:
        ax.plot(XX, CC, "k.-", label = "last")
        m = max(m, max(CC))

    ax.set_ylim(0, m*1.1)
    # ax.set_xlim(0, 1.0)
    
    ax.grid(True)
    
    ax.legend(loc = 'center')
    
    ax.set_xlabel("gap [m]")
    ax.set_ylabel("capacitance [pF]")

    D.exportfigure("PLOT")
    D.closedocument()
