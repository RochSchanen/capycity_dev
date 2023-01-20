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

SKIP_NUMERICS, SKIP_RESULTS = False, False

if not SKIP_NUMERICS:

    D = Document()
    D.opendocument("./local/maps.pdf")

    w, h, g, sh = 0.15, 0.03, 0.03, 0.0

    sr = 0.95 # shield radius

    NN = []
    CC = []
    XX = []

    # for sh in arange(0.0, w/2, w/4/16):
    # for r in range(1, 9):    
    # for sr in arange(0.10, 0.95, 0.025):    
    # for sh in list(linspace(0.0, +w/2*1.25, 25)):    
    # for g in list(linspace(0.03, 0.30, 10)):    
    for h in list(linspace(0.01, w, 30)):    
        NN.append(
            ([sr],                   # shield (radius)
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
        x = C1[3]
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
    #  

    D1 = [
        ]

    PLOTS = [
        # (D1, "2048X2048, width = 0.150", 1.0), 
        ]

    #################################################################

    fg, ax = selectfigure("PLOT")

    m = 0.0 # (capacitance is always positive)
    if PLOTS:
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
    
    ax.set_xlabel("height [m]")
    ax.set_ylabel("capacitance [pF]")

    D.exportfigure("PLOT")
    D.closedocument()
