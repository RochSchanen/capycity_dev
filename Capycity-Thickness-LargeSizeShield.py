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

    # width, height, gap, excentricity, shield radius
    w, h, g, e, r = 8.400E-3, 1.500E-3, 0.300E-3, 0.000E-3, 100E-3

    # lists instanciation
    PSS, CC, XX = [], [], []

    # build phase space section
    for i in [0]:
        PSS.append(
            ([r],                   # shield (radius)
            [+e, (+h+g)/2.0, w, h], # top    (x,y,w,h)
            [-e, (-h-g)/2.0, w, h], # bottom (x,y,w,h)
            [100]*8)                # series (steps, ...)
        )

    # run through parameters:
    # SHIELD, ELECTRODE1, ELECTRODE2, SERIES
    for SH, C1, C2, SERIES in PSS:

        n0 = 3

        S = SolverTwoDimensions(n0)
        n = 2**(n0+len(SERIES)-1)
        print(n)

        S.ll = 100E-6*n # 100 um resolution
        print(S.ll)

        S.addMap("SH", DiskAperture(*SH)) # shield
        S.addMap("C1",   PlateSolid(*C1)) # plate 1
        S.addMap("C2",   PlateSolid(*C2)) # plate 2

        K, R = S.jacobiStepSeries(SERIES, "C1", "C2")
        CC.append(R[-1]*1E12) # [pF]
        h = C1[3]
        x = C1[3]
        XX.append(x)

        # create/select new figure
        fg, ax, bx = vplot(S, "FIGURE", "C1")
        ax.set_title(f"V map")

        nx, ny = S.nn

        # HEADER
        txt  = f"width      = {w:.3f}mm\n"
        txt += f"height     = {h:.3f}mm\n"
        txt += f"gap        = {g:.3f}mm\n"
        txt += f"shield rad.= {sr:.3f}mm\n"
        txt += f"\n"
        txt += f"resolution = {nx} X {ny}\n"
        txh = headerText(txt, fg)

        # FOOTER
        txt += f"SERIES: {SERIES}\n"        
        txt  = f"total steps = {K[-1]:6}\n"
        txt += f"capacity  = {R[-1]*1E12:.3f}pF\n"
        txf = footerText(txt, fg)

        # terminal output
        print(f"{x:.6f}, {R[-1]*1E12:.6f},")

        # add new page to document
        D.exportfigure("FIGURE")

        # reset plot and artistry
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
    # 1024

    D1 = [
        0.003000, 37.859640,
        0.015375, 41.130691,
        0.027750, 43.943739,
        0.040125, 44.925251,
        0.052500, 47.508757,
        0.064875, 48.097791,
        0.077250, 48.900367,
        0.089625, 49.293895,
        0.102000, 49.715724,
        0.114375, 50.628950,
        0.126750, 50.904051,
        0.139125, 51.220791,
        0.151500, 51.546386,
        0.163875, 51.781697,
        0.176250, 51.990669,
        0.188625, 52.162298,
        0.201000, 52.380134,
        0.213375, 52.531481,
        0.225750, 52.770311,
        0.238125, 52.954906,
        0.250500, 53.126443,
        0.262875, 53.245252,
        0.275250, 53.347620,
        0.287625, 53.476203,
        0.300000, 53.637300,
        ]

    #################################################################
    # 2048

    D2 = [
        0.003000, 41.906523,
        0.015375, 45.346863,
        0.027750, 48.107187,
        0.040125, 49.028652,
        0.052500, 51.731893,
        0.064875, 52.289262,
        0.077250, 52.857983,
        0.089625, 53.350583,
        0.102000, 53.765926,
        0.114375, 54.655786,
        0.126750, 54.992924,
        0.139125, 55.295213,
        0.151500, 55.511009,
        0.163875, 55.798751,
        0.176250, 56.008431,
        0.188625, 56.169167,
        0.201000, 56.425590,
        0.213375, 56.573636,
        0.225750, 56.737730,
        0.238125, 56.962532,
        0.250500, 57.131702,
        0.262875, 57.242782,
        0.275250, 57.375288,
        0.287625, 57.500818,
        0.300000, 57.606601,
    ]

    #################################################################
    # 4096

    D3 = [
        0.003000, 48.538379,
        0.015375, 51.185378,
        0.027750, 54.080220,
        0.040125, 54.948012,
        0.052500, 57.558844,
        0.064875, 58.187439,
        0.077250, 58.727975,
        0.089625, 59.173037,
        0.102000, 59.640085,
        0.114375, 60.504708,
        0.126750, 60.803406,
        0.139125, 61.140281,
        0.151500, 61.346838,
        0.163875, 61.609423,
        0.176250, 61.850550,
        0.188625, 62.027919,
        0.201000, 62.235305,
        0.213375, 62.405915,
        0.225750, 62.585581,
        0.238125, 62.773826,
        0.250500, 62.964280,
        0.262875, 63.087236,
        0.275250, 63.184047,
        0.287625, 63.326869,
        0.300000, 63.446725,
        ]

    #################################################################

    PLOTS = [
         (D1, "1024X1024", 1.0), 
         (D2, "2048X2048", 1.0), 
         (D3, "4096X4096", 1.0), 
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
