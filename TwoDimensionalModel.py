from capycity.capycity import Document
from capycity.capycity import SolverTwoDimensions
from capycity.capycity import DiskAperture
from capycity.capycity import PlateSolid
from capycity.capycity import vplot
from capycity.capycity import headerText
from capycity.capycity import footerText
from capycity.capycity import selectfigure

from matplotlib.artist import Artist
from matplotlib.pyplot import semilogx

from numpy import arange
from numpy import linspace
from numpy import array
from numpy import append

SKIP_NUMERICS, SKIP_RESULTS = False, False

# width, height, gap, excentricity, shield radius
w, h, g, e, r = 8.400E-3, 1.500E-3, 1.000E-3, 0.000E-3, 50.000E-3

if not SKIP_NUMERICS:

    D = Document()
    D.opendocument("./local/maps.pdf")

    # lists instanciation
    PSS, CC, XX = [], [], []

    # build phase space section
    for g in linspace(0.300E-3, 1.200E-3, 10):
        PSS.append(
            ([r],                   # shield (radius)
            [+e, (+h+g)/2.0, w, h], # top    (x,y,w,h)
            [-e, (-h-g)/2.0, w, h], # bottom (x,y,w,h)
            [100]*11)               # series (steps, ...)
        )

    # run through parameters:
    # SHIELD, ELECTRODE1, ELECTRODE2, SERIES
    for SH, C1, C2, SERIES in PSS:

        n0 = 3

        S = SolverTwoDimensions(n0)
        n = 2**(n0+len(SERIES)-1)
        S.ll = 25E-6 * n

        # print(f"size = {S.ll*1E3}mm x {S.ll*1E3}mm")
        # print(f"resolution = {n} points = {S.ll/n*1E6}um x {S.ll/n*1E6}um")

        S.addMap("SH", DiskAperture(*SH)) # shield
        S.addMap("C1",   PlateSolid(*C1)) # plate 1
        S.addMap("C2",   PlateSolid(*C2)) # plate 2

        K, R = S.jacobiStepSeries(SERIES, "C1", "C2")

        # select plot figure
        fg, ax = selectfigure("Capacitance convergence")
        ax.semilogx(K, array(R)*1E12, ".-")
        ax.grid(True)
        
        # export and reset figure
        D.exportfigure("Capacitance convergence")
        ax.clear()

        # record last data point
        CC.append(R[-1]*1E12) # [pF]
        
        # record phase space section parameter
        _w, _h, _g, _e, _r = C1[2], C1[3], 2*C1[1]-h, C1[0], SH[0]

        x = _g

        XX.append(x)

        # select map figure
        fg, ax, bx = vplot(S, "Potential map", "C1")
        ax.set_title(f"V map")

        nx, ny = S.nn

        # HEADER
        txt  = f"width      = {_w*1E3:.3f}mm\n"
        txt += f"height     = {_h*1E3:.3f}mm\n"
        txt += f"gap        = {_g*1E3:.3f}mm\n"
        txt += f"excentr.   = {_e*1E3:.3f}mm\n"
        txt += f"shield rad.= {_r*1E3:.3f}mm\n"
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

        # export and reset figure
        D.exportfigure("Potential map")
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
    # 2048

    D1 = [
        0.000300, 215.541566,
        0.000413, 173.314126,
        0.000525, 145.490106,
        0.000638, 125.691685,
        0.000750, 111.018014,
        0.000862, 90.427512,
        0.000975, 82.840607,
        0.001088, 76.621922,
        0.001200, 71.191746,
        0.001312, 66.676491,
        0.001425, 62.681802,
        0.001538, 59.270301,
        0.001650, 56.119605,
        0.001762, 50.852179,
        0.001875, 48.737732,
        0.001988, 46.678013,
        0.002100, 44.902400,
        0.002213, 43.147926,
        0.002325, 41.683160,
        0.002437, 40.191834,
        0.002550, 39.007014,
        0.002662, 36.493071,
        0.002775, 35.420513,
        0.002888, 34.451420,
        0.003000, 33.453212,
        ]

    #################################################################
    # 4096

    D2 = [
        0.000300, 239.405924,
        0.000413, 187.974302,
        0.000525, 155.526830,
        0.000638, 124.330119,
        0.000750, 109.967487,
        0.000862, 98.709603,
        0.000975, 89.915650,
        0.001088, 79.089782,
        0.001200, 73.432762,
        0.001312, 68.709508,
        0.001425, 64.484970,
        0.001538, 59.150653,
        0.001650, 55.958258,
        0.001762, 53.371929,
        0.001875, 50.924880,
        0.001988, 47.674320,
        0.002100, 45.728990,
        0.002213, 44.017664,
        0.002325, 42.586138,
        0.002437, 40.239214,
        0.002550, 38.952309,
        0.002662, 37.690146,
        0.002775, 36.693690,
        0.002888, 34.967758,
        0.003000, 33.970595,
        ]

    #################################################################

    PLOTS = [
         (D1, "2048x2048", 1.0), 
         (D2, "4096x4096", 1.0), 
        ]

    #################################################################

    fg, ax = selectfigure("PLOT")
    units_x = 1E3  # [mm]
    # units_y = 1E12 # [pF] 

    m = 0.0 # (capacitance is always positive)
    if PLOTS:
        for DT, LB, RT in PLOTS:
            X, Y = array(DT[0::2]), array(DT[1::2])/RT
            ax.plot(X*units_x, Y,
                ".-",
                label = LB)
            m = max(Y)

    if not SKIP_NUMERICS:
        ax.plot(array(XX)*units_x, CC,
            "k.-",
            label = "last")
        m = max(m, max(CC))

    C = 8.854E-12*w/array(X)
    ax.plot(X*units_x, C*1E12, 
        '-.r', linewidth = 1.0,
        label = "Theoretical (no stray field)")
    m = max(m, max(C))

    ax.set_ylim(0, m*1.1)
    ax.grid(True)    
    
    ax.legend(loc = 'center')
    
    ax.set_xlabel("gap [mm]")
    ax.set_ylabel("capacitance [pF/m]")

    # HEADER
    txt  = f"width      = {w*1E3:.3f}mm\n"
    txt += f"height     = {h*1E3:.3f}mm\n"
    # txt += f"gap        = {g*1E3:.3f}mm\n"
    txt += f"excentr.   = {e*1E3:.3f}mm\n"
    txt += f"shield rad.= {r*1E3:.3f}mm\n"
    # txt += f"\n"
    # txt += f"resolution = {nx} X {ny}\n"
    txh = headerText(txt, fg)

    D.exportfigure("PLOT")
    D.closedocument()
