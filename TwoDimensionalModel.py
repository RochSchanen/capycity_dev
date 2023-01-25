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

SKIP_NUMERICS, SKIP_RESULTS = True, False

# units
unit_length         = (1E+03, "mm")
unit_capacitance    = (1E+12, "pF")

w = 8.400   # width
h = 1.500   # height
g = 1.000   # gap
e = 0.000   # excentricity
r = 50.000  # shield radius


ll = 100.0  # grid size
dd = 0.200  # resolution

N0 = 

if not SKIP_NUMERICS:

    D = Document()
    D.opendocument("./local/maps.pdf")

    # lists instanciation
    PSS, CC, XX = [], [], []

    # build phase space section
    for g in linspace(0.300, 1.500, 5):
        PSS.append(
            ([r],                   # shield (radius)
            [+e, (+h+g)/2.0, w, h], # top    (x,y,w,h)
            [-e, (-h-g)/2.0, w, h], # bottom (x,y,w,h)
            [100]*7)               # series (steps, ...)
        )

    # run through parameters:
    # SHIELD, ELECTRODE1, ELECTRODE2, SERIES
    for SH, C1, C2, SERIES in PSS:

        P0 = 3

        S = SolverTwoDimensions(P0)
        n = 2**(P0+len(SERIES)-1)
        S.ll = dxx * n

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
        
        # record parameters
        _w, _h, _g, _e, _r = C1[2], C1[3], 2*C1[1]-h, C1[0], SH[0]

        # set variable
        x = _g
        XX.append(x)

        # select map figure
        fg, ax, bx = vplot(S, "VMAP", "C1")
        ax.set_title(f"VMAP")

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
        D.exportfigure("VMAP")

        ax.clear()
        Artist.remove(txh)
        Artist.remove(txf)

    fg, ax = selectfigure("CPLOT")
    ax.plot(XX, CC, "r.-")
    ax.set_ylim(0, max(CC)*1.1)

    D.exportfigure("CPLOT")

    D.closedocument()

if not SKIP_RESULTS:

    D = Document()
    D.opendocument("./local/results.pdf")

    #################################################################
    # 2048, 50um

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
    # 4096, 25um

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
    # 8192, 50um

    D3 = [
        0.000300, 239.431792,
        0.000400, 187.991955,
        0.000500, 155.477476,
        0.000600, 132.961551,
        0.000700, 116.508285,
        0.000800, 103.915673,
        0.000900, 93.961525,
        0.001000, 85.818308,
        0.001100, 79.117618,
        0.001200, 73.479731,
        ]

    #################################################################
    # 1024, 100um

    D4 = [
        0.000300, 224.201314,
        0.000393, 149.592865,
        0.000486, 149.583064,
        0.000579, 113.216398,
        0.000672, 113.162972,
        0.000766, 91.654012,
        0.000859, 91.662060,
        0.000952, 77.436828,
        0.001045, 77.342553,
        0.001138, 67.211142,
        0.001231, 67.205979,
        0.001324, 59.623251,
        0.001417, 59.570289,
        0.001510, 53.658518,
        0.001603, 53.535728,
        0.001697, 53.535728,
        0.001790, 48.806020,
        0.001883, 48.773233,
        0.001976, 44.901341,
        0.002069, 44.896196,
        0.002162, 41.685919,
        0.002255, 41.590883,
        0.002348, 38.843673,
        0.002441, 38.854736,
        0.002534, 36.523652,
        0.002628, 36.407352,
        0.002721, 34.348294,
        0.002814, 34.345270,
        0.002907, 32.569139,
        0.003000, 32.541634,
        ]

    #################################################################
    # 8192, 12.5um

    D5 = [
        0.000300, 253.290925,
        0.000341, 221.027917,
        0.000383, 208.113070,
        0.000424, 186.411810,
        0.000466, 168.884841,
        0.000507, 161.115733,
        0.000548, 148.090665,
        0.000590, 137.157100,
        0.000631, 132.399831,
        0.000672, 123.679426,
        0.000714, 116.178944,
        0.000755, 112.747578,
        0.000797, 106.487956,
        0.000838, 100.927712,
        0.000879, 98.653388,
        0.000921, 93.832082,
        0.000962, 91.769374,
        0.001003, 87.556536,
        0.001045, 83.897449,
        0.001086, 82.395926,
        0.001128, 79.116708,
        0.001169, 76.214856,
        0.001210, 74.760346,
        0.001252, 72.261885,
        0.001293, 69.766471,
        0.001334, 68.703664,
        0.001376, 66.678018,
        0.001417, 64.504128,
        0.001459, 63.632070,
        0.001500, 61.710195,
        ]

    #################################################################

    PLOTS = [
         (D5, r"8192x8192, 12.5$\mu$m", 1.0),
         (D2, r"4096x4096, 25$\mu$m", 1.0),
         # (D1, r"2048x2048, 50$\mu$m", 1.0), 
         # (D3, r"8192x8192, 25$\mu$m", 1.0),
         # (D4, r"1024x1024, 100$\mu$m", 1.0),
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

    fg, ax = selectfigure("NORM")
    units_x = 1E3  # [mm]
    # units_y = 1E12 # [pF] 

    m = 0.0 # (capacitance is always positive)
    if PLOTS:
        for DT, LB, RT in PLOTS:
            X, Y = array(DT[0::2]), array(DT[1::2])/RT
            ax.plot(X*units_x, Y*X,
                ".-",
                label = LB)
            m = max(Y*X)

    C = 8.854E-12*w/array(X)
    ax.plot(X*units_x, C*1E12*X, 
        '-.r', linewidth = 1.0,
        label = "Theoretical (no stray field)")
    m = max(m, max(C*1E12*X))

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

    D.exportfigure("NORM")

    D.closedocument()
