
# from numpy import log, ceil, logspace, diff, array

################################################

# S = [100, 100, 100, 100, 100, 100]
# m = ceil(log(sum(S)) / log(10))
# P = list(logspace(1, m, 30, dtype = 'int'))
# P.append(sum(S))
# P.sort()

# print(m)
# print(P)

################################################

# print([10, 20]+[30, 40])
# print([10, 20]+[])
# print([1, 2, 3, 4, 5][:-1])
# print(array([[1, 2], [3, 4]])[:, 0])
# print([0, 1, 2, 3, 4][-1])
# print([0, 1, 2, 3, 4][0:-1])

#################################################

# from numpy import exp, log
# from numpy import array
# from numpy import ceil, floor

# def sep():
# 	print(f"{'-'*50}")
# 	return

# sep()

# X = 1.23456*exp((array(range(13))-6)*log(10.0))
# print(X)

# sep()

# X = 1.23456*exp((array(range(13))-6)*log(10.0))
# print(floor(log(X)/log(10)))
# print(log(X)/log(10)//3)
# print(floor(log(X)/log(1000)))

# sep()

# X = 1.23456*exp((array(range(13))-6)*log(10.0))
# D = -floor(log(X)/log(1000))
# print([x*1000**d for x, d in zip(X, D)])
# print([f"E{int(-3*d):+03d}" for d in D])

# sep()

# X = 1.23456*exp((array(range(13))-6)*log(10.0))
# print([f"{x*1000**d:.3f}E{int(-3*d):+03d}" for x, d in zip(X, -floor(log(X)/log(1000)))])

# from numpy import floor, log
# def unf(x):
#     d = -floor(log(x)/log(1000))
#     return f"{x*1000**d:3.2f}{'yzafpnµm kMGTPEZY'[int(8-d)]}"




######################################################################

# SELECT = [
#     "SERIES-CONVERGENCE",  # 0
#     "SERIES-COMPARISON",   # 1
#     "MULTI-GRID",          # 2
#     ][2]

# if __name__ == "__main__":

#     time_compute = time()

#     if SELECT == "CONVERGENCE":
        
#         from matplotlib.artist import Artist

#         D = Document()
#         D.opendocument("../local/plot.pdf")

#         w, h, g, e, r = 0.3, 0.05, 0.05, 0.0, 0.45

#         NN = []

#         # for n in [1, 5, 10, 50, 100, 200]:
#         for n in [3, 10, 30, 100, 300]:
#             NN.append(
#                 ([r],                       # shield (radius)
#                 [+e, (+h+g)/2.0, w, h],     # top    (x,y,w,h)
#                 [-e, (-h-g)/2.0, w, h],     # bottom (x,y,w,h)
#                 [n]*8)                      # series (steps, ...)
#             )

#         X, Y = [], []

#         for SH, C1, C2, N in NN:

#             S = SolverTwoDimensions()

#             S.addMap("SH", DiskAperture(*SH)) # SHIELD
#             S.addMap("C1",   PlateSolid(*C1)) # PLATE1
#             S.addMap("C2",   PlateSolid(*C2)) # PLATE2

#             K, C = S.jacobiStepSeries(N, "C1", "C2")

#             fg, ax = selectfigure("PLOT")
#             ax.semilogx(K, array(C)*1E12, ".-")
#             ax.set_xlabel("iterations")
#             ax.set_ylabel("Capacitance computed [pF/m]")

#             X.append(K[-1])
#             Y.append(C[-1])
#             print(K[-1], C[-1])

#         fg, ax = selectfigure("RESULTS")
#         ax.semilogx(X, array(Y)*1E12, ".-")

#         ax.set_xlabel("iterations")
#         ax.set_ylabel("Capacitance computed [pF/m]")
#         ax.set_ylim(0.0, 100.0)

#         D.exportfigure("PLOT") 
#         D.exportfigure("RESULTS") 

#         D.closedocument()

#     if SELECT == "COMPARISON":

#         DATA = array([
#             [
#             3, 3.0436270563981022e-12,
#             10, 5.1216678036887545e-12,
#             30, 5.388440551520647e-12,
#             100, 5.384585983847195e-12,
#             300, 5.3845859348425295e-12,
#             ],
#             [
#             7, 1.3545469497385262e-11,
#             21, 1.4124039507603287e-11,
#             61, 1.3673979676305713e-11,
#             200, 1.3576555499625977e-11,
#             600, 1.3576477010394754e-11,
#             ],
#             [
#             10, 2.2063907663026516e-11,
#             32, 2.3459003810195627e-11,
#             92, 2.3688678652095264e-11,
#             300, 2.386603938645882e-11,
#             900, 2.3916223989714314e-11,
#             ],
#             [
#             15, 3.4069942354760546e-11,
#             43, 3.634845120054825e-11,
#             123, 3.638886271302465e-11,
#             401, 3.6355243645236484e-11,
#             1200, 3.637257009179239e-11,
#             ],
#             [
#             19, 5.0055143884454275e-11,
#             54, 5.148642724515242e-11,
#             154, 5.182449505402855e-11,
#             500, 5.176298878944347e-11,
#             1500, 5.1766003171494677e-11,
#             ],
#             [
#             23, 5.893946748795823e-11,
#             65, 5.780413511747902e-11,
#             184, 5.78213629295024e-11,
#             600, 5.784758349382356e-11,
#             1800, 5.785180109663687e-11,
#             ],
#             [
#             27, 6.04384253286411e-11,
#             76, 5.851355714867663e-11,
#             215, 5.817602289962751e-11,
#             700, 5.7952018021336755e-11,
#             2100, 5.7888629596266815e-11,
#             ],
#             [
#             31, 6.332311209754091e-11,
#             87, 6.048215759261286e-11,
#             247, 5.992852197701082e-11,
#             801, 5.971289191392915e-11,

#             2400, 5.967733448169437e-11,
#             ],
#             [
#             35, 6.832368567458584e-11,
#             98, 6.308824560138485e-11,
#             275, 6.205336337416299e-11,
#             907, 6.169605168399679e-11,
#             2700, 6.158444576691678e-11,
#             ],
#             [
#             39, 7.113133410465773e-11,
#             100, 6.415771689528503e-11,
#             308, 6.27469380545123e-11,
#             1000, 6.227336547622629e-11,
#             3000, 6.209861777818127e-11,
#             ],
#             ])

#         from numpy import size, shape

#         D = Document()
#         D.opendocument("../local/plot.pdf")

#         fg, ax = selectfigure("RESULTS1")
#         N = 8
#         for j in range(shape(DATA)[0]):
#             X, Y = DATA[j, 0::2], DATA[j, 1::2]
#             ax.semilogx(X, Y*1E12, ".-", label = f"{N}x{N}")
#             N <<= 1

#         ax.set_xlabel("iterations")
#         ax.set_ylabel("Computed capacitance[pF/m]")
#         ax.legend()
#         D.exportfigure("RESULTS1") 

#         fg, ax = selectfigure("RESULTS2")
#         X, Y = DATA[:, -2], DATA[:, -1]
#         ax.semilogx(X, Y*1E12, ".-")

#         ax.set_xlabel("iterations")
#         ax.set_ylabel("Computed capacitance[pF/m]")
#         D.exportfigure("RESULTS2") 

#         D.closedocument()

#     if SELECT == "MULTI-GRID":

#         # electrodes position and geometry
#         Width   = 8.400
#         Height  = 1.500
#         Gap     = 0.300
#         Offset  = 0.000

#         # main grid 
#         Size, Length = 8, 102.4

#         # shield geometry
#         Radius = (Length*0.95)/2.0

#         # initial map series
#         MainSeries = [100]*5

#         # subnet map series
#         SubnetSeries = [
#             (80, 100),
#             (60, 100),
#             (40, 100),
#             (20, 100),
#             ]

#         # INIT SOLVER
#         S = SolverTwoDimensions(n = Size, l = Length)
#         S.addPart("C1", PlateSolid(0.0, +(Height+Gap)/2, Width, Height))
#         S.addPart("C2", PlateSolid(0.0, -(Height+Gap)/2, Width, Height))
#         S.addPart("SH", DiskAperture(Radius))
        
#         def Step():
#             S.M["C1"].jacobiStep()
#             S.M["C2"].jacobiStep()
#             # S.M["SH"].jacobiStep()

#         # main series
#         for N in MainSeries[:-1]:
#             for n in range(N): Step()
#             S.incrementResolution()
#         # last of the main series
#         for n in range(MainSeries[-1]): Step()

#         # subnet series
#         for L, N in SubnetSeries:
#             # ADD SUBNET
#             S.addSubnet(L)
#             for n in range(N): Step()

#         # DISPLAY

#         PP = ["P1"]

#         fg, ax, bx = mplot(S, "P1", "C1")
#         if S.M["C1"].S: splot(ax, S.M["C1"])

#         s = S.M["C1"].S
#         while s:
#             PN = f"P{len(PP)+1}"
#             PP.append(PN)
#             fg, ax, bx = mapPlot(s, PN)
#             if s.S: splot(ax, s)
#             s = s.S

#         # BUILD PDF
#         D = Document()
#         D.opendocument("../local/plot.pdf")
#         for p in PP:
#             D.exportfigure(p)
#         D.closedocument()

#     # done
#     time_done = time()
#     print(f"load in {time_compute - time_start}s")
#     print(f"compute in {time_done - time_compute}s")

