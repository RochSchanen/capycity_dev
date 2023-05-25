text = """
	(Col("Angle")+6.5)		(Col("Capacitance")/1.87)
Angle	Angle_new	Capacitance	Capacitance_new
deg	deg	pF	pF
	Ref to Roch		Dielectric constant=1.87
0	6.5	43.1	23.048128342246
-1	5.5	41.97	22.4438502673797
-2	4.5	39.81	21.2887700534759
-3	3.5	37.86	20.2459893048128
-4	2.5	34.89	18.6577540106952
-5	1.5	33.11	17.7058823529412
-6	0.5	31.36	16.7700534759358
-7	-0.5	28.73	15.3636363636364
-8	-1.5	26.196	14.0085561497326
-9	-2.5	23.679	12.6625668449198
-10	-3.5	21.66	11.5828877005348
-11	-4.5	18.937	10.1267379679144
-12	-5.5	15.906	8.50588235294118
-13	-6.5	14.65	7.83422459893048
-14	-7.5	13.876	7.42032085561497
-15	-8.5	13.59	7.26737967914438
-16	-9.5	13.77	7.36363636363636
-17	-10.5	14.01	7.49197860962567
-18	-11.5	15.215	8.13636363636364
-19	-12.5	17.187	9.19090909090909
-20	-13.5	19.18	10.2566844919786
-21	-14.5	21.72	11.6149732620321
-22	-15.5	23.22	12.4171122994652
-23	-16.5	26.133	13.9748663101604
1	7.5	41.597	22.244385026738
2	8.5	39.75	21.2566844919786
3	9.5	37.9	20.2673796791444
4	10.5	34.8	18.6096256684492
5	11.5	33.2	17.7540106951872
6	12.5	30.77	16.4545454545455
7	13.5	27.97	14.9572192513369
8	14.5	25.92	13.8609625668449
9	15.5	23.08	12.3422459893048
10	16.5	21.15	11.3101604278075
11	17.5	18.18	9.72192513368984
12	18.5	15.74	8.41711229946524
13	19.5	14.74	7.88235294117647
14	20.5	14.02	7.49732620320856
15	21.5	13.76	7.35828877005348
16	22.5	13.9	7.4331550802139
17	23.5	14.61	7.81283422459893
18	24.5	15.72	8.40641711229947
19	25.5	17.76	9.49732620320856
20	26.5	20.61	11.0213903743315
21	27.5	23.9	12.7807486631016
22	28.5	25.75	13.7700534759358
23	29.5	27.7	14.8128342245989
24	30.5	29.52	15.7860962566845
25	31.5	32.49	17.3743315508021
26	32.5	34.56	18.4812834224599
27	33.5	37.64	20.1283422459893
28	34.5	39.79	21.2780748663102
29	35.5	41.06	21.9572192513369
30	36.5	41.86	22.3850267379679
31	37.5	41.56	22.2245989304813
32	38.5	39.9	21.3368983957219
33	39.5	37.92	20.2780748663102
34	40.5	34.75	18.5828877005348
35	41.5	33.08	17.6898395721925
36	42.5	30.33	16.2192513368984
37	43.5	28.19	15.0748663101604
"""

# get list of line from text
A, C = [], []
for l in text.split(f"\n")[5:]:
	V = l.split("\t")
	if len(V) == 4:
		A.append(float(V[0]))
		C.append(float(V[2]))

# convert lists into numpy arrays
from numpy import array
A, C = array(A), array(C)

####################
### USEFUL TOOLS ###
####################

def selectfigure(name):
    from matplotlib.pyplot import figure, fignum_exists
    if not fignum_exists(name):
        # create figure
        fg = figure(name)
        # set A4 paper dimensions
        fg.set_size_inches(8.2677, 11.6929)
        # create square axis
        w, h = array([1, 1 / 1.4143])*0.7
        x, y = (1-w)/2, (1-h)/2
        ax = fg.add_axes([x, y, w, h])
    else:
        # select figure
        # (here the figure can be of any type)
        fg = figure(name)
        # get axes
        ax = fg.get_axes()[0]
    # done
    return fg, ax

class Document():

    def __init__(self, pathname = None):
        if pathname is not None:
            self._DOC = self.opendocument(pathname)
        return

    def opendocument(self, pathname):
        from matplotlib.backends.backend_pdf import PdfPages
        self._DOC = PdfPages(pathname)
        return self._DOC

    def exportfigure(self, name):
        args = selectfigure(name)
        self._DOC.savefig(args[0])
        return

    def closedocument(self):
        self._DOC.close()
        return

#########################
### IDEAL CAPACITANCE ###
#########################

# the following function compute the overlap angle
# from the angle of rotation (the origin is defined
# when the overlap is maximum)
def overlap(a):
    # a: angle in degrees (origin at maximum)
    
    from numpy import mod, fabs
    # mod(x,1.0) returns the fractional part of x
    # fabs(x) returns the absolute value of x

    # define boundaries, M: maximum overlap, T: period
    M, T = 12.0, 30.0
    
    # compute overlap from angle
    O = M-T*abs(mod(a/T+0.5, 1.0)-0.5)
    
    # nullify negative values
    O[O<0] = 0
    
    # done
    return O 

def capacitance(angle, gap):
    # angle: angle in degrees
    # gap: gap in metres

    from scipy.constants import epsilon_0 as E0
    from numpy import pi, square, fabs

    # overlap "eta" in radiants
    eta = overlap(angle)*pi/180.0

    # define star shape capacitance parameters
    R1, R2 = 40E-3/2, 15E-3/2
    N, L, l = 12, R1-R2, (R1+R2)/2.0

    # compute star shape capacitance
    C0 = N*E0*L*l/gap*eta

    # define ring capacitance parameters
    R3 = 13E-3/2

    # compute ring capacitance
    C1 = E0*(pi*square(R2)-pi*square(R3))/gap

    # done
    return (C0+C1)*1E12 # [pF]

########################
### PDF PLOT RESULTS ###
########################

# 40µm plastic sheet parameters:
# e: measured dielectric constant relative to vacuum (epsilon)
# g: measured thickness (capacitance gap)
e, g = 1.818, 210E-6

fg, ax = selectfigure("capa-vs-angle")

# plot C versus A (shifted by 6.5 degree to set the origin)
ax.plot(A, C/e, "k.-")

# same plot shifed by "nearly" one period (30.5 degerees)
# ax.plot(A-6.5 + 30.5, C/e, "k.-")

# analytical plot
ax.plot(A, capacitance(A, g), "r-")
# ax.plot(A + 30.5, capacitance(A, g), "r-")

# add labels
ax.set_xlabel("angle [degrees]")
ax.set_ylabel("Capacitance [pF]")

# add vertical and horizontal guide lines
ax.vlines(  [list(array([-12, -6, 0, +6, +12]))],
            -10, 110, "b", "dashed", linewidth = 0.50)
ax.vlines(  [list(array([-15, +15]))],
            -10, 110, "r", "dashed", linewidth = 0.50)

# modify ticks
# from numpy import linspace
# ax.set_xticks(linspace(-24, 48, 1+(48+24)//12))
# ax.set_xticks(linspace(-24, 48, 1+(48+24)//6), minor = True)
# ax.set_yticks(linspace(-20, 120, 1+(120+20)//10))
# ax.set_yticks(linspace(-20, 120, 1+(120+20)//5), minor = True)

# show grid on both axes
# ax.grid("True", axis = "y", which = "both")
ax.grid("True", axis = "both", which = "both")

# export figure to pdf document
doc = Document("../result.pdf") # file path + name
doc.exportfigure("capa-vs-angle") # figure name
doc.closedocument()
