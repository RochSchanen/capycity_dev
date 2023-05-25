text = """
	(Col("Angle")+6.5)		(Col("Capacitance")/1.87)
Angle	Angle	Capacitance	Capacitance_new
deg	deg	pF	pF
	ref to Roch		Dielectric constant=1.87
0	6.5	76.4	40.855614973262
-1	5.5	73.52	39.3155080213904
-2	4.5	69.44	37.1336898395722
-3	3.5	64.45	34.4652406417112
-4	2.5	60.29	32.2406417112299
-5	1.5	54.9	29.3582887700535
-6	0.5	50.489	26.9994652406417
-7	-0.5	46.113	24.6593582887701
-8	-1.5	41.13	21.9946524064171
-9	-2.5	36.579	19.5609625668449
-10	-3.5	32.01	17.1176470588235
-11	-4.5	27.09	14.4866310160428
-12	-5.5	22.419	11.9887700534759
-13	-6.5	20.269	10.8390374331551
-14	-7.5	19.45	10.4010695187166
-15	-8.5	19.176	10.2545454545455
-16	-9.5	19.418	10.3839572192513
-17	-10.5	20.348	10.8812834224599
-18	-11.5	21.659	11.5823529411765
-19	-12.5	26.2759	14.0512834224599
-20	-13.5	30.733	16.4347593582888
-21	-14.5	37.472	20.0385026737968
-22	-15.5	40.177	21.4850267379679
-23	-16.5	46.1096	24.6575401069519
1	7.5	75.545	40.3983957219251
2	8.5	70.22	37.5508021390374
3	9.5	65.789	35.1812834224599
4	10.5	61.97	33.1390374331551
5	11.5	58.02	31.0267379679144
6	12.5	53.13	28.4117647058824
7	13.5	47.635	25.4732620320856
8	14.5	41.437	22.1588235294118
9	15.5	37.118	19.8491978609626
10	16.5	32.67	17.4705882352941
11	17.5	27.595	14.7566844919786
12	18.5	22.228	11.8866310160428
13	19.5	20.652	11.0438502673797
14	20.5	19.76	10.5668449197861
15	21.5	19.537	10.4475935828877
16	22.5	19.69	10.5294117647059
17	23.5	20.3	10.855614973262
18	24.5	22.134	11.8363636363636
19	25.5	25.69	13.7379679144385
20	26.5	31.51	16.8502673796791
21	27.5	36.24	19.379679144385
22	28.5	41.11	21.9839572192513
23	29.5	46.32	24.7700534759358
24	30.5	50.1	26.7914438502674
25	31.5	55.05	29.4385026737968
26	32.5	60.125	32.1524064171123
27	33.5	64.6	34.5454545454545
28	34.5	69.6	37.2192513368984
29	35.5	74.337	39.7524064171123
30	36.5	76.91	41.1283422459893
31	37.5	75.366	40.3026737967914
32	38.5	72.25	38.6363636363636
33	39.5	67.42	36.0534759358289
34	40.5	62.12	33.2192513368984
35	41.5	57.4	30.6951871657754
36	42.5	53.19	28.4438502673797
37	43.5	46.95	25.1069518716578
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
e, g = 1.87, 100E-6

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
