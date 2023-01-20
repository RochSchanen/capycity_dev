

from numpy import arange

print(arange(0.10, 0.45, 0.025))

print([100]*10)

############################################################# PARTIAL AVERAGING

# from numpy import arange, vsplit, hsplit

# A = arange(0, 16).reshape(4, 4)
# print(A)

# V = vsplit(A, 2)
# print(V)
# print()

# H = hsplit(A, 2)
# print(H)
# print()

############################################################# PARTIAL AVERAGING

# from numpy import array, average

# A = array([
# 	[1, 2, 3, 4, 5, 6],
# 	[1, 2, 3, 4, 5, 6],
# 	[1, 2, 3, 4, 5, 6],
# 	[1, 2, 3, 4, 5, 6],
# 	[1, 2, 3, 4, 5, 6],
# 	[1, 2, 3, 4, 5, 6],
# 	])

# print(f"orginal array:\n{A}\n")

# B = A.reshape(6, 3, 2)
# print(f"split array:\n{B}\n")

# C = average(B, axis = 2)
# print(f"averaged lines:\n{C}\n")

############################################################ CLASS DICTIONARIES
	
# myclasscounter = 0

# class myclass:
	
# 	""" text """
	
# 	n = 0
# 	L = []

# 	def __init__(self):
# 		print(f"create '{__class__.__name__}'")
# 		self.L.append("hello")
# 		self.n += 1
# 		self.m = 99
# 		global myclasscounter
# 		myclasscounter += 1
# 		self.MYNAME = f"{__class__.__name__}_{myclasscounter}"
# 		print(f"__class__.__dict__ :\n{__class__.__dict__}\n")
# 		print(f"instance.__dict__ :\n{self.__dict__}\n")
# 		return

# m = myclass()
# print(myclasscounter)

# m.n = 100
# m.L.append("world")

# n = myclass()

# print(m.n, n.n)

# print(myclasscounter)


# !!!!!!!!!!!!!! note that the second instance inherit the content of the
# list self.L from the first instance. but not the value of self.n
# is there any logic to that...

########################################################## PATCH VERTICAL LINES

# def section_decor(self):

#     l = vlines(
#         [-self.r, +self.r],             # horizontal positions
#         0.0, 1.0,                       # vertical span
#         edgecolor = (0.7, 0.7, 0.7),    # color
#         linestyle = "--",               # line style
#         linewidth = 0.25,               # line thickness
#         )

#     # done
#     return l
