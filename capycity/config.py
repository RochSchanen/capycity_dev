# file: config.py

# from package: "https://numpy.org/"
from numpy import array

# define data type
from numpy import uint32 as DATATYPE
""" the precision does not influence the error when storing the data.
the idea that the result precision is 1/256 for the uint8 data type
turns out to be naive. Even when the number of iteration has defintely
passed the test of stability for a laplace solution. The final result
can still be largely influenced by the integer precision. This can
easily be checked empirically. As it turns out, a size of 32 bits seems
to be the minimum required to provide a reliable solution. (This is
higher than I personaly expected. I thought uint16 would be largely
sufficient but it is not...). For development purposes, using 8 bits
is a good ideas when checking operations on arrays. just uncomment the
line next to use 8 bits range data type """
# from numpy import uint8 as DATATYPE

# define ZER and MAX values
ZER, MAX = array([0, -1], DATATYPE)
""" ZER is always represented by a chain of "zero" bits.
MAX is always represented by a chain of "one" bits.
The size of the chain depends on the integer size. """
