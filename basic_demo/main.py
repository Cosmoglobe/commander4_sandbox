import numpy as np
import ctypes as ct

# open the shared library in question
fortlib = ct.CDLL('demo_lib.so')
# define the arguent types of the function
fortlib.demo_sub.argtypes = [ct.POINTER(ct.c_double), ct.c_int64, ct.c_int64]

# demo array to call the function with
arr = np.zeros((2,3), order="F")

# call the function; it should print a few things
fortlib.demo_sub(arr.ctypes.data_as(ct.POINTER(ct.c_double)), arr.shape[0], arr.shape[1])

# demonstrate that the function actually changed the array contents
print(arr)
