import numpy as np
import ctypes as ct

# NB Most of the ugly details of the function calls and will be hidden from
# the "end user" later on

# open the shared library in question
fortlib = ct.CDLL('demo_lib.so')
# define the argument types of the function
fortlib.demo_sub.argtypes = [ct.POINTER(ct.c_double), ct.c_int64, ct.c_int64]
# define the argument types of the function
fortlib.set_mod_arr.argtypes = [ct.POINTER(ct.c_double), ct.c_int64, ct.c_int64]
# define the argument types of the function
fortlib.get_mod_arr.argtypes = [ct.POINTER(ct.c_double), ct.c_int64, ct.c_int64]

# demo array to call the function with
arr = np.zeros((2,3), order="F")

# call the function; it should print a few things
fortlib.demo_sub(arr.ctypes.data_as(ct.POINTER(ct.c_double)), arr.shape[0], arr.shape[1])

fortlib.set_mod_arr(arr.ctypes.data_as(ct.POINTER(ct.c_double)), arr.shape[0], arr.shape[1])

arr[()]=42

fortlib.get_mod_arr(arr.ctypes.data_as(ct.POINTER(ct.c_double)), arr.shape[0], arr.shape[1])

# demonstrate that the function actually changed the array contents
print(arr)
