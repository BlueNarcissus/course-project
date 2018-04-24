import scipy.integrate as integrate
import ctypes
import time

# c library for function
lib = ctypes.CDLL('/Users/liangjin/Desktop/computer/2_python/project/Scipy/integrate/ctype/testlib.so')
func = lib.func
func.restype = ctypes.c_double
func.argtypes = (ctypes.c_int, ctypes.c_double)

tic = time.time()
result = integrate.nquad(func, [[0,10], [-10,0], [-1,1]])
toc = time.time()
print(result)
print('time elapsed: %f s' %(toc-tic))

# normal lambda function
tic = time.time()
f = lambda x0, x1, x2: x0-x1*x2
result = integrate.nquad(f, [[0,10], [-10,0], [-1,1]])
toc = time.time()
print(result)
print('time elapsed: %f s' %(toc-tic))

