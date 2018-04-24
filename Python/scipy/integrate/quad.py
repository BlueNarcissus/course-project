import scipy.integrate as integrate
import scipy.special as special
import numpy as np

# bessel function
result = integrate.quad(lambda x:special.jv(2.5,x), 0, 4.5)
print(result) # first element is integration, second is error.

# polynomial
func = lambda x, a, b: a*x**2+b
a = 2
b = 1
result = integrate.quad(func, 0, 1, args=(a,b))
print(result)


# exponential, infinite limit
# step1: integrate over t
func = lambda t, n, x: np.exp(-t*x)/t**n
expint = lambda n, x: integrate.quad(func, 1, np.inf, args=(n, x))[0]

vec_expint = np.vectorize(expint) # vectorize a function

n = 3
x = np.arange(1.0, 4.0, 0.5)
result = vec_expint(n, x)
print(result)

# evaluate with the predefined function
result_true = special.expn(n, x)
print(result_true)

# step2: integrate over x
func_n = lambda x: expint(n, x)
result_n = integrate.quad(func_n, 0, np.inf)
print(result_n)


















