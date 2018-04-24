import scipy.integrate as integrate
import scipy.special as special
import numpy as np

# gassian quadrature ???
func = lambda x: np.exp(-x)/x**2
result = integrate.quadrature(func, 1, 5)
print(result)

# Romberg integration
result = integrate.romberg(func, 1, 5)
print(result)
