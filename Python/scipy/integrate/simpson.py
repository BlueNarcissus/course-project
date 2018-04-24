import scipy.integrate as integrate
import numpy as np

# simpson's rule to integrate x**2 in range [1,4]
x = np.array([1,3,4])
func_square = lambda x:x**2
func_triple = lambda x:x**3

y1 = func_square(x)
result1 = integrate.simps(y1, x)
print(result1)

# simpson's rule to integrate x**3 in range [1,4]
# not a precise estimate, because order of polynomial is larger than 2.
y2 = func_triple(x)
result2 = integrate.simps(y2, x)
print(result2)


