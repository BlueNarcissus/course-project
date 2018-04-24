import scipy.integrate as integrate
import scipy.special as special
import numpy as np

""" double quad """
# inner integrals need to be defined as functions
def I(n):
    return integrate.dblquad(lambda t, x: np.exp(-x*t)/t**n, 0, np.inf, lambda t:1, lambda t:np.inf)

print('I(3):', I(3))

# non-constant quad
area = integrate.dblquad(lambda x, y: x*y, 0, 0.5, lambda x:0, lambda x: 1-2*x)
print('Area:', area)



""" n-fold quad """
# arguments from innermost integral to outermost
def J(n):
    return integrate.nquad(lambda t, x: np.exp(-x*t)/t**n, [[1, np.inf], [0, np.inf]])

print('J(3):', J(3))

# non-constant quad
bounds_y = [0, 0.5]
def bounds_x(y):
    return [0, 1-2*y]

area = integrate.nquad(lambda x,y:x*y, [bounds_x, bounds_y])
print('Area:', area)

area = integrate.nquad(lambda x,y:x*y, [lambda x: [0, 1-2*x], [0,0.5]])
print('Area:', area)
