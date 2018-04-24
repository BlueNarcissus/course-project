from scipy.integrate import odeint
from scipy.special import gamma, airy
import numpy as np

y0 = [-1.0/(3**(1.0/3.0)*gamma(1.0/3.0)), 1.0/(3**(2.0/3.0)*gamma(2.0/3.0))]

def gradient(y, t):
    return [[0,t],[1,0]]

def func(y, t):
    return [t*y[1], y[0]]

x = np.arange(0, 4.0, 0.01)
t = x

# ode integration
y = odeint(func, y0, t) # the derivative function, initial values, variable points
y2 = odeint(func, y0, t, Dfun=gradient)

# exact result, airy function for check
y_check = airy(x)[0]

print('airy: ', y_check[:36:6]) # slice [0:36], pick one number every six number
print('ode: ', y[:36:6, 1])
print('ode: ', y2[:36:6, 1])

