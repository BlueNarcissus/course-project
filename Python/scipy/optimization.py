import numpy as np
from scipy.optimize import minimize

""" nelder-mead method """
def rosenbrock(x):
    return sum(100*(x[1:]-x[:-1]**2)**2 + (1-x[:-1])**2)

x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2])

print('nelder-mead ...')
result = minimize(rosenbrock, x0, method='nelder-mead', options={'xtol':1e-8, 'disp': True})
print(result.x)


""" BFGS method """
# function to compute derivative of rosenbrock function
def rosen_der(x):
    der = np.zeros_like(x)
    xm = x[1:-1]
    xm_p1 = x[2:]
    xm_m1 = x[:-2]
    der[1:-1] = 200*(xm-xm_m1**2) - 400*xm*(xm_p1-xm**2) - 2*(1-xm)
    der[0] = -400*x[0]*(x[1]-x[0]**2) - 2*(1-x[0])
    der[-1] = 200*(x[-1]-x[-2]**2)
    return der

print('\nBFGS ...')
result = minimize(rosenbrock, x0, method='BFGS', jac=rosen_der, options={'disp': True})
print(result.x)


""" Newton-CG """
def rosen_hess(x):
    diagonal = np.zeros_like(x)
    diagonal[0] = 2 - 400*x[1] + 1200*x[0]**2
    diagonal[-1] = 200
    diagonal[1:-1] = 202 - 400*x[2:] + 1200*x[1:-1]**2
    
    H = np.diag(-400*x[:-1], 1) + np.diag(-400*x[:-1], -1)
    H += np.diag(diagonal)
    return H

print('\nNewton-CG ...')
result = minimize(rosenbrock, x0, method='Newton-CG', jac=rosen_der, hess=rosen_hess,
                  options={'xtol':1e-8, 'disp': True})
print(result.x)


""" Newton-CG, Hessian product """
def rosen_hess_prod(x, p):
    Hp = np.zeros_like(x)
    Hp[1:-1] = -400*x[:-2]*p[:-2] - 400*x[1:-1]*p[2:] + (202+1200*x[1:-1]**2-400*x[2:])*p[1:-1]
    Hp[0] = (1200*x[0]**2 - 400*x[1] + 2)*p[0] - 400*x[0]*p[1]
    Hp[-1] = -400*x[-2]*p[-2]+200*p[-1]
    return Hp

print('\nHessian product ...')
result = minimize(rosenbrock, x0, method='Newton-CG', jac=rosen_der, hessp=rosen_hess_prod,
                  options={'xtol':1e-8, 'disp': True})
print(result.x)
