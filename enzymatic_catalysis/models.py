import pyabc
import numpy as np
import scipy as sp

# VARIABLES

# MODEL

n_timepoints = 50
timepoints = np.linspace(0, 5, n_timepoints)
y0 = 


def f(x, t, th0, th1, th2, th3):
    x0, x1, x2, x3 = x
    dx0 = - th0*x0*x1 + th1*x2
    dx1 = - th0*x0*x1 + (th1+th2)*x2 - th3*x1*x3
    dx2 = th0*x0*x1 - (th1+th2)*x2 + th3*x1*x3
    dx3 = th2*x2 - th3*x1*x3
    return dx0, dx1, dx2, dx3


def y(p):
    th0 = p['th0']
    th1 = p['th1']
    th2 = p['th2']
    th3 = p['th3']
    sol = sp.integrate.odeint(f, y0, timepoints, args=(th0, th1, th2, th3))
