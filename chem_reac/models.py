import pyabc
import numpy as np
import scipy as sp

# VARIABLES

noise = 0.04

# MODEL

n_timepoints = 11
timepoints = sp.arange(n_timepoints)
init = sp.array([1, 0])


def f(x, t0, th0, th1):
    x0, x1 = x
    dx0 = - th0 * x0 + th1 * x1
    dx1 = th0 * x0 - th1 * x1
    return dx0, dx1


def y(p):
    th0 = p['th0']
    th1 = p['th1']
    sol = sp.integrate.odeint(f, init, timepoints, args=(th0, th1))
    return sol


def model(p):
    y0 = y(p)[:, 1]
    return {'y0': y0} 


def model_random(p):
    y0 = y(p)[:, 1] + noise * np.random.rand(n_timepoints)
    return {'y0': y0}


def distance(x, y):
    return sp.absolute(x['y0'] - y['y0']).sum()


class ArrayPNormDistance(pyabc.PNormDistance):

    def __init__(self):
        super().__init__(p=1)

    def initialize(self, t, sample_from_prior):
        sum_stats = []
        for sum_stat in sample_from_prior:
            sum_stats.append(normalize_sum_stat(sum_stat))
        super().initialize(t, sum_stats)

    def __call__(self, t, x, y):
        x = normalize_sum_stat(x)
        y = normalize_sum_stat(y)
        return super().__call__(t, x, y)


def normalize_sum_stat(x):
    x_flat = {}
    for key, value in x.items():
        for j in range(len(value)):
            x_flat[(key, j)] = value[j]
    return x_flat


# TRUE VALUES

th0_true, th1_true = sp.exp([-2.5, -2])
th_true = {'th0': th0_true, 'th1': th1_true}
data_true = model(th_true)

# MEASURED DATA


def f_data_meas():
    data_meas = np.zeros(n_timepoints)
    for _ in range(1000):
        data_meas += model(th_true)['y0'] + 0.01*np.random.rand(n_timepoints)
    data_meas /= 1000

    return {'y0': data_meas}


data_meas = {'y0': data_true['y0'] + noise * np.random.rand(n_timepoints)}

