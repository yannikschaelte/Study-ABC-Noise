from gillespie import gillespie
import numpy as np
import scipy as sp

# CONSTANTS

MAX_T = 0.1
N_TEST_TIMES = 20


class Model1:
    x0 = sp.array([40, 3])
    pre = sp.array([[1, 1]], dtype=int)
    post = sp.array([0, 2])

    def __call__(self, par):
        t, X = gillespie(self.x0,
                         sp.array([float(par['rate'])]),
                         self.pre, self.post,
                         MAX_T)
        return {'t': t, 'X': X}


class Model2(Model1):
    pre = sp.array([[1, 0]], dtype=int)
    post = sp.array([[0, 1]])


true_rate = 2.3
obs = [Model1()({'rate': true_rate})]
t_test_times = sp.linspace(0, MAX_T, N_TEST_TIMES)

def distance(x, y):
    xt_ind = sp.searchsorted(x['t'], t_test_times) - 1
    yt_ind = sp.searchsorted(y['t'], t_test_times) - 1
    error = (sp.absolute(x['X'][:, 1][xt_ind]
                       - y['X'][:, 1][yt_ind]).sum()
             / t_test_times.size)
    return error
