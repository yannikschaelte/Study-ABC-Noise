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
    true_rate = {'rate': 2.3}
    _obs = None

    def extract_rates(self, par):
        return sp.array([float(par['rate'])])

    def __call__(self, par):
        t, X = gillespie(self.x0,
                         self.extract_rates(par),
                         self.pre, self.post,
                         MAX_T)
        return {'t': t, 'X': X}
    
    def obs(self):
        if self._obs is not None:
            return self._obs

        obs_file = self.__class__.__name__ + ".dat"
        try:
            obs = pickle.load(open(obs_file, 'rb'))
        except Exception:
            obs = self(self.extract_rates(self.true_rate))
            pickle.dump(obs, open(obs_file, 'wb'))

        self._obs = obs
        return obs


class Model2(Model1):
    pre = sp.array([[1, 0]], dtype=int)
    post = sp.array([[0, 1]])


class MRNAModel(Model1):
    """
    MRNA Protein
    """
    x0 = sp.array([0, 0])
    pre = sp.array([[0, 0], [1, 0], [1, 0], [0, 1]], dtype=int)
    post = sp.array([[1, 0], [0, 0], [1, 1], [0, 0]], dtype=int)


class BirthDeathModel(Model1):
    x0 = sp.array([20])
    pre = sp.array([[0], [1]])
    post = sp.array([[1], [0]])
    true_rate = {'r0': 0, 'r1': 0.1}

obs = [Model1()({'rate': true_rate})]
t_test_times = sp.linspace(0, MAX_T, N_TEST_TIMES)


def distance(x, y):
    xt_ind = sp.searchsorted(x['t'], t_test_times) - 1
    yt_ind = sp.searchsorted(y['t'], t_test_times) - 1
    error = (sp.absolute(x['X'][:, 1][xt_ind]
                       - y['X'][:, 1][yt_ind]).sum()
             / t_test_times.size)
    return error


