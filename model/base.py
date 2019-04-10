import pyabc
from abc import ABC

class Model(ABC):

    def __init__(self):
        self.pop_size = 1000
        self.n_pop = 20
        self.eps_min = 0.0
        self.min_acc_rate = 0.0

    def get_p_true(self):
        return None

    def get_prior(self):
        raise NotImplementedError()

    def get_transition(self):
        return pyabc.MultivariateNormalTransition()

    def get_distance(self):
        """
        Distance to use for deterministic acceptance.
        """
        return pyabc.PNormDistance(p=2)

    def get_kernel(self):
        """
        Kernel to use as for stochastic acceptance.
        """
        raise NotImplementedError()

    def get_eps(self):
        return pyabc.MedianEpsilon()

    def get_id(self):
        raise NotImplementedError()

    def call(self):
        raise NotImplementedError()

    def call_noisy(self):
        raise NotImplementedError()

    def get_true_unnormalized_pdf_fun(y_obs):
        raise NotImplementedError()
