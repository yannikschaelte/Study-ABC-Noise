from .base import Model
import numpy as np


class Gaussian1DModel(Model):

    def __init__(prior_lb=0, prior_ub=5, noise_std=0.05):
        self.prior_lb = prior_lb
        self.prior_ub = prior_ub
        self.noise_std = noise_std
        self.noise_model = np.random.randn

    def get_prior(self):
        return pyabc.Distribution(
            **{key: pyabc.RV('uniform', self.prior_lb,
                             self.prior_ub - self.prior_lb)
               for key in self.p})
    
    def get_id(self):
        return f"gaussian1d_{noise_std}"

    def call(self, p):
        return {'y0': p['p0']}

    def call_noisy(self, p):
        return {'y0': p['p0'] + self.noise_std * self.noise_model()}
