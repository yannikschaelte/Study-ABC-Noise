from ..vars import ModelVars
import numpy as np


class Gaussian1DModelVars(ModelVars):

    def __init__(prior_lb=0, prior_ub=5, noise_std=0.05):
        super().__init__(p_true = {'p0': 2.5})
        self.prior_lb = prior_lb
        self.prior_ub = prior_ub
        self.noise_std = noise_std
        self.noise_model = np.random.randn

    def get_id(self):
        return f"gaussian1d_{noise_std}"

    def get_prior(self):
        return pyabc.Distribution(
            **{key: pyabc.RV('uniform', self.prior_lb,
                             self.prior_ub - self.prior_lb)
               for key in self.p})
    
    def get_distance(self):
        return pyabc.PNormDistance(p=2)

    def get_kernel(self):
        return pyabc.distance.IndependentNormalKernel(
            mean=[0], var=[self.noise_std**2])
    
    def get_model(self):
        def model(p):
            return {'y0': p['p0']}
        return model

    def get_model_noisy(self):
        def model_noisy(p):
            return {'y0': p['p0'] + self.noise_std * self.noise_model()}
        return model_noisy

    def generate_data(self):
        return call_noisy(self.p_true)
