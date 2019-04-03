from .base import Model


class GaussianModel(Model):

    def __init__(prior_lb, prior_ub, noise_std):
        self.prior_lb = prior_lb
        self.prior_ub = prior_ub
        self.noise_std = noise_std
        self.noise_model = np.random.randn

    def prior(self):
        return pyabc.Distribution(
            **{key: pyabc.RV('uniform', self.prior_lb,
                             self.prior_ub - self.prior_lb)
               for key in self.p})
 
    def call(p):
        return {'y0': p['p0']}

    def call_noisy(p):
        return {'y0': p['p0'] + self.noise_std * self.noise_model()}


