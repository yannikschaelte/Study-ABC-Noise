from .base import Model
import numpy as np


class ConversionReactionModel(Model):

    def __init__(self):
        super().__init__()
        self.limits = {'p0': (0, 0.4), 'p1': (0, 0.4)}
        self.noise_std = 0.02
        self.noise_model = np.random.randn
        self.p_true = {'p0': 0.06, 'p1': 0.08}
        self.ts = np.linspace(0, 30, n_t)
        self.x0 = np.array([1., 0.])

    def get_prior(self):
        return pyabc.Distribution(
            **{key: pyabc.RV('uniform', bounds[0], bounds[1])
               for key, bounds in self.limits.items()})
    
    def get_distance(self):
        def l2(x, y):
            return np.sum(np.power( (x['y'] - y['y']) / self.noise_std, 2))
        return l2

    def call(self, p):
        y = x(p, self.x0, self.ts)[1, :]
        return {'y': y.flatten()}
    
    def call_noisy(self, p):
        y = x(p, self.x0, self.ts)[1, :] \
            + self.noise_std * self.noise_model(1, len(self.ts))
        return {'y': y.flatten()}

def x(p, x0, ts):
    """
    States via analytic solution of ODE.
    Returns an array of shape n_x * n_t.
    """
    p0 = p['p0']
    p1 = p['p1']
    n_t = len(ts)
    sol = np.zeros((2, n_t))
    for ix, t in enumerate(ts):
        e = np.exp(- (p0 + p1) * t)
        A = 1 / (- p0 - p1) * np.array([[- p1 - p0 * e, - p1 + p1 * e],
                                        [- p0 + p0 * e, - p0 - p1 * e]])
        sol[:, ix] = np.dot(A, x0).flatten()
