from ..vars import ModelVars
import numpy as np
import pyabc
import matplotlib.pyplot as plt


class NonIdAbleModelVars(ModelVars):

    def __init__(self):
        super().__init__(p_true = {'p0': 0.4, 'p1': 0.5})
        self.limits = {'p0': (0, 1), 'p1': (0, 1)}
        self.noise_std = 0.2
        self.noise_model = np.random.randn
        self.n_t = 5
        self.t_max = 30
        self.x0 = 1.0

    def get_ts(self):
        return np.linspace(0, self.t_max, self.n_t)

    def get_id(self):
        return f"nia_{self.n_t}_{self.noise_std}"

    def get_prior(self):
        return pyabc.Distribution(
            **{key: pyabc.RV('uniform', bounds[0], bounds[1])
               for key, bounds in self.limits.items()})

    def get_distance(self):
        def l2(x, y):
            return np.sum(np.power( (x['y'] - y['y']) / self.noise_std, 2))
        return l2

    def get_kernel(self):
        kernel = pyabc.distance.IndependentNormalKernel(
            mean=np.zeros(self.n_t),
            var=self.noise_std**2 * np.ones(self.n_t))
        return kernel

    def get_model(self):
        def model(p):
            y = x(p, self.x0, self.get_ts())
            return {'y': y.flatten()}
        return model

    def get_model_noisy(self):
        def model_noisy(p):
            y = x(p, self.x0, self.get_ts()) \
                + self.noise_std * self.noise_model(1, self.n_t)
            return {'y': y.flatten()}
        return model_noisy

    def generate_data(self):
        y = self.get_model_noisy()(self.p_true)
        return y

    def viz_x(self, ax=None):
        if ax is None:
            _, ax = plt.subplots()
        ts = self.get_ts()
        xs = x(self.p_true, self.x0, ts)
        ax.plot(ts, xs, 'x-', color='C2', label="Species A")
        ax.set_xlabel("Time [au]")
        ax.set_ylabel("Concentration [au]")
        ax.legend()
        return ax

    def viz_y(self, y, label="", ax=None):
        if ax is None:
            _, ax = plt.subplots()
        ax.plot(self.get_ts(), y['y'], 'x-', color='C2', label="Species A")
        ax.set_xlabel("Time [au]")
        ax.set_ylabel("Concentration [au]")
        ax.legend()
        return ax

    def viz_data_and_sim(self, y, label="", ax=None):
        if ax is None:
            _, ax = plt.subplots()
        ts = self.get_ts()
        y_true = self.get_model()(self.p_true)
        ax.plot(ts, y_true['y'], 'x-', color='C0', label="Observable")
        ax.plot(ts, y['y'], 'x-', color='C2', label="Measured data")
        ax.set_xlabel("Time [au]")
        ax.set_ylabel("Concentration [au]")
        ax.legend()
        return ax


class NonIdAblePrioredModelVars(NonIdAbleModelVars):

    def get_id(self):
        return f"niap_{self.n_t}_{self.noise_std}"

    def get_prior(self):
        return pyabc.Distribution(
            **{key: pyabc.RV('norm', 0.5, 0.4)
               for key in self.limits})


def x(p, x0, ts):
    """
    States via analytic solution of ODE.
    Returns an array of shape n_x * n_t.
    """
    p0 = p['p0']
    p1 = p['p1']
    n_t = len(ts)
    sol = np.zeros((n_t))
    for ix, t in enumerate(ts):
        sol[ix] = np.exp( (p0 - p1) * t ) * x0
    return sol
