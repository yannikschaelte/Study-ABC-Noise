from ..vars import ModelVars
import numpy as np
import scipy as sp
import pyabc
import matplotlib.pyplot as plt


P0_TRUE = 0.06
P1_TRUE = 0.08


class ConversionReactionModelVars(ModelVars):

    def __init__(self, p_true = None, pdf_max = None, n_t: int = 10, t_max: float = 30,
                 n_pop: int = None, n_acc: int = None):
        if p_true is None:
            p_true = {'p0': P0_TRUE, 'p1': P1_TRUE}
        super().__init__(p_true = p_true, pdf_max = pdf_max, n_pop = n_pop, n_acc = n_acc)
        self.limits = {'p0': (0, 0.4), 'p1': (0, 0.4)}
        self.noise_std = 0.02
        self.noise_model = np.random.randn
        self.n_t = n_t
        self.t_max = t_max
        self.x0 = np.array([1., 0.])

    def get_ts(self):
        return np.linspace(0, self.t_max, self.n_t)

    def get_id(self):
        return f"cr_{self.n_t}_{self.noise_std}"

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
            var=self.noise_std**2 * np.ones(self.n_t),
            pdf_max=self.pdf_max)
        return kernel

    def get_model(self):
        def model(p):
            y = x(p, self.x0, self.get_ts())[1, :]
            return {'y': y.flatten()}
        return model

    def get_model_noisy(self):
        def model_noisy(p):
            y = x(p, self.x0, self.get_ts())[1, :] \
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
        ax.plot(ts, xs[0, :], 'x-', color='C1', label="Species A")
        ax.plot(ts, xs[1, :], 'x-', color='C2', label="Species B")
        ax.set_xlabel("Time [au]")
        ax.set_ylabel("Concentration [au]")
        ax.legend()
        return ax

    def viz_y(self, y, label="", ax=None):
        if ax is None:
            _, ax = plt.subplots()
        ax.plot(self.get_ts(), y['y'], 'x-', color='C2', label="Species B")
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


class ConversionReactionLaplaceModelVars(ConversionReactionModelVars):

    def get_distance(self):
        def l1(x, y):
            return np.sum(np.abs( (x['y'] - y['y']) / self.noise_std))
        return l2

    def get_kernel(self):
        kernel = pyabc.distance.IndependentLaplaceKernel(
            scale=self.noise_std * np.ones(self.n_t),
            pdf_max=self.pdf_max)
        return kernel

    def get_model(self):
        def model(p):
            y = x(p, self.x0, self.get_ts())[1, :]
            return {'y': y.flatten()}
        return model
    
    def get_model_noisy(self):
        def model_noisy(p):
            y = x(p, self.x0, self.get_ts())[1, :] \
                + np.random.laplace(0, self.noise_std, self.n_t)
            return {'y': y.flatten()}
        return model_noisy

    def generate_data(self):
        y = self.get_model_noisy()(self.p_true)
        return y


class ConversionReaction1dModelVars(ConversionReactionModelVars):

    def __init__(self, n_acc = 10000, p_true = None, pdf_max = None,
                 n_t: int = 10, t_max: float = 30,
                 n_pop: int = 15):
        if p_true is None:
            p_true = {'p0': P0_TRUE}
        super().__init__(p_true = p_true, pdf_max = pdf_max, n_pop = n_pop, n_acc = n_acc)
        # just assume a narrow prior
        self.limits = {'p0': (0.053, 0.066)}


class ConversionReactionUVarModelVars(ConversionReactionModelVars):

    def __init__(self, n_acc=1000, pdf_max=None):
        super().__init__(p_true = {'p0': 0.06, 'p1': 0.08, 'std': 0.02},
                         pdf_max = pdf_max, n_acc=n_acc)
        self.limits = {'p0': (0, 0.15), 'p1': (0, 0.15), 'std': (0.01, 0.05)}
        # used in distance, just for normalization
        self.noise_std = self.p_true['std']
        self.noise_model = np.random.randn
        self.n_t = 10
        self.t_max = 30
        self.x0 = np.array([1., 0.])

    def get_id(self):
        return f"cr_uvar_{self.n_t}_{self.noise_std}"

    def get_model_noisy(self):
        def model_noisy(p):
            y = x(p, self.x0, self.get_ts())[1, :] \
                + p['std'] * self.noise_model(1, self.n_t)
            return {'y': y.flatten()}
        return model_noisy

    def get_kernel(self):
        def compute_var(p):
            return p['std']**2 * np.ones(self.n_t)
        kernel = pyabc.distance.IndependentNormalKernel(
            var=compute_var,
            pdf_max=self.get_pdf_max()
        )
        return kernel
    
    def get_pdf_max(self):
        if self.pdf_max is not None:
            return self.pdf_max
        return - 0.5 * self.n_t * np.log(2 * np.pi * self.limits['std'][0]**2)

def x(p, x0, ts):
    """
    States via analytic solution of ODE.
    Returns an array of shape n_x * n_t.
    """
    p0 = p.get('p0', 0.06)
    p1 = p.get('p1', 0.08)
    n_t = len(ts)
    sol = np.zeros((2, n_t))
    for ix, t in enumerate(ts):
        e = np.exp(- (p0 + p1) * t)
        A = 1 / (- p0 - p1) * np.array([[- p1 - p0 * e, - p1 + p1 * e],
                                        [- p0 + p0 * e, - p0 - p1 * e]])
        sol[:, ix] = np.dot(A, x0).flatten()
    return sol


def normal_dty_1d(y_bar, y, sigma):
    dty = ( 1 / np.sqrt( 2 * np.pi * sigma**2 ) 
            * np.exp( - ( (y_bar - y) / sigma)**2 / 2) )
    return dty


def normal_dty(y_bar, y, sigma):
    """
    Uncorrelated multivariate Gaussian density.

    y_bar: size dim
        point at which to evaluate the density
    y, sigma: size dim
        For N(y, sigma).
    """
    dim = len(y_bar)
    dties = np.zeros(dim)
    for j in range(dim):
        dties[j] = normal_dty_1d(y_bar[j], y[j], sigma[j])
    dty = np.prod(dties)
    return dty


def get_acceptance_probability_integrand_from_prior(model_vars, y_obs, pdf_max=1.0):
    model = model_vars.get_model()
    y_obs = y_obs['y'].flatten()
    prior = model_vars.get_prior()

    def acceptance_probability_integrand_from_prior(p):
        if type(p) is not dict:
            p = {key: p[i] for i, key in enumerate(model_vars.p_true)}

        # data
        y = model(p)['y'].flatten()
        
        sigma = model_vars.noise_std * np.ones(model_vars.n_t)

        # acceptance probability
        likelihood_val = normal_dty(y_obs, y, sigma)
        acceptance_probability = likelihood_val / pdf_max

        # prior
        prior_val = prior.pdf(p)

        # integrand
        integrand = acceptance_probability * prior_val

        return integrand

    return acceptance_probability_integrand_from_prior


def get_posterior_unscaled(model_vars, y_obs):
    model = model_vars.get_model()
    y_obs = y_obs['y'].flatten()
    prior = model_vars.get_prior()

    def posterior_unscaled(p):
        if type(p) is not dict:
            p = {key: p[i] for i, key in enumerate(model_vars.p_true)}

        # data
        y = model(p)['y'].flatten()

        sigma = model_vars.noise_std * np.ones(model_vars.n_t)

        # likelihood
        likelihood_val = normal_dty(y_obs, y, sigma)

        # prior
        prior_val = prior.pdf(p)

        # posterior value
        unscaled_posterior = likelihood_val * prior_val

        return unscaled_posterior

    return posterior_unscaled


def get_posterior_normalization(posterior_unscaled, model_vars):
    limits = model_vars.limits
    posterior_normalization = sp.integrate.dblquad(
        lambda x, y: posterior_unscaled([x, y]),
        limits['p0'][0], limits['p0'][1],
        lambda x: limits['p1'][0], lambda x: limits['p1'][1])[0]
    return posterior_normalization


def get_posterior_normalization_1d(posterior_unscaled, model_vars):
    limits = model_vars.limits
    posterior_normalization = sp.integrate.quad(
        lambda x: posterior_unscaled([x]),
        limits['p0'][0], limits['p0'][1])[0]
    print(f"Normalization: {posterior_normalization}")
    return posterior_normalization


def get_posterior_scaled(model_vars, y_obs):
    posterior_unscaled = get_posterior_unscaled(model_vars, y_obs)
    posterior_normalization = get_posterior_normalization(
        posterior_unscaled, model_vars)
    
    def posterior_scaled(p):
        return posterior_unscaled(p) / posterior_normalization

    return posterior_scaled


def get_posterior_scaled_1d(model_vars, y_obs):
    posterior_unscaled = get_posterior_unscaled(model_vars, y_obs)
    posterior_normalization = get_posterior_normalization_1d(
        posterior_unscaled, model_vars)

    def posterior_scaled(p):
        return posterior_unscaled(p) / posterior_normalization

    return posterior_scaled


def get_acceptance_probability_integrand_from_posterior(
        model_vars, y_obs, posterior_scaled, pdf_max=1.0):
    model = model_vars.get_model()
    y_obs = y_obs['y'].flatten()

    def acceptance_probability_integrand_from_posterior(p):
        if type(p) is not dict:
            p = {key: p[i] for i, key in enumerate(model_vars.p_true)}

        # data
        y = model(p)['y'].flatten()

        sigma = model_vars.noise_std * np.ones(model_vars.n_t)

        # acceptance probability
        likelihood_val = normal_dty(y_obs, y, sigma)
        acceptance_probability = likelihood_val / pdf_max

        # posterior
        posterior_val = posterior_scaled(p)

        # integrand
        integrand = acceptance_probability * posterior_val

        return integrand

    return acceptance_probability_integrand_from_posterior
