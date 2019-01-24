import numpy as np
import scipy as sp
from models import *


def pdf(x_0, x):
    return stats.multivariate_normal.pdf(np.array(list(x_0.values())) - np.array(list(x.values())), mean=[0], cov=[noise**2])


noise = 0.01


def integrand_from_prior(th):
    y = model({'th0': th})['y0']
    y_obs = get_y_meas()['y0']
    return 1 / (prior_ub - prior_lb) * np.exp(- (y - y_obs)**2 / (2 * noise**2))


def true_posterior_unscaled(th):
    def uniform_dty(th):
        if th < prior_lb or th > prior_ub:
            return 0
        else:
            return 1 / (prior_ub - prior_lb)

    def normal_dty(y, y_obs):
        return 1 / np.sqrt(2*np.pi*noise**2) \
               * np.exp( - (y - y_obs)**2 / (2 * noise**2) )

    return normal_dty(model({'th0': th})['y0'], get_y_meas()['y0']) * uniform_dty(th)


true_posterior_normalization = integrate.quad(true_posterior_unscaled, prior_lb, prior_ub)[0]


def true_posterior(th):
    return true_posterior_unscaled(th) / true_posterior_normalization


def integrand_from_posterior(th):
    y = model({'th0': th})['y0']
    y_obs = get_y_meas()['y0']
    return np.exp(- (y - y_obs)**2 / (2 * noise**2)) * true_posterior(th)


integral_prior= sp.integrate.quad(integrand_from_prior, prior_lb, prior_ub)[0]
integral_posterior = sp.integrate.quad(integrand_from_posterior, prior_lb, prior_ub)[0]

print(integral_prior, integral_posterior)



