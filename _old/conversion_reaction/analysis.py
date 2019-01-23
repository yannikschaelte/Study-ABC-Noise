import numpy as np
import scipy as sp
from models import *


print(y_obs)

def integrand_from_prior(th):
    if type(th) is list:
        th = {'th0': th[0], 'th1': th[1]}

    # data
    y = model(th)['y'].flatten()
    y_obs = get_y_obs()['y'].flatten()

    sigma = noise * np.ones(len(y))

    # likelihood
    likelihood_val = normal_dty(y_obs, y, sigma)

    # c
    c = normal_dty(y_obs, y_obs, sigma)

    # prior
    prior_val = prior.pdf(th)

    return (likelihood_val / c) * prior_val


def true_posterior_unscaled(th):
    if type(th) is list:
        th = {'th0': th[0], 'th1': th[1]}

    # data
    y = model(th)['y'].flatten()
    y_obs = get_y_obs()['y'].flatten()
    
    sigma = noise * np.ones(len(y))

    # likelihood
    likelihood_val = normal_dty(y_obs, y, sigma)

    # prior
    prior_val = prior.pdf(th)

    # posterior value
    unscaled_posterior = likelihood_val * prior_val

    return unscaled_posterior


true_posterior_normalization = integrate.dblquad(lambda x, y: true_posterior_unscaled([x,y]), limits['th0'][0], limits['th0'][1], lambda x: limits['th1'][0], lambda x: limits['th1'][1])[0]


def true_posterior(th):
    return true_posterior_unscaled(th) / true_posterior_normalization


def integrand_from_posterior(th):
    if type(th) is list:
        th = {'th0': th[0], 'th1': th[1]}

    # data 
    y = model(th)['y'].flatten()
    y_obs = get_y_obs()['y'].flatten()
    
    sigma = noise * np.ones(len(y))

    # likelihood
    likelihood_val = normal_dty(y_obs, y, sigma)

    # c
    c = normal_dty(y_obs, y_obs, sigma)

    # posterior
    posterior_val = true_posterior(th)

    return (likelihood_val / c)  * posterior_val


integral_prior= sp.integrate.dblquad(lambda x, y: integrand_from_prior([x,y]), limits['th0'][0], limits['th0'][1], lambda x: limits['th1'][0], lambda x: limits['th1'][1])[0]
integral_posterior = sp.integrate.dblquad(lambda x, y: integrand_from_posterior([x,y]), limits['th0'][0], limits['th0'][1], lambda x: limits['th1'][0], lambda x: limits['th1'][1])[0]

print(integral_prior, integral_posterior)



