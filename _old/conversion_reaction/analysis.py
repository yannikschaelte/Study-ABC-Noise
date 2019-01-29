import numpy as np
import scipy as sp
import pickle
from models import *


print(y_obs)

def integrand_from_prior(th, c):
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

    return (likelihood_val / c) * prior_val


c_max = normal_dty(get_y_obs()['y'].flatten(), get_y_obs()['y'].flatten(), noise * np.ones(n_timepoints))
def integrand_from_prior_c_max(th):
    return integrand_from_prior(th, c_max)
c_found = pickle.load(open("best_found_pdf_" + str(noise) + "_" + str(n_timepoints) + ".dat", 'rb'))
def integrand_from_prior_c_found(th):
    return integrand_from_prior(th, c_found)

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


def integrand_from_posterior(th, c):
    if type(th) is list:
        th = {'th0': th[0], 'th1': th[1]}

    # data 
    y = model(th)['y'].flatten()
    y_obs = get_y_obs()['y'].flatten()
    
    sigma = noise * np.ones(len(y))

    # likelihood
    likelihood_val = normal_dty(y_obs, y, sigma)

    # posterior
    posterior_val = true_posterior(th)

    return (likelihood_val / c)  * posterior_val

def integrand_from_posterior_c_max(th):
    return integrand_from_posterior(th, c_max)
def integrand_from_posterior_c_found(th):
    return integrand_from_posterior(th, c_found)


integral_prior_max = sp.integrate.dblquad(lambda x, y: integrand_from_prior_c_max([x,y]), limits['th0'][0], limits['th0'][1], lambda x: limits['th1'][0], lambda x: limits['th1'][1])[0]
integral_posterior_max = sp.integrate.dblquad(lambda x, y: integrand_from_posterior_c_max([x,y]), limits['th0'][0], limits['th0'][1], lambda x: limits['th1'][0], lambda x: limits['th1'][1])[0]


integral_prior_found = sp.integrate.dblquad(lambda x, y: integrand_from_prior_c_found([x,y]), limits['th0'][0], limits['th0'][1], lambda x: limits['th1'][0], lambda x: limits['th1'][1])[0]
integral_posterior_found = sp.integrate.dblquad(lambda x, y: integrand_from_posterior_c_found([x,y]), limits['th0'][0], limits['th0'][1], lambda x: limits['th1'][0], lambda x: limits['th1'][1])[0]

print("prior_c_max, posterior_c_max, prior_c_found, posterior_c_found", integral_prior_max, integral_posterior_max, integral_prior_found, integral_posterior_found)


pickle.dump((integral_prior_max, integral_posterior_max, integral_prior_found, integral_posterior_found), open("analysis.dat", 'wb'))
