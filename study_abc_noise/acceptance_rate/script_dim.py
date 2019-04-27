"""
Compare the behavior of the accepatance rate with the data dimension.
"""

from study_abc_noise.model.conversion_reaction import *
from study_abc_noise.model import ConversionReactionModelVars as ModelVars
from study_abc_noise.vars import get_data
from study_abc_noise.optimize import get_optimal_kernel_value
import pyabc
import numpy as np
import scipy as sp
import cloudpickle as pickle


list_n_t = [4, 7, 10, 13, 16, 19, 22]
n_data = 1

list_pdf_max = []
list_pdf_max_min = []
list_acc_rate_prior = []
list_acc_rate_prior_min = []
list_acc_rate_posterior = []
list_acc_rate_posterior_min = []

for n_t in list_n_t:
    for i_data in range(n_data):
        mv = ModelVars(n_t=n_t)
        # dummy task
        #task = Task.from_vars(AnalysisVars(), mv)
        limits = mv.limits
        # generate data
        y_obs = get_data(mv, i_data)
        acc_prob_integrand_from_prior = get_acceptance_probability_integrand_from_prior(
            mv, y_obs, 1.0)

        pdf_max = normal_dty(y_obs['y'].flatten(), y_obs['y'].flatten(), mv.noise_std * np.ones(mv.n_t))
        list_pdf_max.append(pdf_max)

        pdf_max_min = np.exp(get_optimal_kernel_value(mv, y_obs)[1])
        list_pdf_max_min.append(pdf_max_min)

        print("pdf_max: ", pdf_max, pdf_max_min)

        posterior_scaled = get_posterior_scaled(mv, y_obs)
        acc_prob_integrand_from_posterior = get_acceptance_probability_integrand_from_posterior(
            mv, y_obs, posterior_scaled, 1.0)

        integral_prior = sp.integrate.dblquad(
            lambda x, y: acc_prob_integrand_from_prior([x, y]),
            limits['p0'][0], limits['p0'][1],
            lambda x: limits['p1'][0], lambda x: limits['p1'][1])[0]
        integral_posterior = sp.integrate.dblquad(
            lambda x, y: acc_prob_integrand_from_posterior([x, y]),
            limits['p0'][0], limits['p0'][1],
            lambda x: limits['p1'][0], lambda x: limits['p1'][1])[0]
        
        acc_rate_prior = integral_prior / pdf_max
        acc_rate_prior_min = integral_prior / pdf_max_min
        acc_rate_posterior = integral_posterior / pdf_max
        acc_rate_posterior_min = integral_posterior / pdf_max_min

        list_acc_rate_prior.append(acc_rate_prior)
        list_acc_rate_prior_min.append(acc_rate_prior_min)
        list_acc_rate_posterior.append(acc_rate_posterior)
        list_acc_rate_posterior_min.append(acc_rate_posterior_min)

        print(acc_rate_prior, acc_rate_prior_min, acc_rate_posterior, acc_rate_posterior_min)

with open("dim.dat", 'wb') as f:
    pickle.dump([list_pdf_max, list_pdf_max_min,
                 list_acc_rate_prior, list_acc_rate_prior_min,
                 list_acc_rate_posterior, list_acc_rate_posterior_min], f)
