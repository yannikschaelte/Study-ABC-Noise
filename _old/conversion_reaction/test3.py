"""
This is for testing the stochastic acceptor.
"""

import pyabc
import numpy as np
import scipy as sp
from models import *
from scipy import stats

# VARIABLES

db_path = "sqlite:///db3.db"

# ACCEPTOR

def pdf(x_0, x):
    v_0 = np.array(list(x_0.values()))
    v = np.array(list(x.values()))
    mean = np.zeros(n_timepoints)
    cov = noise**2 * np.eye(n_timepoints)
    return stats.multivariate_normal.pdf(v_0 - v, mean=mean, cov=cov)

acceptor = pyabc.StochasticAcceptor(pdf=pdf, temp_decay_exp=4, target_acceptance_rate=0.5, use_target_acceptance_rate = False)

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior,
                   distance_function=pyabc.NoDistance(),
                   population_size=pop_size,
                   transitions=transition,
                   eps=pyabc.NoEpsilon(),
                   acceptor=acceptor,
                   sampler=sampler)

abc.new(db_path, y_obs)
h = abc.run(minimum_epsilon=1, max_nr_populations=max_nr_populations, min_acceptance_rate=min_acceptance_rate)
h = pyabc.History(db_path)

# PLOT
#viz("test3", h)
