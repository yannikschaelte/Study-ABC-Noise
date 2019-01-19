"""
This is for testing the stochastic acceptor.
"""

import pyabc
import numpy as np
import scipy as sp
from models import *
from scipy import stats

# VARIABLES

db_path = "sqlite:///db4.db"

# ACCEPTOR

distr = stats.multivariate_normal(np.zeros(n_timepoints), noise**2 * np.eye(n_timepoints))
nr_pops = 1
acceptor = pyabc.StochasticAcceptor(distribution=distr, nr_populations=nr_pops)
acceptor.max_temp = 200

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior,
                   distance_function=distance_l2,
                   population_size=pop_size,
                   transitions=transition,
                   eps=eps,
                   acceptor=acceptor,
                   sampler=sampler)

abc.new(db_path, y_obs)
h = abc.run(minimum_epsilon=0, max_nr_populations=nr_pops, min_acceptance_rate=min_acceptance_rate)
h = pyabc.History(db_path)

# PLOT
#viz("test3", h)
