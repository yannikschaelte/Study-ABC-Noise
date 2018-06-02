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

distr = stats.multivariate_normal(np.zeros(n_timepoints), noise**2 * np.eye(n_timepoints))
acceptor = pyabc.StochasticAcceptor(distr=distr)

# PERFORM ABC ANALYSIS WITH 0,1-THRESHOLD AS PREPARATION

max_nr_populations = max_nr_populations - 1

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior,
                   distance_function=distance,
                   population_size=pop_size,
                   transitions=transition,
                   eps=eps,
                   acceptor=acceptor,
                   sampler=sampler)

abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations)

# PERFORM ABC ANALYSIS WITH STOCHASTIC ACCEPTOR

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior,
                   distance_function=pyabc.distance_functions.NoDistance(),
                   population_size=pop_size,
                   transitions=transition,
                   eps=pyabc.epsilon.NoEpsilon(),
                   acceptor=acceptor,
                   sampler=sampler)

abc.load(db_path)

h = abc.run(minimum_epsilon=0, max_nr_populations=1)

# PLOT

visualize("test3", h)
