"""
This is for testing the stochastic acceptor
with known variance.

Attention:
    Set models.noise = 0.01 to get nicer results.
"""

import pyabc
from models import *
from scipy import stats

import logging
logger = logging.getLogger("Acceptor")
logger.setLevel(logging.DEBUG)

# VARIABLEs

db_path = "sqlite:///db7.db"
n_pops = 8

acceptor = pyabc.StochasticAcceptor()
distance = pyabc.distance_functions.NormalKernel(mean=[0], cov=[noise**2])

# PERFORM ABC ANALYSIS

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior_uvar,
                   distance_function=distance,
                   population_size=100,
                   transitions=transition,
                   eps=pyabc.NoEpsilon(),
                   acceptor=acceptor)

abc.new(db_path, get_y_meas())

h = abc.run(minimum_epsilon=1, max_nr_populations=n_pops)

# PLOT

visualize("test7", h)
